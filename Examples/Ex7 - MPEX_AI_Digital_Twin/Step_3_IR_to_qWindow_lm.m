%% Step_3_IR_to_qWindow_lm_autoFamily.m
% Build q_window_lm from actual flux-tube intersections with the cylindrical window
% and automatically select the source-tube family closest to a reference
% cylindrical crossing z-location.
%
% PURPOSE
%   1) Read IR heat-flux map q^W(theta,z)
%   2) Trace a dense set of seed tubes to the cylindrical window r = R_window
%   3) Sample q^W at the traced hit locations
%   4) Automatically select the connected seed family whose hits are closest
%      to a desired reference z_hit on the cylinder
%   5) Build q_window_lm on flux-tube identity coordinates (seed r, seed theta)
%   6) Build a true theta-z map from selected flux-tube hits
%
% IMPORTANT CONSISTENCY
%   The annulus indices (l,m) are defined by FLUX-TUBE IDENTITY:
%
%       l <- seed radial identity      (hit_seed_r)
%       m <- seed azimuthal identity   (hit_seed_theta)
%
%   The cylindrical window hit location is stored separately as:
%
%       r_window_lm
%       theta_window_deg_lm
%       z_window_lm
%
% INPUT
%   1) CompositeData_ShotSeries_5.mat
%   2) bfield_protoMPEX.nc
%
% OUTPUT
%   window_bundle_qWindow_fromIntersections.mat
%
% -------------------------------------------------------------------------

clear; close all; clc;

%% ------------------------------------------------------------------------
% USER INPUTS
% -------------------------------------------------------------------------
irFile      = 'CompositeData_ShotSeries_1.mat';
profileFile = 'bfield_protoMPEX_shotSeries_1.nc';
outFile     = 'window_bundle_qWindow_fromIntersections.mat';

% IR frame selection
frameMode   = 'max';      % 'single', 'mean', 'max'
singleFrame = 34;

% IR axial coordinate conversion
convert_cm_to_m = false;

% Physical center of experimental window in machine coordinates
z_window_center_machine = 1.75915;   % [m]

% How to align the IR axial coordinate
%   'center_on_window'  -> shift IR z so its midpoint is at z_window_center_machine
%   'use_offset'        -> use a fixed offset z0_machine
zAlignmentMode = 'center_on_window';

% Used only if zAlignmentMode = 'use_offset'
z0_machine = 1.6;         % [m]

% Cylindrical helicon window radius
R_window = 0.0626;        % [m]

% Broad allowed cylindrical window range before automatic family selection
z_window_min_keep = 1.60;   % [m]
z_window_max_keep = 1.9;   % [m]
enforceTargetWindowCoverage = true;
useIRZRangeAsWindowMask = true;

% Automatic family selection target from contour-vs-traced diagnostic
z_ref_target = 1.75915;   % [m]   % no longer used directly

% Primary keep tolerance in cylinder-hit z
dz_keep = 0.015;          % [m]   % no longer used directly

% Optional seed-radius refinement around best family
useSeedRadiusRefinement = false;
dr_seed_keep = 0.006;     % [m]   % no longer used directly

% Final annular indexing
nL = 12;
nM = 36;

% Dense candidate tube seeding
nSeedR = 80;
nSeedTheta = 240;

% Seed plane
z_seed = 1.745;           % [m]
r_seed_min = 0.00;        % [m]
r_seed_max = 0.20;        % [m]

% Tracing controls
ds_trace = 2.0e-4;        % [m]
nStepsMax = 50000;
stop_if_outside_grid = true;
try_reverse_if_forward_fails = true;

% Optional bin cleanup
minTubeCountPerBin = 3;

% Parallel
useParallel = true;

makePlots = true;

%% ------------------------------------------------------------------------
% START PARALLEL POOL
% -------------------------------------------------------------------------
if useParallel
    p = gcp('nocreate');
    if isempty(p)
        parpool;
    end
end

%% ------------------------------------------------------------------------
% STEP 0. LOAD IR DATA q^W(theta,z)
% -------------------------------------------------------------------------
fprintf('Reading IR data: %s ...\n', irFile);

S_ir = load(irFile);
if ~isfield(S_ir, 'f')
    error('Expected structure f in %s.', irFile);
end

f = S_ir.f;
requiredFields = {'dT','phi_2D','s_2D'};
for k = 1:numel(requiredFields)
    if ~isfield(f, requiredFields{k})
        error('Missing field f.%s in IR file.', requiredFields{k});
    end
end

IR     = double(f.dT);
phi_2D = double(f.phi_2D);
s_2D   = double(f.s_2D);

[~, ~, nf] = size(IR);

theta_deg_ir = mod(double(phi_2D(1,:)), 360);
z_like = double(s_2D(:,1));

if convert_cm_to_m
    z_local_m = z_like * 1e-2;
else
    z_local_m = z_like;
end

switch lower(zAlignmentMode)
    case 'center_on_window'
        z_local_center = 0.5 * (min(z_local_m) + max(z_local_m));
        z_ir_m = z_local_m - z_local_center + z_window_center_machine;
    case 'use_offset'
        z_ir_m = z_local_m + z0_machine;
    otherwise
        error('Unknown zAlignmentMode: %s', zAlignmentMode);
end

[theta_deg_ir, idxTheta] = sort(theta_deg_ir);
IR = IR(:, idxTheta, :);

[z_ir_m, idxZ] = sort(z_ir_m);
IR = IR(idxZ, :, :);

switch lower(frameMode)
    case 'single'
        if singleFrame < 1 || singleFrame > nf
            error('singleFrame out of range.');
        end
        q_window_frame = IR(:,:,singleFrame);
    case 'mean'
        q_window_frame = mean(IR, 3, 'omitnan');
    case 'max'
        q_window_frame = max(IR, [], 3, 'omitnan');
    otherwise
        error('Unknown frameMode: %s', frameMode);
end

q_window_frame(~isfinite(q_window_frame)) = NaN;

[theta_deg_ir_unique, ia] = unique(theta_deg_ir, 'stable');
q_window_frame_unique = q_window_frame(:, ia);

[theta_deg_ir_unique, isrt] = sort(theta_deg_ir_unique);
q_window_frame_unique = q_window_frame_unique(:, isrt);

z_window_min = min(z_ir_m);
z_window_max = max(z_ir_m);

fprintf('IR window z-range after alignment = [%g, %g] m\n', z_window_min, z_window_max);
fprintf('IR window center after alignment  = %g m\n', 0.5 * (z_window_min + z_window_max));
fprintf('Broad kept cylindrical z-range    = [%g, %g] m\n', z_window_min_keep, z_window_max_keep);
fprintf('Auto-family reference z target    = %.8f m\n', z_ref_target);

%% ------------------------------------------------------------------------
% STEP 1. LOAD MAGNETIC FIELD
% -------------------------------------------------------------------------
fprintf('Reading magnetic field: %s ...\n', profileFile);

r  = double(ncread(profileFile, 'r')); r = r(:);
z  = double(ncread(profileFile, 'z')); z = z(:);

Br = double(ncread(profileFile, 'br'));
Bz = double(ncread(profileFile, 'bz'));

hasBt = true;
try
    Bt = double(ncread(profileFile, 'bt'));
catch
    hasBt = false;
    Bt = zeros(size(Bz));
end

nR = numel(r);
nZ = numel(z);

if ~isequal(size(Br), [nR, nZ]); Br = permute(Br, [2 1]); end
if ~isequal(size(Bz), [nR, nZ]); Bz = permute(Bz, [2 1]); end
if ~isequal(size(Bt), [nR, nZ]); Bt = permute(Bt, [2 1]); end

fprintf('Magnetic domain: r=[%g,%g], z=[%g,%g]\n', min(r), max(r), min(z), max(z));
fprintf('Toroidal field present in file: %d\n', hasBt);

if z_seed < min(z) || z_seed > max(z)
    error('z_seed is outside magnetic grid.');
end

if R_window < min(r) || R_window > max(r)
    error('R_window is outside magnetic radial grid.');
end

FBr = griddedInterpolant({r, z}, Br, 'linear', 'none');
FBt = griddedInterpolant({r, z}, Bt, 'linear', 'none');
FBz = griddedInterpolant({r, z}, Bz, 'linear', 'none');

FqIR = griddedInterpolant({z_ir_m, theta_deg_ir_unique}, q_window_frame_unique, 'linear', 'none');

%% ------------------------------------------------------------------------
% STEP 2. SEED A DENSE CANDIDATE TUBE BUNDLE
% -------------------------------------------------------------------------
fprintf('Seeding dense candidate tube bundle ...\n');

r_seed_edges = linspace(r_seed_min, r_seed_max, nSeedR+1);
r_seed_cent  = 0.5 * (r_seed_edges(1:end-1) + r_seed_edges(2:end));

theta_seed_edges_deg = linspace(0, 360, nSeedTheta+1);
theta_seed_cent_deg  = 0.5 * (theta_seed_edges_deg(1:end-1) + theta_seed_edges_deg(2:end));

[THETA_SEED_DEG, R_SEED] = meshgrid(theta_seed_cent_deg, r_seed_cent);

seed_r_vec = R_SEED(:);
seed_theta_vec = THETA_SEED_DEG(:);

nSeedTotal = numel(seed_r_vec);
fprintf('Total candidate seed tubes = %d\n', nSeedTotal);

%% ------------------------------------------------------------------------
% STEP 3. TRACE EACH TUBE TO CYLINDRICAL WINDOW r = R_window
% -------------------------------------------------------------------------
fprintf('Tracing seed tubes to cylindrical window r = %.4f m ...\n', R_window);

hit_r_all          = NaN(nSeedTotal,1);
hit_theta_all      = NaN(nSeedTotal,1);
hit_z_all          = NaN(nSeedTotal,1);
hit_qw_all         = NaN(nSeedTotal,1);
hit_seed_r_all     = seed_r_vec;
hit_seed_theta_all = seed_theta_vec;
hit_valid          = false(nSeedTotal,1);

if useParallel
    parfor ii = 1:nSeedTotal
        r0 = seed_r_vec(ii);
        th0 = seed_theta_vec(ii);

        [ok, r_hit, th_hit, z_hit] = trace_to_cyl_window_fast( ...
            r0, th0, z_seed, +1, FBr, FBt, FBz, R_window, ds_trace, nStepsMax, ...
            r, z, stop_if_outside_grid);

        if ~ok && try_reverse_if_forward_fails
            [ok, r_hit, th_hit, z_hit] = trace_to_cyl_window_fast( ...
                r0, th0, z_seed, -1, FBr, FBt, FBz, R_window, ds_trace, nStepsMax, ...
                r, z, stop_if_outside_grid);
        end

        if ~ok
            continue;
        end

        if enforceTargetWindowCoverage
            if z_hit < z_window_min_keep || z_hit > z_window_max_keep
                continue;
            end
        elseif useIRZRangeAsWindowMask
            if z_hit < z_window_min || z_hit > z_window_max
                continue;
            end
        end

        th_mod = mod(th_hit, 360);
        qv = FqIR(z_hit, th_mod);

        if ~isfinite(qv)
            [~, ith] = min(abs(theta_deg_ir_unique - th_mod));
            [~, izh] = min(abs(z_ir_m - z_hit));
            qv = q_window_frame_unique(izh, ith);
        end

        if ~isfinite(qv)
            continue;
        end

        hit_r_all(ii)     = r_hit;
        hit_theta_all(ii) = th_mod;
        hit_z_all(ii)     = z_hit;
        hit_qw_all(ii)    = qv;
        hit_valid(ii)     = true;
    end
else
    for ii = 1:nSeedTotal
        r0 = seed_r_vec(ii);
        th0 = seed_theta_vec(ii);

        [ok, r_hit, th_hit, z_hit] = trace_to_cyl_window_fast( ...
            r0, th0, z_seed, +1, FBr, FBt, FBz, R_window, ds_trace, nStepsMax, ...
            r, z, stop_if_outside_grid);

        if ~ok && try_reverse_if_forward_fails
            [ok, r_hit, th_hit, z_hit] = trace_to_cyl_window_fast( ...
                r0, th0, z_seed, -1, FBr, FBt, FBz, R_window, ds_trace, nStepsMax, ...
                r, z, stop_if_outside_grid);
        end

        if ~ok
            continue;
        end

        if enforceTargetWindowCoverage
            if z_hit < z_window_min_keep || z_hit > z_window_max_keep
                continue;
            end
        elseif useIRZRangeAsWindowMask
            if z_hit < z_window_min || z_hit > z_window_max
                continue;
            end
        end

        th_mod = mod(th_hit, 360);
        qv = FqIR(z_hit, th_mod);

        if ~isfinite(qv)
            [~, ith] = min(abs(theta_deg_ir_unique - th_mod));
            [~, izh] = min(abs(z_ir_m - z_hit));
            qv = q_window_frame_unique(izh, ith);
        end

        if ~isfinite(qv)
            continue;
        end

        hit_r_all(ii)     = r_hit;
        hit_theta_all(ii) = th_mod;
        hit_z_all(ii)     = z_hit;
        hit_qw_all(ii)    = qv;
        hit_valid(ii)     = true;
    end
end

% Keep all broad-valid hits first
hit_r_all_valid          = hit_r_all(hit_valid);
hit_theta_all_valid      = hit_theta_all(hit_valid);
hit_z_all_valid          = hit_z_all(hit_valid);
hit_qw_all_valid         = hit_qw_all(hit_valid);
hit_seed_r_all_valid     = hit_seed_r_all(hit_valid);
hit_seed_theta_all_valid = hit_seed_theta_all(hit_valid);

nHitsBroad = numel(hit_r_all_valid);
fprintf('Broad-valid selected tube/window intersections = %d\n', nHitsBroad);

if nHitsBroad == 0
    error('No valid tube/window intersections found in broad requested z-range.');
end

fprintf('Broad-valid hit z range = [%g, %g] m\n', min(hit_z_all_valid), max(hit_z_all_valid));

%% ------------------------------------------------------------------------
% STEP 3b. KEEP THE FULL FLUX-TUBE BUNDLE THAT HITS THE WINDOW z-RANGE
% -------------------------------------------------------------------------
fprintf('\nSelecting full bundle that hits requested window z-range...\n');

keepMask_family = ...
    (hit_z_all_valid >= z_window_min_keep) & ...
    (hit_z_all_valid <= z_window_max_keep);

fprintf('Hits kept in requested z-window = %d / %d\n', ...
    nnz(keepMask_family), numel(keepMask_family));

hit_r          = hit_r_all_valid(keepMask_family);
hit_theta      = hit_theta_all_valid(keepMask_family);
hit_z          = hit_z_all_valid(keepMask_family);
hit_qw         = hit_qw_all_valid(keepMask_family);
hit_seed_r     = hit_seed_r_all_valid(keepMask_family);
hit_seed_theta = hit_seed_theta_all_valid(keepMask_family);

nHits = numel(hit_r);
if nHits == 0
    error('No hits remain inside requested z-window. Broaden z_window_min_keep / z_window_max_keep.');
end

fprintf('Final selected bundle z range      = [%g, %g] m\n', min(hit_z), max(hit_z));
fprintf('Final selected bundle seed-r range = [%g, %g] m\n', min(hit_seed_r), max(hit_seed_r));
fprintf('Final selected bundle theta range  = [%g, %g] deg\n', min(hit_theta), max(hit_theta));

if nHits < nL*nM
    warning('Fewer selected hits than final annulus bins. Consider increasing nSeedR/nSeedTheta or widening dz_keep.');
end

%% ------------------------------------------------------------------------
% STEP 4. BUILD ANNULUS FROM FLUX-TUBE IDENTITY
% -------------------------------------------------------------------------
fprintf('Building annulus from tube identity ...\n');

r_bundle_min = min(hit_seed_r);
r_bundle_max = max(hit_seed_r);

r_edges = linspace(r_bundle_min, r_bundle_max, nL+1);
r_cent  = 0.5 * (r_edges(1:end-1) + r_edges(2:end));

theta_edges_deg = linspace(0, 360, nM+1);
theta_cent_deg  = 0.5 * (theta_edges_deg(1:end-1) + theta_edges_deg(2:end));

dtheta_deg = theta_edges_deg(2) - theta_edges_deg(1);
dtheta_rad = deg2rad(dtheta_deg);

q_window_lm = NaN(nL, nM);
q_window_std_lm = NaN(nL, nM);
tube_count_lm = zeros(nL, nM);

theta_window_deg_lm = NaN(nL, nM);
z_window_lm         = NaN(nL, nM);
r_window_lm         = NaN(nL, nM);

seed_r_lm           = NaN(nL, nM);
seed_theta_deg_lm   = NaN(nL, nM);

for il = 1:nL
    rmask = (hit_seed_r >= r_edges(il)) & (hit_seed_r < r_edges(il+1));
    if il == nL
        rmask = (hit_seed_r >= r_edges(il)) & (hit_seed_r <= r_edges(il+1));
    end

    for im = 1:nM
        thmask = (hit_seed_theta >= theta_edges_deg(im)) & (hit_seed_theta < theta_edges_deg(im+1));
        if im == nM
            thmask = (hit_seed_theta >= theta_edges_deg(im)) & (hit_seed_theta <= theta_edges_deg(im+1));
        end

        mask = rmask & thmask;

        if any(mask)
            q_window_lm(il,im) = mean(hit_qw(mask), 'omitnan');
            q_window_std_lm(il,im) = std(hit_qw(mask), 'omitnan');
            tube_count_lm(il,im) = nnz(mask);

            theta_window_deg_lm(il,im) = mean(hit_theta(mask), 'omitnan');
            z_window_lm(il,im)         = mean(hit_z(mask), 'omitnan');
            r_window_lm(il,im)         = mean(hit_r(mask), 'omitnan');

            seed_r_lm(il,im)         = mean(hit_seed_r(mask), 'omitnan');
            seed_theta_deg_lm(il,im) = mean(hit_seed_theta(mask), 'omitnan');
        end
    end
end

lowCountMask = tube_count_lm < minTubeCountPerBin;

q_window_lm(lowCountMask) = NaN;
q_window_std_lm(lowCountMask) = NaN;
theta_window_deg_lm(lowCountMask) = NaN;
z_window_lm(lowCountMask) = NaN;
r_window_lm(lowCountMask) = NaN;
seed_r_lm(lowCountMask) = NaN;
seed_theta_deg_lm(lowCountMask) = NaN;

valid = isfinite(q_window_lm);
if any(valid(:)) && any(~valid(:))
    [MM, LL] = meshgrid(1:nM, 1:nL);
    valid_pts = [LL(valid), MM(valid)];
    valid_vals = q_window_lm(valid);

    missing = ~valid;
    missing_pts = [LL(missing), MM(missing)];

    for k = 1:size(missing_pts,1)
        d2 = sum((valid_pts - missing_pts(k,:)).^2, 2);
        [~, idx] = min(d2);
        q_window_lm(missing_pts(k,1), missing_pts(k,2)) = valid_vals(idx);
    end
end

%% ------------------------------------------------------------------------
% STEP 5. BUILD ANNULAR PATCH AREAS
% -------------------------------------------------------------------------
A_window_lm = zeros(nL, nM);
for il = 1:nL
    for im = 1:nM
        A_window_lm(il,im) = 0.5 * (r_edges(il+1)^2 - r_edges(il)^2) * dtheta_rad;
    end
end

%% ------------------------------------------------------------------------
% STEP 6. DIAGNOSTICS
% -------------------------------------------------------------------------
fprintf('Final annulus map finite q_window_lm = %d / %d\n', ...
    sum(isfinite(q_window_lm(:))), numel(q_window_lm));
fprintf('Min q_window_lm = %g\n', min(q_window_lm(:), [], 'omitnan'));
fprintf('Max q_window_lm = %g\n', max(q_window_lm(:), [], 'omitnan'));

fprintf('Bundle radial range used for annulus = [%g, %g] m\n', ...
    min(hit_seed_r), max(hit_seed_r));
fprintf('Cylindrical hit radius range         = [%g, %g] m\n', ...
    min(hit_r), max(hit_r));
fprintf('Binned z_window_lm range             = [%g, %g] m\n', ...
    min(z_window_lm(:), [], 'omitnan'), max(z_window_lm(:), [], 'omitnan'));

%% ------------------------------------------------------------------------
% STEP 7. BUILD TRUE THETA-Z MAP FROM FLUX-TUBE HITS
% -------------------------------------------------------------------------
theta_scatter = hit_theta(:);
z_scatter     = hit_z(:);
q_scatter     = hit_qw(:);

valid_sc = isfinite(theta_scatter) & isfinite(z_scatter) & isfinite(q_scatter);

theta_scatter = theta_scatter(valid_sc);
z_scatter     = z_scatter(valid_sc);
q_scatter     = q_scatter(valid_sc);

if isempty(theta_scatter)
    warning('No valid selected hits available for theta-z interpolation.');
    theta_grid_plot = linspace(0, 360, 360);
    z_grid_plot     = linspace(z_window_min_keep, z_window_max_keep, 240);
    q_window_theta_z_map = NaN(numel(z_grid_plot), numel(theta_grid_plot));
else
    theta_ext = [theta_scatter; theta_scatter-360; theta_scatter+360];
    z_ext     = [z_scatter; z_scatter; z_scatter];
    q_ext     = [q_scatter; q_scatter; q_scatter];

    valid_ext = isfinite(theta_ext) & isfinite(z_ext) & isfinite(q_ext);
    theta_ext = theta_ext(valid_ext);
    z_ext     = z_ext(valid_ext);
    q_ext     = q_ext(valid_ext);

    theta_grid_plot = linspace(0, 360, 360);
    z_grid_plot     = linspace(min(z_scatter), max(z_scatter), 240);
    [THG, ZG] = meshgrid(theta_grid_plot, z_grid_plot);

    if numel(theta_ext) < 3
        warning('Too few valid points for griddata interpolation.');
        q_window_theta_z_map = NaN(size(THG));
    else
        q_window_theta_z_map = griddata(theta_ext, z_ext, q_ext, THG, ZG, 'linear');

        bad = ~isfinite(q_window_theta_z_map);
        if any(bad(:))
            q_window_theta_z_map(bad) = griddata(theta_ext, z_ext, q_ext, ...
                THG(bad), ZG(bad), 'nearest');
        end
    end
end

%% ------------------------------------------------------------------------
% STEP 8. PLOTS
% -------------------------------------------------------------------------
if makePlots
    figure('Color','w');
    imagesc(theta_deg_ir_unique, z_ir_m, q_window_frame_unique);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('machine z [m]');
    title('Measured IR map q^W(\theta,z) on cylindrical window');
    colorbar;

    figure('Color','w','Position',[100 100 1200 420]);

    subplot(1,3,1);
    scatter(hit_theta_all_valid, hit_z_all_valid, 10, hit_qw_all_valid, 'filled'); hold on;
    yline(z_window_min_keep, 'k--', 'LineWidth', 1.2);
    yline(z_window_max_keep, 'k--', 'LineWidth', 1.2);
    yline(z_ref_target, 'r--', 'LineWidth', 1.5);
    xlabel('\theta [deg]');
    ylabel('z_{hit} [m]');
    title('All broad-valid traced hits');
    xlim([0 360]);
    ylim([z_window_min_keep z_window_max_keep]);
    colorbar;
    grid on;
    
    subplot(1,3,2);
    scatter(hit_seed_r_all_valid, hit_z_all_valid, 10, hit_qw_all_valid, 'filled'); hold on;
    yline(z_window_min_keep, 'k--', 'LineWidth', 1.2);
    yline(z_window_max_keep, 'k--', 'LineWidth', 1.2);
    xlabel('seed r [m]');
    ylabel('z_{hit} [m]');
    title('Seed-radius coverage of kept window band');
    colorbar;
    grid on;
    
    subplot(1,3,3);
    scatter(hit_theta, hit_z, 14, hit_qw, 'filled');
    xlabel('\theta [deg]');
    ylabel('z_{hit} [m]');
    title('Final selected bundle');
    xlim([0 360]);
    ylim([z_window_min_keep z_window_max_keep]);
    colorbar;
    grid on;

    figure('Color','w','Position',[100 100 1000 420]);

    subplot(1,2,1);
    scatter(hit_theta, hit_z, 18, hit_qw, 'filled');
    xlabel('\theta [deg]');
    ylabel('window hit z [m]');
    title('Selected flux-tube hits on cylindrical window');
    xlim([0 360]);
    ylim([z_window_min_keep z_window_max_keep]);
    colorbar;

    subplot(1,2,2);
    imagesc(theta_grid_plot, z_grid_plot, q_window_theta_z_map);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('window hit z [m]');
    title('q^W(\theta,z) from flux-tube mapping');
    xlim([0 360]);
    ylim([z_window_min_keep z_window_max_keep]);
    cmax_theta_z = max(q_window_theta_z_map(:), [], 'omitnan');

    if isempty(cmax_theta_z) || ~isscalar(cmax_theta_z) || ~isfinite(cmax_theta_z) || cmax_theta_z <= 0
        cmax_theta_z = 1.0;
    end
    
    clim([0 cmax_theta_z]);
    colorbar;

    figure('Color','w');
    scatter(hit_theta, hit_r, 25, hit_qw, 'filled');
    xlabel('\theta [deg]');
    ylabel('cylindrical window hit radius [m]');
    title('Actual tube/window intersections on cylindrical surface');
    colorbar;

    figure('Color','w');
    scatter(hit_seed_theta, hit_seed_r, 25, hit_qw, 'filled');
    xlabel('seed \theta [deg]');
    ylabel('seed r [m]');
    title('Flux-tube identity used to build annulus');
    colorbar;

    figure('Color','w');
    imagesc(theta_cent_deg, r_cent, q_window_lm);
    set(gca,'YDir','normal');
    xlabel('seed \theta [deg]');
    ylabel('annulus r_{bundle} [m]');
    title('q_{window,lm}^{W} built from tube/window intersections');
    colorbar;

    [THC, RC] = meshgrid(deg2rad(theta_cent_deg), r_cent);
    XX = RC .* cos(THC);
    YY = RC .* sin(THC);

    figure('Color','w');
    scatter(XX(:), YY(:), 60, q_window_lm(:), 'filled');
    axis equal;
    xlabel('x [m]');
    ylabel('y [m]');
    title('Annulus built from flux-tube identity');
    colorbar;

    figure('Color','w','Position',[120 120 1000 420]);

    subplot(1,2,1);
    imagesc(theta_cent_deg, r_cent, tube_count_lm);
    set(gca,'YDir','normal');
    xlabel('seed \theta [deg]');
    ylabel('r_{bundle} [m]');
    title('Number of traced tubes per annulus bin');
    colorbar;

    subplot(1,2,2);
    imagesc(theta_cent_deg, r_cent, z_window_lm);
    set(gca,'YDir','normal');
    xlabel('seed \theta [deg]');
    ylabel('r_{bundle} [m]');
    title('Mapped z_{window,lm}');
    colorbar;
end

%% ------------------------------------------------------------------------
% STEP 9. SAVE
% -------------------------------------------------------------------------
z_window_nominal = 0.5 * (min(z_ir_m) + max(z_ir_m));

save(outFile, ...
    'q_window_lm', ...
    'q_window_std_lm', ...
    'tube_count_lm', ...
    'theta_cent_deg', ...
    'theta_edges_deg', ...
    'r_cent', ...
    'r_edges', ...
    'A_window_lm', ...
    'R_window', ...
    'z_window_nominal', ...
    'z_window_center_machine', ...
    'z_window_min', ...
    'z_window_max', ...
    'z_window_min_keep', ...
    'z_window_max_keep', ...
    'enforceTargetWindowCoverage', ...
    'zAlignmentMode', ...
    'r_window_lm', ...
    'theta_window_deg_lm', ...
    'z_window_lm', ...
    'seed_r_lm', ...
    'seed_theta_deg_lm', ...
    'hit_r', ...
    'hit_theta', ...
    'hit_z', ...
    'hit_qw', ...
    'hit_seed_r', ...
    'hit_seed_theta', ...
    'q_window_theta_z_map', ...
    'theta_grid_plot', ...
    'z_grid_plot', ...
        'minTubeCountPerBin', ...
    'z_ref_target', ...
    'dz_keep', ...
    'useSeedRadiusRefinement', ...
    'dr_seed_keep');

fprintf('Saved %s\n', outFile);

%% ========================================================================
% LOCAL FUNCTION
% ========================================================================
function [ok, r_hit, th_hit, z_hit] = trace_to_cyl_window_fast( ...
    r0, th0_deg, z0, zsign, FBr, FBt, FBz, R_window, ds_trace, nStepsMax, ...
    r_grid, z_grid, stop_if_outside_grid)

    ok = false;
    r_hit = NaN; th_hit = NaN; z_hit = NaN;

    r_now = r0;
    th_now_deg = th0_deg;
    z_now = z0;

    for is = 2:nStepsMax
        br_now = FBr(r_now, z_now);
        bt_now = FBt(r_now, z_now);
        bz_now = FBz(r_now, z_now);

        if any(~isfinite([br_now, bt_now, bz_now]))
            return;
        end

        bmag_now = sqrt(br_now^2 + bt_now^2 + bz_now^2);
        if bmag_now < 1e-12
            return;
        end

        br_u = br_now / bmag_now;
        bt_u = bt_now / bmag_now;
        bz_u = bz_now / bmag_now;

        if sign(bz_u) ~= 0
            sgn = zsign * sign(bz_u);
        else
            sgn = zsign;
        end

        dr_ds = sgn * br_u;
        dtheta_ds = sgn * bt_u / max(r_now, 1e-8);
        dz_ds = sgn * bz_u;

        r_next = r_now + dr_ds * ds_trace;
        th_next_deg = mod(th_now_deg + rad2deg(dtheta_ds * ds_trace), 360);
        z_next = z_now + dz_ds * ds_trace;

        crossed = (r_now - R_window) * (r_next - R_window) <= 0 && abs(r_next - r_now) > 1e-12;
        if crossed
            f = (R_window - r_now) / (r_next - r_now);
            r_hit = R_window;
            th_hit = mod(th_now_deg + f * (th_next_deg - th_now_deg), 360);
            z_hit = z_now + f * (z_next - z_now);
            ok = true;
            return;
        end

        if stop_if_outside_grid
            if r_next < min(r_grid) || r_next > max(r_grid) || ...
               z_next < min(z_grid) || z_next > max(z_grid)
                return;
            end
        end

        r_now = r_next;
        th_now_deg = th_next_deg;
        z_now = z_next;
    end
end