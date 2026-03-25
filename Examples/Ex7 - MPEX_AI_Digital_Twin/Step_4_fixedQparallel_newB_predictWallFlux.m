%% Step_4_fixedQparallel_newB_predictSurfaceFlux.m
% Predict surface heat flux for a new B-field using fixed Q_parallel(l,m),
% and compare against reference window and target maps.
%
% WHAT THIS SCRIPT DOES
%   1) Reads the fixed tube-carried source Q_parallel(l,m)
%   2) Traces each tube in a NEW magnetic field
%   3) Detects the first-hit surface:
%         1 = cylindrical window  (r = R_window)
%         2 = plane surface       (z = z_plane)
%   4) Computes the correct local incidence angle:
%         window: sin(alpha) = |Br| / |B|
%         plane : sin(alpha) = |Bz| / |B|
%   5) Builds predicted new window and plane maps
%   6) Loads REFERENCE window and target maps
%   7) Makes window and target comparison figures
%
% REQUIRED INPUT FILES
%   1) fullSource_Qparallel.mat
%   2) bfield_protoMPEX_shotSeries_5.nc                  (new B-field)
%   3) window_bundle_qWindow_fromIntersections.mat       (reference window)
%   4) trueAnnulusSource_toTarget_fromFullQparallel_output.nc  (reference target)
%
% OUTPUT
%   fixedQparallel_newB_predictSurfaceFlux.nc
%
% -------------------------------------------------------------------------

clear; close all; clc;

%% ------------------------------------------------------------------------
% USER INPUTS
% -------------------------------------------------------------------------
sourceFile = 'fullSource_Qparallel.mat';
newBFile   = 'bfield_protoMPEX_shotSeries_5.nc';
refWindowFile = 'window_bundle_qWindow_fromIntersections.mat';
refTargetFile = 'trueAnnulusSource_toTarget_fromFullQparallel_output.nc';
outFile    = 'fixedQparallel_newB_predictSurfaceFlux.nc';



% Surfaces
R_window = 0.0626;   % [m] cylindrical window radius
z_plane  = 3.900;    % [m] plane surface (target / limiter / dump)

% Tracing
ds_trace = 0.001;    % [m]
nStepsMax = 7000;
stop_if_outside_grid = true;
try_reverse_if_forward_fails = true;

% Rasterization
nThetaGrid = 320;
nZGrid     = 320;
nXGrid     = 320;
nYGrid     = 320;

% Optional window z range for display / clipping
useCustomWindowZRange = false;
z_window_min_custom = 1.60;
z_window_max_custom = 1.90;

% Plot controls
makePlots = true;
makeTrajectoryPlot = true;
nTrajToPlot = 60;

% Common plotting controls for comparison
forceColorAxisFromZero = true;
forceSameWindowAxis = true;
forceSameTargetAxis = true;

%% ------------------------------------------------------------------------
% STEP 0. LOAD REFERENCE WINDOW DATA
% -------------------------------------------------------------------------
fprintf('Reading reference window file: %s ...\n', refWindowFile);

RefW = load(refWindowFile);

reqRefW = {'q_window_lm','theta_cent_deg','r_cent'};
for k = 1:numel(reqRefW)
    if ~isfield(RefW, reqRefW{k})
        error('Missing variable %s in %s.', reqRefW{k}, refWindowFile);
    end
end

q_window_ref_theta_r = double(RefW.q_window_lm);
theta_cent_deg_ref   = double(RefW.theta_cent_deg(:)).';
r_cent_ref           = double(RefW.r_cent(:));

if isfield(RefW,'z_window_lm')
    z_window_lm_ref = double(RefW.z_window_lm);
else
    error('Reference window file must contain z_window_lm for theta-z comparison.');
end

% Reference window theta-z representation
q_window_ref_theta_z = q_window_ref_theta_r;

%% ------------------------------------------------------------------------
% STEP 1. LOAD REFERENCE TARGET DATA
% -------------------------------------------------------------------------
fprintf('Reading reference target file: %s ...\n', refTargetFile);

xt_grid_ref = double(ncread(refTargetFile, 'xt_grid')); xt_grid_ref = xt_grid_ref(:).';
yt_grid_ref = double(ncread(refTargetFile, 'yt_grid')); yt_grid_ref = yt_grid_ref(:).';
q_target_map_ref = double(ncread(refTargetFile, 'q_target_map'));

%% ------------------------------------------------------------------------
% STEP 2. LOAD FIXED Q_parallel(l,m)
% -------------------------------------------------------------------------
fprintf('Reading fixed source file: %s ...\n', sourceFile);

S = load(sourceFile);

requiredVars = {'Q_parallel_full_lm','r_cent_full','r_edges_full', ...
                'theta_cent_deg','theta_edges_deg','A_full_lm'};
for k = 1:numel(requiredVars)
    if ~isfield(S, requiredVars{k})
        error('Missing variable %s in %s.', requiredVars{k}, sourceFile);
    end
end

Q_parallel_lm = double(S.Q_parallel_full_lm);
r_cent = double(S.r_cent_full(:));
r_edges = double(S.r_edges_full(:));
theta_cent_deg = double(S.theta_cent_deg(:)).';
theta_edges_deg = double(S.theta_edges_deg(:)).';
A_source_lm = double(S.A_full_lm);

if isfield(S, 'z_source')
    z_source = double(S.z_source);
else
    z_source = 1.745;
end

[nL, nM] = size(Q_parallel_lm);

if numel(r_cent) ~= nL
    error('Length of r_cent_full does not match rows of Q_parallel_full_lm.');
end
if numel(theta_cent_deg) ~= nM
    error('Length of theta_cent_deg does not match cols of Q_parallel_full_lm.');
end

fprintf('Loaded fixed Q_parallel size = [%d x %d]\n', nL, nM);
fprintf('Source radial range = [%g, %g] m\n', min(r_edges), max(r_edges));
fprintf('z_source = %.6f m\n', z_source);

%% ------------------------------------------------------------------------
% STEP 3. LOAD NEW B-FIELD
% -------------------------------------------------------------------------
fprintf('Reading new B-field file: %s ...\n', newBFile);

r = double(ncread(newBFile, 'r')); r = r(:);
z = double(ncread(newBFile, 'z')); z = z(:);

Br = double(ncread(newBFile, 'br'));
Bt = double(ncread(newBFile, 'bt'));
Bz = double(ncread(newBFile, 'bz'));

nR = numel(r);
nZ = numel(z);

if ~isequal(size(Br), [nR, nZ])
    Br = permute(Br, [2 1]);
    Bt = permute(Bt, [2 1]);
    Bz = permute(Bz, [2 1]);
end

if z_source < min(z) || z_source > max(z)
    error('z_source is outside B-field domain.');
end
if z_plane < min(z) || z_plane > max(z)
    error('z_plane is outside B-field domain.');
end
if R_window < min(r) || R_window > max(r)
    error('R_window is outside B-field radial domain.');
end

FBr = griddedInterpolant({r, z}, Br, 'linear', 'none');
FBt = griddedInterpolant({r, z}, Bt, 'linear', 'none');
FBz = griddedInterpolant({r, z}, Bz, 'linear', 'none');

fprintf('B-field domain: r=[%g,%g] m, z=[%g,%g] m\n', min(r), max(r), min(z), max(z));

%% ------------------------------------------------------------------------
% STEP 4. BUILD FIXED SOURCE ANNULUS PATCHES
% -------------------------------------------------------------------------
fprintf('Building fixed source annulus patches ...\n');

[THETA_LM_DEG, R_LM] = meshgrid(theta_cent_deg, r_cent);
THETA_LM_RAD = deg2rad(THETA_LM_DEG);

Xsrc_lm = R_LM .* cos(THETA_LM_RAD);
Ysrc_lm = R_LM .* sin(THETA_LM_RAD);

xsrc_c = NaN(nL, nM, 4);
ysrc_c = NaN(nL, nM, 4);
zsrc_c = NaN(nL, nM, 4);

for il = 1:nL
    r1 = r_edges(il);
    r2 = r_edges(il+1);

    for im = 1:nM
        th1 = theta_edges_deg(im);
        th2 = theta_edges_deg(im+1);

        rCorners  = [r1, r2, r2, r1];
        thCorners = [th1, th1, th2, th2];

        for ic = 1:4
            xsrc_c(il,im,ic) = rCorners(ic) * cosd(thCorners(ic));
            ysrc_c(il,im,ic) = rCorners(ic) * sind(thCorners(ic));
            zsrc_c(il,im,ic) = z_source;
        end
    end
end

%% ------------------------------------------------------------------------
% STEP 5. TRACE PATCH CENTERS TO FIRST-HIT SURFACE
% -------------------------------------------------------------------------
fprintf('Tracing patch centers to first-hit surface ...\n');

surface_type_lm = zeros(nL, nM);   % 1 = window, 2 = plane

x_path = NaN(nL, nM, nStepsMax);
y_path = NaN(nL, nM, nStepsMax);
z_path = NaN(nL, nM, nStepsMax);

x_hit_lm = NaN(nL, nM);
y_hit_lm = NaN(nL, nM);
z_hit_lm = NaN(nL, nM);
r_hit_lm = NaN(nL, nM);
theta_hit_deg_lm = NaN(nL, nM);

for il = 1:nL
    for im = 1:nM
        r0 = R_LM(il,im);
        th0 = THETA_LM_DEG(il,im);

        [ok, surface_type, xhit, yhit, zhit, xpath, ypath, zpath_] = ...
            trace_point_to_first_surface(r0, th0, z_source, ...
            FBr, FBt, FBz, r, z, R_window, z_plane, ds_trace, ...
            nStepsMax, stop_if_outside_grid, +1);

        if ~ok && try_reverse_if_forward_fails
            [ok, surface_type, xhit, yhit, zhit, xpath, ypath, zpath_] = ...
                trace_point_to_first_surface(r0, th0, z_source, ...
                FBr, FBt, FBz, r, z, R_window, z_plane, ds_trace, ...
                nStepsMax, stop_if_outside_grid, -1);
        end

        if ok
            surface_type_lm(il,im) = surface_type;
            x_hit_lm(il,im) = xhit;
            y_hit_lm(il,im) = yhit;
            z_hit_lm(il,im) = zhit;
            r_hit_lm(il,im) = hypot(xhit, yhit);
            theta_hit_deg_lm(il,im) = mod(atan2d(yhit, xhit), 360);

            ns = min(numel(xpath), nStepsMax);
            x_path(il,im,1:ns) = xpath(1:ns);
            y_path(il,im,1:ns) = ypath(1:ns);
            z_path(il,im,1:ns) = zpath_(1:ns);
        end
    end
end

fprintf('Hit window count = %d\n', nnz(surface_type_lm == 1));
fprintf('Hit plane  count = %d\n', nnz(surface_type_lm == 2));
fprintf('Missed count     = %d\n', nnz(surface_type_lm == 0));

%% ------------------------------------------------------------------------
% STEP 6. TRACE PATCH CORNERS TO FIRST-HIT SURFACE
% -------------------------------------------------------------------------
fprintf('Tracing patch corners to first-hit surface ...\n');

xhit_c = NaN(nL, nM, 4);
yhit_c = NaN(nL, nM, 4);
zhit_c = NaN(nL, nM, 4);
surface_type_patch = zeros(nL, nM);
patch_ok = false(nL, nM);
fail_corner = NaN(nL, nM);

for il = 1:nL
    for im = 1:nM
        st_center = surface_type_lm(il,im);

        if st_center == 0
            patch_ok(il,im) = false;
            continue;
        end

        ok_all = true;

        for ic = 1:4
            rc  = hypot(xsrc_c(il,im,ic), ysrc_c(il,im,ic));
            thc = mod(atan2d(ysrc_c(il,im,ic), xsrc_c(il,im,ic)), 360);
            zc  = zsrc_c(il,im,ic);

            [ok, surface_type, xhit, yhit, zhit] = ...
                trace_point_to_first_surface(rc, thc, zc, ...
                FBr, FBt, FBz, r, z, R_window, z_plane, ds_trace, ...
                nStepsMax, stop_if_outside_grid, +1);

            if ~ok && try_reverse_if_forward_fails
                [ok, surface_type, xhit, yhit, zhit] = ...
                    trace_point_to_first_surface(rc, thc, zc, ...
                    FBr, FBt, FBz, r, z, R_window, z_plane, ds_trace, ...
                    nStepsMax, stop_if_outside_grid, -1);
            end

            if ~ok || surface_type ~= st_center
                ok_all = false;
                fail_corner(il,im) = ic;
                break;
            end

            xhit_c(il,im,ic) = xhit;
            yhit_c(il,im,ic) = yhit;
            zhit_c(il,im,ic) = zhit;
        end

        patch_ok(il,im) = ok_all;
        if ok_all
            surface_type_patch(il,im) = st_center;
        end
    end
end

fprintf('Valid window patches = %d\n', nnz(patch_ok & surface_type_patch == 1));
fprintf('Valid plane patches  = %d\n', nnz(patch_ok & surface_type_patch == 2));

%% ------------------------------------------------------------------------
% STEP 7. COMPUTE SURFACE HEAT FLUX USING LOCAL TUBE ORIENTATION
% -------------------------------------------------------------------------
fprintf('Computing surface incidence and heat flux ...\n');

sin_alpha_hit_lm = NaN(nL, nM);
q_surface_lm = NaN(nL, nM);

for il = 1:nL
    for im = 1:nM
        if surface_type_lm(il,im) == 0
            continue;
        end

        rr = r_hit_lm(il,im);
        zz = z_hit_lm(il,im);

        br_h = FBr(rr, zz);
        bt_h = FBt(rr, zz);
        bz_h = FBz(rr, zz);

        if any(~isfinite([br_h, bt_h, bz_h]))
            continue;
        end

        Bmag_h = sqrt(br_h^2 + bt_h^2 + bz_h^2);
        if Bmag_h < 1e-12
            continue;
        end

        switch surface_type_lm(il,im)
            case 1
                sin_alpha_hit_lm(il,im) = abs(br_h) / Bmag_h;
            case 2
                sin_alpha_hit_lm(il,im) = abs(bz_h) / Bmag_h;
        end

        q_surface_lm(il,im) = Q_parallel_lm(il,im) .* sin_alpha_hit_lm(il,im);
    end
end

%% ------------------------------------------------------------------------
% STEP 8. SPLIT WINDOW AND PLANE OUTPUTS
% -------------------------------------------------------------------------
q_window_lm_new = NaN(nL, nM);
q_plane_lm_new  = NaN(nL, nM);

q_window_lm_new(surface_type_lm == 1) = q_surface_lm(surface_type_lm == 1);
q_plane_lm_new(surface_type_lm  == 2) = q_surface_lm(surface_type_lm  == 2);

%% ------------------------------------------------------------------------
% STEP 9. RASTERIZE WINDOW CYLINDRICAL MAP IN (theta,z)
% -------------------------------------------------------------------------
fprintf('Rasterizing window cylindrical map ...\n');

windowMask = (patch_ok & surface_type_patch == 1);

if any(windowMask(:))
    theta_window_vals = theta_hit_deg_lm(windowMask);
    z_window_vals = z_hit_lm(windowMask);

    if useCustomWindowZRange
        zmin_w = z_window_min_custom;
        zmax_w = z_window_max_custom;
    else
        zmin_w = min(z_window_vals, [], 'omitnan');
        zmax_w = max(z_window_vals, [], 'omitnan');
    end

    theta_grid = linspace(0, 360, nThetaGrid);
    z_grid = linspace(zmin_w, zmax_w, nZGrid);
    [THG, ZG] = meshgrid(theta_grid, z_grid);

    q_window_sum   = zeros(nZGrid, nThetaGrid);
    q_window_count = zeros(nZGrid, nThetaGrid);

    for il = 1:nL
        for im = 1:nM
            if ~windowMask(il,im) || ~isfinite(q_window_lm_new(il,im))
                continue;
            end

            thv = squeeze(mod(atan2d(yhit_c(il,im,:), xhit_c(il,im,:)), 360));
            zv  = squeeze(zhit_c(il,im,:));
            thc = mod(theta_hit_deg_lm(il,im), 360);

            thv_unwrap = thv;
            thv_unwrap(thv - thc > 180)  = thv_unwrap(thv - thc > 180) - 360;
            thv_unwrap(thv - thc < -180) = thv_unwrap(thv - thc < -180) + 360;

            THG_local = THG;
            THG_local(THG - thc > 180)  = THG_local(THG - thc > 180) - 360;
            THG_local(THG - thc < -180) = THG_local(THG - thc < -180) + 360;

            in = inpolygon(THG_local, ZG, thv_unwrap, zv);

            q_window_sum(in)   = q_window_sum(in) + q_window_lm_new(il,im);
            q_window_count(in) = q_window_count(in) + 1;
        end
    end

    q_window_map_new = NaN(nZGrid, nThetaGrid);
    maskW = q_window_count > 0;
    q_window_map_new(maskW) = q_window_sum(maskW) ./ q_window_count(maskW);
else
    theta_grid = linspace(0, 360, nThetaGrid);
    z_grid = linspace(z_source-0.1, z_source+0.1, nZGrid);
    q_window_map_new = NaN(nZGrid, nThetaGrid);
end

%% ------------------------------------------------------------------------
% STEP 10. RASTERIZE PLANE MAP IN (x,y)
% -------------------------------------------------------------------------
fprintf('Rasterizing plane map ...\n');

planeMask = (patch_ok & surface_type_patch == 2);

if any(planeMask(:))
    xmin = min(xhit_c(planeMask), [], 'omitnan');
    xmax = max(xhit_c(planeMask), [], 'omitnan');
    ymin = min(yhit_c(planeMask), [], 'omitnan');
    ymax = max(yhit_c(planeMask), [], 'omitnan');

    pad = 0.002;
    x_grid = linspace(xmin-pad, xmax+pad, nXGrid);
    y_grid = linspace(ymin-pad, ymax+pad, nYGrid);
    [XG, YG] = meshgrid(x_grid, y_grid);

    q_plane_sum   = zeros(nYGrid, nXGrid);
    q_plane_count = zeros(nYGrid, nXGrid);

    for il = 1:nL
        for im = 1:nM
            if ~planeMask(il,im) || ~isfinite(q_plane_lm_new(il,im))
                continue;
            end

            xv = squeeze(xhit_c(il,im,:));
            yv = squeeze(yhit_c(il,im,:));

            in = inpolygon(XG, YG, xv, yv);

            q_plane_sum(in)   = q_plane_sum(in) + q_plane_lm_new(il,im);
            q_plane_count(in) = q_plane_count(in) + 1;
        end
    end

    q_plane_map_new = NaN(nYGrid, nXGrid);
    maskP = q_plane_count > 0;
    q_plane_map_new(maskP) = q_plane_sum(maskP) ./ q_plane_count(maskP);
else
    x_grid = linspace(-0.05, 0.05, nXGrid);
    y_grid = linspace(-0.05, 0.05, nYGrid);
    q_plane_map_new = NaN(nYGrid, nXGrid);
end

% Aliases for comparison block
theta_grid_new = theta_grid;
z_grid_new = z_grid;
xt_grid_new = x_grid;
yt_grid_new = y_grid;
q_target_map_new = q_plane_map_new;

%% ------------------------------------------------------------------------
% STEP 11. DIAGNOSTICS
% -------------------------------------------------------------------------
fprintf('Finite q_window_lm_new = %d / %d\n', sum(isfinite(q_window_lm_new(:))), numel(q_window_lm_new));
fprintf('Finite q_plane_lm_new  = %d / %d\n', sum(isfinite(q_plane_lm_new(:))), numel(q_plane_lm_new));

if any(isfinite(q_window_lm_new(:)))
    fprintf('Window: min/max/mean sin(alpha) = %g / %g / %g\n', ...
        min(sin_alpha_hit_lm(surface_type_lm==1), [], 'omitnan'), ...
        max(sin_alpha_hit_lm(surface_type_lm==1), [], 'omitnan'), ...
        mean(sin_alpha_hit_lm(surface_type_lm==1), 'omitnan'));
end

if any(isfinite(q_plane_lm_new(:)))
    fprintf('Plane:  min/max/mean sin(alpha) = %g / %g / %g\n', ...
        min(sin_alpha_hit_lm(surface_type_lm==2), [], 'omitnan'), ...
        max(sin_alpha_hit_lm(surface_type_lm==2), [], 'omitnan'), ...
        mean(sin_alpha_hit_lm(surface_type_lm==2), 'omitnan'));
end

%% ------------------------------------------------------------------------
% STEP 12. PLOTS OF NEW PREDICTION
% -------------------------------------------------------------------------
if makePlots
    figure('Color','w','Position',[100 100 1150 420]);

    subplot(1,2,1);
    imagesc(theta_cent_deg, r_cent, Q_parallel_lm);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('source annulus radius [m]');
    title('Fixed Q_{\parallel,lm}');
    colorbar;

    subplot(1,2,2);
    imagesc(theta_cent_deg, r_cent, surface_type_lm);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('source annulus radius [m]');
    title('First-hit surface type (1=window, 2=plane)');
    colorbar;

    figure('Color','w');
    imagesc(theta_cent_deg, r_cent, sin_alpha_hit_lm);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('source annulus radius [m]');
    title('Incidence factor from local tube orientation');
    colorbar;

    figure('Color','w','Position',[120 120 1150 420]);

    subplot(1,2,1);
    imagesc(theta_cent_deg, r_cent, q_window_lm_new);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('source annulus radius [m]');
    title('Predicted new window heat flux');
    colorbar;

    subplot(1,2,2);
    imagesc(theta_cent_deg, r_cent, q_plane_lm_new);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('source annulus radius [m]');
    title('Predicted new plane heat flux');
    colorbar;

    if any(windowMask(:))
        figure('Color','w');
        imagesc(theta_grid, z_grid, q_window_map_new);
        set(gca,'YDir','normal');
        xlabel('\theta [deg]');
        ylabel('z [m]');
        title('Predicted window heat-flux map');
        colorbar;
    end

    if any(planeMask(:))
        figure('Color','w');
        imagesc(x_grid, y_grid, q_plane_map_new);
        set(gca,'YDir','normal');
        axis equal tight;
        xlabel('x [m]');
        ylabel('y [m]');
        title('Predicted plane heat-flux map');
        colorbar;
    end

    figure('Color','w');
    imagesc(theta_cent_deg, r_cent, patch_ok);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('source annulus radius [m]');
    title('patch\_ok map');
    colorbar;

    if makeTrajectoryPlot
        figure('Color','w');
        hold on;
        count = 0;
        for il = 1:nL
            for im = 1:nM
                xx = squeeze(x_path(il,im,:));
                yy = squeeze(y_path(il,im,:));
                zz = squeeze(z_path(il,im,:));
                good = isfinite(xx) & isfinite(yy) & isfinite(zz);
                if any(good)
                    plot3(xx(good), yy(good), zz(good), 'LineWidth', 1.0);
                    count = count + 1;
                end
                if count >= nTrajToPlot
                    break;
                end
            end
            if count >= nTrajToPlot
                break;
            end
        end
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        title('Source-to-first-surface trajectories');
        grid on;
        axis equal;
        view(3);
    end
end

%% ------------------------------------------------------------------------
% STEP 13. COMPARISON BLOCK: WINDOW
% Build reference window map from scattered (\theta,z) points directly
% -------------------------------------------------------------------------

% Make sure reference arrays are compatible
if ~isequal(size(q_window_ref_theta_r), size(z_window_lm_ref))
    if isequal(size(q_window_ref_theta_r.'), size(z_window_lm_ref))
        q_window_ref_theta_r = q_window_ref_theta_r.';
    elseif isequal(size(z_window_lm_ref.'), size(q_window_ref_theta_r))
        z_window_lm_ref = z_window_lm_ref.';
    else
        error('q_window_ref_theta_r and z_window_lm_ref are incompatible in size.');
    end
end

[nL_ref, nM_ref] = size(q_window_ref_theta_r);

if numel(theta_cent_deg_ref) ~= nM_ref
    if numel(theta_cent_deg_ref) == nL_ref && nL_ref == nM_ref
        warning('theta_cent_deg_ref length is ambiguous; using transpose assumption.');
        q_window_ref_theta_r = q_window_ref_theta_r.';
        z_window_lm_ref = z_window_lm_ref.';
        [nL_ref, nM_ref] = size(q_window_ref_theta_r);
    else
        error('theta_cent_deg_ref length does not match reference theta dimension.');
    end
end

% Build scattered reference coordinates
TH_REF_LM = repmat(theta_cent_deg_ref(:).', nL_ref, 1);
Z_REF_LM  = z_window_lm_ref;
Q_REF_LM  = q_window_ref_theta_r;

% Extend theta periodically to handle the 0/360 seam
th_ref_valid = [];
z_ref_valid  = [];
q_ref_valid  = [];

for shift = [-360, 0, 360]
    TH_tmp = TH_REF_LM + shift;
    valid_tmp = isfinite(TH_tmp) & isfinite(Z_REF_LM) & isfinite(Q_REF_LM);

    th_ref_valid = [th_ref_valid; TH_tmp(valid_tmp)]; %#ok<AGROW>
    z_ref_valid  = [z_ref_valid;  Z_REF_LM(valid_tmp)]; %#ok<AGROW>
    q_ref_valid  = [q_ref_valid;  Q_REF_LM(valid_tmp)]; %#ok<AGROW>
end

if isempty(th_ref_valid)
    error('No finite reference window points available for theta-z interpolation.');
end

[TH_NEW, Z_NEW] = meshgrid(theta_grid_new, z_grid_new);

q_window_ref_on_newgrid = griddata(th_ref_valid, z_ref_valid, q_ref_valid, ...
    TH_NEW, Z_NEW, 'linear');

bad = ~isfinite(q_window_ref_on_newgrid);
if any(bad(:))
    q_window_ref_on_newgrid(bad) = griddata(th_ref_valid, z_ref_valid, q_ref_valid, ...
        TH_NEW(bad), Z_NEW(bad), 'nearest');
end

% For the separate "reference q^W shown in (\theta,z)" plot,
% create a regular display grid from the same scattered data
z_ref_axis_plot = linspace(min(z_ref_valid), max(z_ref_valid), 200);
theta_ref_axis_plot = linspace(0, 360, 361);
[TH_REF_PLOT, Z_REF_PLOT] = meshgrid(theta_ref_axis_plot, z_ref_axis_plot);

q_window_ref_theta_z_plot = griddata(th_ref_valid, z_ref_valid, q_ref_valid, ...
    TH_REF_PLOT, Z_REF_PLOT, 'linear');

bad_plot = ~isfinite(q_window_ref_theta_z_plot);
if any(bad_plot(:))
    q_window_ref_theta_z_plot(bad_plot) = griddata(th_ref_valid, z_ref_valid, q_ref_valid, ...
        TH_REF_PLOT(bad_plot), Z_REF_PLOT(bad_plot), 'nearest');
end

dq_window_theta_z = q_window_map_new - q_window_ref_on_newgrid;

ratio_window_theta_z = NaN(size(q_window_map_new));
mask_ref_w = abs(q_window_ref_on_newgrid) > 1e-12;
ratio_window_theta_z(mask_ref_w) = q_window_map_new(mask_ref_w) ./ q_window_ref_on_newgrid(mask_ref_w);

%% ------------------------------------------------------------------------
% STEP 14. COMPARISON BLOCK: TARGET
% -------------------------------------------------------------------------
[X_REF, Y_REF] = meshgrid(xt_grid_ref, yt_grid_ref);
[X_NEW, Y_NEW] = meshgrid(xt_grid_new, yt_grid_new);

q_target_ref_on_newgrid = interp2(X_REF, Y_REF, q_target_map_ref, X_NEW, Y_NEW, 'linear', NaN);
bad = ~isfinite(q_target_ref_on_newgrid);
if any(bad(:))
    q_target_ref_on_newgrid(bad) = interp2(X_REF, Y_REF, q_target_map_ref, X_NEW(bad), Y_NEW(bad), 'nearest');
end

dq_target = q_target_map_new - q_target_ref_on_newgrid;

ratio_target = NaN(size(q_target_map_new));
mask_ref_t = abs(q_target_ref_on_newgrid) > 1e-12;
ratio_target(mask_ref_t) = q_target_map_new(mask_ref_t) ./ q_target_ref_on_newgrid(mask_ref_t);

%% ------------------------------------------------------------------------
% STEP 14b. COMMON AXES AND COLOR LIMITS FOR COMPARISON PLOTS
% -------------------------------------------------------------------------

% ---------- Window axes ----------
% ---------- Window axes ----------
theta_min_plot = 0;
theta_max_plot = 360;

if forceSameWindowAxis
    zmin_plot_window = min([z_ref_axis_plot(:); z_grid_new(:)], [], 'omitnan');
    zmax_plot_window = max([z_ref_axis_plot(:); z_grid_new(:)], [], 'omitnan');
else
    zmin_plot_window = min(z_grid_new(:), [], 'omitnan');
    zmax_plot_window = max(z_grid_new(:), [], 'omitnan');
end

if ~isfinite(zmin_plot_window) || ~isfinite(zmax_plot_window) || zmax_plot_window <= zmin_plot_window
    zmin_plot_window = min(z_grid_new(:), [], 'omitnan');
    zmax_plot_window = max(z_grid_new(:), [], 'omitnan');
end

% ---------- Target axes ----------
if forceSameTargetAxis
    xmin_plot_target = min([xt_grid_ref(:); xt_grid_new(:)], [], 'omitnan');
    xmax_plot_target = max([xt_grid_ref(:); xt_grid_new(:)], [], 'omitnan');
    ymin_plot_target = min([yt_grid_ref(:); yt_grid_new(:)], [], 'omitnan');
    ymax_plot_target = max([yt_grid_ref(:); yt_grid_new(:)], [], 'omitnan');
else
    xmin_plot_target = min(xt_grid_new, [], 'omitnan');
    xmax_plot_target = max(xt_grid_new, [], 'omitnan');
    ymin_plot_target = min(yt_grid_new, [], 'omitnan');
    ymax_plot_target = max(yt_grid_new, [], 'omitnan');
end

% ---------- Common color limits ----------
if forceColorAxisFromZero
    wmin = 0;
    tmin = 0;
else
    wmin = min([q_window_ref_on_newgrid(:); q_window_map_new(:)], [], 'omitnan');
    tmin = min([q_target_ref_on_newgrid(:); q_target_map_new(:)], [], 'omitnan');
end

wmax = max([q_window_ref_on_newgrid(:); q_window_map_new(:)], [], 'omitnan');
tmax = max([q_target_ref_on_newgrid(:); q_target_map_new(:)], [], 'omitnan');

% Safety in case maps are empty
if ~isfinite(wmax) || wmax <= wmin
    wmax = wmin + 1;
end
if ~isfinite(tmax) || tmax <= tmin
    tmax = tmin + 1;
end

%% ------------------------------------------------------------------------
% STEP 15. COMPARISON FIGURES
% -------------------------------------------------------------------------
if makePlots

    % Window comparison
    figure('Color','w','Position',[80 80 1250 900]);

    subplot(2,2,1);
    imagesc(theta_grid_new, z_grid_new, q_window_ref_on_newgrid);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('window hit z [m]');
    title('Reference window heat flux in (\theta,z)');
    xlim([theta_min_plot theta_max_plot]);
    ylim([zmin_plot_window zmax_plot_window]);
    clim([wmin wmax]);
    colorbar;

    subplot(2,2,2);
    imagesc(theta_grid_new, z_grid_new, q_window_map_new);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('window hit z [m]');
    title('Predicted new window heat flux in (\theta,z)');
    xlim([theta_min_plot theta_max_plot]);
    ylim([zmin_plot_window zmax_plot_window]);
    clim([wmin wmax]);
    colorbar;

    subplot(2,2,3);
    imagesc(theta_grid_new, z_grid_new, dq_window_theta_z);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('window hit z [m]');
    title('\Delta q_{window} = q_{new} - q_{ref}');
    xlim([theta_min_plot theta_max_plot]);
    ylim([zmin_plot_window zmax_plot_window]);
    colorbar;

    subplot(2,2,4);
    imagesc(theta_grid_new, z_grid_new, ratio_window_theta_z);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('window hit z [m]');
    title('q_{new}/q_{ref} at window');
    xlim([theta_min_plot theta_max_plot]);
    ylim([zmin_plot_window zmax_plot_window]);
    colorbar;

    sgtitle('Window heat-flux comparison');

    % Keep original reference style
    figure('Color','w','Position',[100 100 1150 420]);

    subplot(1,2,1);
    imagesc(theta_ref_axis_plot, z_ref_axis_plot, q_window_ref_theta_z_plot);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('window hit z [m]');
    title('Reference q^W shown in (\theta,z)');
    xlim([theta_min_plot theta_max_plot]);
    ylim([zmin_plot_window zmax_plot_window]);
    clim([wmin wmax]);
    colorbar;

    subplot(1,2,2);
    imagesc(theta_cent_deg_ref, r_cent_ref, q_window_ref_theta_r);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('source annulus r [m]');
    title('Reference q^W shown in (\theta,r)');
    xlim([theta_min_plot theta_max_plot]);
    clim([wmin wmax]);
    colorbar;

    % Matching predicted style
    figure('Color','w','Position',[120 120 1150 420]);

    subplot(1,2,1);
    imagesc(theta_grid_new, z_grid_new, q_window_map_new);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('window hit z [m]');
    title('Predicted q^W shown in (\theta,z)');
    xlim([theta_min_plot theta_max_plot]);
    ylim([zmin_plot_window zmax_plot_window]);
    clim([wmin wmax]);
    colorbar;

    subplot(1,2,2);
    imagesc(theta_cent_deg, r_cent, q_window_lm_new);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('source annulus r [m]');
    title('Predicted q^W shown in (\theta,r)');
    xlim([theta_min_plot theta_max_plot]);
    clim([wmin wmax]);
    colorbar;

    % Target comparison
    figure('Color','w','Position',[100 100 1250 900]);

    subplot(2,2,1);
    imagesc(xt_grid_new, yt_grid_new, q_target_ref_on_newgrid);
    set(gca,'YDir','normal');
    axis equal;
    xlim([xmin_plot_target xmax_plot_target]);
    ylim([ymin_plot_target ymax_plot_target]);
    xlabel('x target [m]');
    ylabel('y target [m]');
    title('Reference target heat-flux map');
    clim([tmin tmax]);
    colorbar;

    subplot(2,2,2);
    imagesc(xt_grid_new, yt_grid_new, q_target_map_new);
    set(gca,'YDir','normal');
    axis equal;
    xlim([xmin_plot_target xmax_plot_target]);
    ylim([ymin_plot_target ymax_plot_target]);
    xlabel('x target [m]');
    ylabel('y target [m]');
    title('Predicted new target heat-flux map');
    clim([tmin tmax]);
    colorbar;

    subplot(2,2,3);
    imagesc(xt_grid_new, yt_grid_new, dq_target);
    set(gca,'YDir','normal');
    axis equal;
    xlim([xmin_plot_target xmax_plot_target]);
    ylim([ymin_plot_target ymax_plot_target]);
    xlabel('x target [m]');
    ylabel('y target [m]');
    title('\Delta q_{target} = q_{new} - q_{ref}');
    colorbar;

    subplot(2,2,4);
    imagesc(xt_grid_new, yt_grid_new, ratio_target);
    set(gca,'YDir','normal');
    axis equal;
    xlim([xmin_plot_target xmax_plot_target]);
    ylim([ymin_plot_target ymax_plot_target]);
    xlabel('x target [m]');
    ylabel('y target [m]');
    title('q_{new}/q_{ref} at target');
    colorbar;

    sgtitle('Target heat-flux comparison');

    fprintf('\nWINDOW COMPARISON:\n');
    fprintf('Reference window max = %.6e\n', max(q_window_ref_on_newgrid(:), [], 'omitnan'));
    fprintf('Predicted window max = %.6e\n', max(q_window_map_new(:), [], 'omitnan'));
    fprintf('Window diff max      = %.6e\n', max(dq_window_theta_z(:), [], 'omitnan'));
    fprintf('Window diff min      = %.6e\n', min(dq_window_theta_z(:), [], 'omitnan'));

    fprintf('\nTARGET COMPARISON:\n');
    fprintf('Reference target max = %.6e\n', max(q_target_ref_on_newgrid(:), [], 'omitnan'));
    fprintf('Predicted target max = %.6e\n', max(q_target_map_new(:), [], 'omitnan'));
    fprintf('Target diff max      = %.6e\n', max(dq_target(:), [], 'omitnan'));
    fprintf('Target diff min      = %.6e\n', min(dq_target(:), [], 'omitnan'));
end

%% ------------------------------------------------------------------------
% STEP 16. SAVE OUTPUT
% -------------------------------------------------------------------------
fprintf('Writing %s ...\n', outFile);

if exist(outFile, 'file')
    delete(outFile);
end

nccreate(outFile, 'theta_cent_deg', 'Dimensions', {'m', nM});
nccreate(outFile, 'r_cent', 'Dimensions', {'l', nL});

nccreate(outFile, 'Q_parallel_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'surface_type_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'sin_alpha_hit_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'q_window_lm_new', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'q_plane_lm_new', 'Dimensions', {'l', nL, 'm', nM});

nccreate(outFile, 'r_hit_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'z_hit_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'theta_hit_deg_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'x_hit_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'y_hit_lm', 'Dimensions', {'l', nL, 'm', nM});

nccreate(outFile, 'q_window_map_new', 'Dimensions', {'z_w', numel(z_grid), 'theta_w', numel(theta_grid)});
nccreate(outFile, 'theta_grid', 'Dimensions', {'theta_w', numel(theta_grid)});
nccreate(outFile, 'z_grid', 'Dimensions', {'z_w', numel(z_grid)});

nccreate(outFile, 'q_plane_map_new', 'Dimensions', {'y_p', numel(y_grid), 'x_p', numel(x_grid)});
nccreate(outFile, 'x_grid', 'Dimensions', {'x_p', numel(x_grid)});
nccreate(outFile, 'y_grid', 'Dimensions', {'y_p', numel(y_grid)});

ncwrite(outFile, 'theta_cent_deg', theta_cent_deg);
ncwrite(outFile, 'r_cent', r_cent);

ncwrite(outFile, 'Q_parallel_lm', Q_parallel_lm);
ncwrite(outFile, 'surface_type_lm', surface_type_lm);
ncwrite(outFile, 'sin_alpha_hit_lm', sin_alpha_hit_lm);
ncwrite(outFile, 'q_window_lm_new', q_window_lm_new);
ncwrite(outFile, 'q_plane_lm_new', q_plane_lm_new);

ncwrite(outFile, 'r_hit_lm', r_hit_lm);
ncwrite(outFile, 'z_hit_lm', z_hit_lm);
ncwrite(outFile, 'theta_hit_deg_lm', theta_hit_deg_lm);
ncwrite(outFile, 'x_hit_lm', x_hit_lm);
ncwrite(outFile, 'y_hit_lm', y_hit_lm);

ncwrite(outFile, 'q_window_map_new', q_window_map_new);
ncwrite(outFile, 'theta_grid', theta_grid);
ncwrite(outFile, 'z_grid', z_grid);

ncwrite(outFile, 'q_plane_map_new', q_plane_map_new);
ncwrite(outFile, 'x_grid', x_grid);
ncwrite(outFile, 'y_grid', y_grid);

fprintf('Done.\n');

%% ========================================================================
% LOCAL FUNCTION
% ========================================================================
function [ok, surface_type, xhit, yhit, zhit, xpath, ypath, zpath_] = ...
    trace_point_to_first_surface(r0, th0_deg, z0, ...
    FBr, FBt, FBz, r_grid, z_grid, R_window, z_plane, ds_trace, ...
    nStepsMax, stop_if_outside_grid, trace_sign)

    ok = false;
    surface_type = 0;
    xhit = NaN; yhit = NaN; zhit = NaN;

    xpath = NaN(1, nStepsMax);
    ypath = NaN(1, nStepsMax);
    zpath_ = NaN(1, nStepsMax);

    r_now = r0;
    th_now_deg = th0_deg;
    z_now = z0;

    xpath(1) = r_now * cosd(th_now_deg);
    ypath(1) = r_now * sind(th_now_deg);
    zpath_(1) = z_now;

    for is = 2:nStepsMax
        br_now = FBr(r_now, z_now);
        bt_now = FBt(r_now, z_now);
        bz_now = FBz(r_now, z_now);

        if any(~isfinite([br_now, bt_now, bz_now]))
            return;
        end

        Bmag_now = sqrt(br_now^2 + bt_now^2 + bz_now^2);
        if Bmag_now < 1e-12
            return;
        end

        br_u = br_now / Bmag_now;
        bt_u = bt_now / Bmag_now;
        bz_u = bz_now / Bmag_now;

        dr_ds = trace_sign * br_u;
        dtheta_ds = trace_sign * bt_u / max(r_now, 1e-8);
        dz_ds = trace_sign * bz_u;

        r_next = r_now + dr_ds * ds_trace;
        th_next_deg = mod(th_now_deg + rad2deg(dtheta_ds * ds_trace), 360);
        z_next = z_now + dz_ds * ds_trace;

        has_window = ((r_now - R_window) * (r_next - R_window) <= 0) && abs(r_next - r_now) > 1e-12;
        has_plane  = ((z_now - z_plane)  * (z_next - z_plane)  <= 0) && abs(z_next - z_now) > 1e-12;

        f_window = inf;
        f_plane  = inf;

        if has_window
            f_window = (R_window - r_now) / (r_next - r_now);
            if f_window < 0 || f_window > 1
                f_window = inf;
            end
        end

        if has_plane
            f_plane = (z_plane - z_now) / (z_next - z_now);
            if f_plane < 0 || f_plane > 1
                f_plane = inf;
            end
        end

        if isfinite(f_window) || isfinite(f_plane)
            if f_window < f_plane
                f = f_window;
                surface_type = 1;
                r_hit = R_window;
                th_hit = mod(th_now_deg + f * (th_next_deg - th_now_deg), 360);
                z_hit = z_now + f * (z_next - z_now);
            else
                f = f_plane;
                surface_type = 2;
                r_hit = r_now + f * (r_next - r_now);
                th_hit = mod(th_now_deg + f * (th_next_deg - th_now_deg), 360);
                z_hit = z_plane;
            end

            xhit = r_hit * cosd(th_hit);
            yhit = r_hit * sind(th_hit);
            zhit = z_hit;

            xpath(is) = xhit;
            ypath(is) = yhit;
            zpath_(is) = zhit;
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

        xpath(is) = r_now * cosd(th_now_deg);
        ypath(is) = r_now * sind(th_now_deg);
        zpath_(is) = z_now;
    end
end