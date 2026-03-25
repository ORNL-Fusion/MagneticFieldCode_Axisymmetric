%% compare_contour_touch_vs_traced_hit.m
% Compare:
%   1) contour-based window touch point from constant-flux geometry
%   2) traced cylindrical crossing point from B-field integration
%
% PURPOSE
%   Reproduce the logic of the older flux-mapping script and compare it
%   directly with the tracing approach used in the newer scripts.
%
% ASSUMPTIONS
%   - Axisymmetric field
%   - NetCDF file contains r, z, br, bz
%   - bt may be absent; if absent, it is set to zero
%
% INPUT
%   bfield_protoMPEX.nc
%
% OUTPUT
%   Diagnostic plots and printed comparison numbers
%
% -------------------------------------------------------------------------

clear; close all; clc;

%% ------------------------------------------------------------------------
% USER INPUTS
% -------------------------------------------------------------------------
profileFile = 'bfield_protoMPEX.nc';

% Window cylinder
R_window = 0.0626;              % [m]

% Seed point for traced field line
r0 = 0.062865;                  % [m]
z0 = 1.75915;                    % [m]
theta0_deg = 0.0;               % [deg]

% Tracing controls
ds_trace = 2.0e-4;              % [m]
nStepsMax = 50000;
stop_if_outside_grid = true;
try_both_directions = true;

% Flux contour level to compare
% 'seed_point'  -> contour through (r0,z0)
% 'manual'      -> use xi_manual
xiMode = 'seed_point';
xi_manual = 1.0;

% Plot controls
makePlots = true;

%% ------------------------------------------------------------------------
% STEP 1. LOAD B FIELD
% -------------------------------------------------------------------------
fprintf('Reading %s ...\n', profileFile);

r = double(ncread(profileFile,'r')); r = r(:);
z = double(ncread(profileFile,'z')); z = z(:);

Br = double(ncread(profileFile,'br'));
Bz = double(ncread(profileFile,'bz'));

hasBt = true;
try
    Bt = double(ncread(profileFile,'bt'));
catch
    hasBt = false;
    Bt = zeros(size(Bz));
end

nR = numel(r);
nZ = numel(z);

if ~isequal(size(Br), [nR, nZ])
    Br = permute(Br, [2 1]);
end
if ~isequal(size(Bz), [nR, nZ])
    Bz = permute(Bz, [2 1]);
end
if ~isequal(size(Bt), [nR, nZ])
    Bt = permute(Bt, [2 1]);
end

fprintf('Grid: nR = %d, nZ = %d\n', nR, nZ);
fprintf('Domain: r = [%g, %g], z = [%g, %g]\n', min(r), max(r), min(z), max(z));

FBr = griddedInterpolant({r,z}, Br, 'linear', 'none');
FBz = griddedInterpolant({r,z}, Bz, 'linear', 'none');
FBt = griddedInterpolant({r,z}, Bt, 'linear', 'none');

%% ------------------------------------------------------------------------
% STEP 2. BUILD FLUX FUNCTION PSI(r,z) FROM Bz
%
% Axisymmetry:
%   Bz = (1/r) dpsi/dr
% so
%   dpsi/dr = r * Bz
%
% We integrate radially from r = 0 outward at each z.
% This reconstructs psi up to an additive function of z; for contour shape
% comparison this is sufficient in practice on a single connected domain.
% -------------------------------------------------------------------------
fprintf('Building approximate axisymmetric flux function psi(r,z) ...\n');

psi = zeros(nR, nZ);

for iz = 1:nZ
    integrand = r .* Bz(:,iz);
    psi(:,iz) = cumtrapz(r, integrand);
end

% Normalize by seed-point value if desired
Fpsi = griddedInterpolant({r,z}, psi, 'linear', 'none');

psi_seed = Fpsi(r0, z0);
if ~isfinite(psi_seed)
    error('Seed point (r0,z0) is outside interpolation support.');
end

switch lower(xiMode)
    case 'seed_point'
        psi_level = psi_seed;
    case 'manual'
        psi_level = xi_manual * psi_seed;
    otherwise
        error('Unknown xiMode.');
end

fprintf('psi(seed)   = %.6e\n', psi_seed);
fprintf('psi(level)  = %.6e\n', psi_level);

%% ------------------------------------------------------------------------
% STEP 3. EXTRACT CONTOUR AND FIND TOUCH/CROSSING WITH r = R_window
% -------------------------------------------------------------------------
fprintf('Extracting contour and intersecting with r = R_window ...\n');

% contourc expects x = z, y = r, data = psi(r,z)
C = contourc(z, r, psi, [psi_level psi_level]);

segments = parseContourMatrix(C);

if isempty(segments)
    error('No contour segments found at requested level.');
end

% Choose the contour segment closest to the seed point
bestSeg = [];
bestDist = inf;

for k = 1:numel(segments)
    zk = segments{k}(1,:);
    rk = segments{k}(2,:);

    d2 = (zk - z0).^2 + (rk - r0).^2;
    dmin = min(d2);

    if dmin < bestDist
        bestDist = dmin;
        bestSeg = segments{k};
    end
end

z_cont = bestSeg(1,:);
r_cont = bestSeg(2,:);

% Find all crossings with r = R_window
z_cross_all = [];
for i = 1:numel(r_cont)-1
    r1 = r_cont(i);
    r2 = r_cont(i+1);
    z1 = z_cont(i);
    z2 = z_cont(i+1);

    crossed = ((r1 - R_window) * (r2 - R_window) <= 0) && abs(r2 - r1) > 1e-12;
    if crossed
        f = (R_window - r1) / (r2 - r1);
        if f >= 0 && f <= 1
            zc = z1 + f * (z2 - z1);
            z_cross_all(end+1,1) = zc; %#ok<SAGROW>
        end
    end
end

if isempty(z_cross_all)
    warning('No contour crossing found with r = R_window.');
    z_contour_hit = NaN;
else
    [~, idxBestCross] = min(abs(z_cross_all - z0));
    z_contour_hit = z_cross_all(idxBestCross);
end

r_contour_hit = R_window;

%% ------------------------------------------------------------------------
% STEP 4. TRACE FIELD LINE TO CYLINDER
% -------------------------------------------------------------------------
fprintf('Tracing field line to cylinder ...\n');

[ok_fwd, r_hit_fwd, th_hit_fwd, z_hit_fwd, r_path_fwd, z_path_fwd] = ...
    trace_to_cyl_window_fast(r0, theta0_deg, z0, +1, ...
    FBr, FBt, FBz, R_window, ds_trace, nStepsMax, r, z, stop_if_outside_grid);

ok = ok_fwd;
r_hit = r_hit_fwd;
th_hit = th_hit_fwd;
z_hit = z_hit_fwd;
r_path = r_path_fwd;
z_path = z_path_fwd;
trace_direction = +1;

if (~ok_fwd) && try_both_directions
    [ok_bwd, r_hit_bwd, th_hit_bwd, z_hit_bwd, r_path_bwd, z_path_bwd] = ...
        trace_to_cyl_window_fast(r0, theta0_deg, z0, -1, ...
        FBr, FBt, FBz, R_window, ds_trace, nStepsMax, r, z, stop_if_outside_grid);

    if ok_bwd
        ok = true;
        r_hit = r_hit_bwd;
        th_hit = th_hit_bwd;
        z_hit = z_hit_bwd;
        r_path = r_path_bwd;
        z_path = z_path_bwd;
        trace_direction = -1;
    end
end

if ok
    fprintf('Traced hit found:\n');
    fprintf('  direction = %d\n', trace_direction);
    fprintf('  r_hit     = %.8f m\n', r_hit);
    fprintf('  z_hit     = %.8f m\n', z_hit);
    fprintf('  theta_hit = %.8f deg\n', th_hit);
else
    fprintf('No traced cylindrical hit found.\n');
end

%% ------------------------------------------------------------------------
% STEP 5. REPORT COMPARISON
% -------------------------------------------------------------------------
fprintf('\n===== COMPARISON =====\n');
fprintf('Seed point:\n');
fprintf('  r0 = %.8f m\n', r0);
fprintf('  z0 = %.8f m\n', z0);

fprintf('\nContour-based crossing with r = R_window:\n');
fprintf('  r_contour_hit = %.8f m\n', r_contour_hit);
fprintf('  z_contour_hit = %.8f m\n', z_contour_hit);

if ok && isfinite(z_contour_hit)
    fprintf('\nDifference:\n');
    fprintf('  dz = z_traced - z_contour = %.8e m\n', z_hit - z_contour_hit);
    fprintf('  dr = r_traced - r_contour = %.8e m\n', r_hit - r_contour_hit);
end

%% ------------------------------------------------------------------------
% STEP 6. PLOTS
% -------------------------------------------------------------------------
if makePlots
    figure('Color','w','Position',[100 100 1100 450]);

    subplot(1,2,1);
    hold on;

    % B magnitude background
    Bmag = sqrt(Br.^2 + Bz.^2 + Bt.^2);
    imagesc(z, r, Bmag);
    set(gca,'YDir','normal');

    % contour segment
    plot(z_cont, r_cont, 'k-', 'LineWidth', 2);

    % cylinder line
    plot([min(z) max(z)], R_window*[1 1], 'r--', 'LineWidth', 1.5);

    % seed point
    plot(z0, r0, 'wo', 'MarkerFaceColor','w', 'MarkerSize',7);

    % contour crossing
    if isfinite(z_contour_hit)
        plot(z_contour_hit, r_contour_hit, 'co', 'MarkerFaceColor','c', 'MarkerSize',8);
    end

    % traced path
    if ok
        plot(z_path, r_path, 'm-', 'LineWidth', 1.5);
        plot(z_hit, r_hit, 'yo', 'MarkerFaceColor','y', 'MarkerSize',8);
    end

    xlabel('z [m]');
    ylabel('r [m]');
    title('Contour vs traced cylindrical crossing');
    colorbar;
    legend({'psi contour','r = R_{window}','seed point','contour hit','traced path','traced hit'}, ...
        'Location','best');
    grid on;
    axis tight;

    subplot(1,2,2);
    hold on;
    plot(z_cont, r_cont - R_window, 'k-', 'LineWidth', 2);
    yline(0, 'r--', 'LineWidth', 1.5);

    if isfinite(z_contour_hit)
        plot(z_contour_hit, 0, 'co', 'MarkerFaceColor','c', 'MarkerSize',8);
    end

    if ok
        plot(z_path, r_path - R_window, 'm-', 'LineWidth', 1.5);
        plot(z_hit, r_hit - R_window, 'yo', 'MarkerFaceColor','y', 'MarkerSize',8);
    end

    xlabel('z [m]');
    ylabel('r - R_{window} [m]');
    title('Crossing diagnostic');
    grid on;
    axis tight;
end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================
function segments = parseContourMatrix(C)
    segments = {};
    idx = 1;
    while idx < size(C,2)
        npts = C(2,idx);
        block = C(:, idx+1:idx+npts);
        segments{end+1} = block; %#ok<AGROW>
        idx = idx + npts + 1;
    end
end

function [ok, r_hit, th_hit, z_hit, r_path, z_path] = trace_to_cyl_window_fast( ...
    r0, th0_deg, z0, zsign, FBr, FBt, FBz, R_window, ds_trace, nStepsMax, ...
    r_grid, z_grid, stop_if_outside_grid)

    ok = false;
    r_hit = NaN;
    th_hit = NaN;
    z_hit = NaN;

    r_path = NaN(1,nStepsMax);
    z_path = NaN(1,nStepsMax);

    r_now = r0;
    th_now_deg = th0_deg;
    z_now = z0;

    r_path(1) = r_now;
    z_path(1) = z_now;

    for is = 2:nStepsMax
        br_now = FBr(r_now, z_now);
        bt_now = FBt(r_now, z_now);
        bz_now = FBz(r_now, z_now);

        if any(~isfinite([br_now, bt_now, bz_now]))
            r_path = r_path(1:is-1);
            z_path = z_path(1:is-1);
            return;
        end

        bmag_now = sqrt(br_now^2 + bt_now^2 + bz_now^2);
        if bmag_now < 1e-12
            r_path = r_path(1:is-1);
            z_path = z_path(1:is-1);
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

        crossed = ((r_now - R_window) * (r_next - R_window) <= 0) && abs(r_next - r_now) > 1e-12;
        if crossed
            f = (R_window - r_now) / (r_next - r_now);
            r_hit = R_window;
            th_hit = mod(th_now_deg + f * (th_next_deg - th_now_deg), 360);
            z_hit = z_now + f * (z_next - z_now);

            r_path(is) = r_hit;
            z_path(is) = z_hit;

            r_path = r_path(1:is);
            z_path = z_path(1:is);
            ok = true;
            return;
        end

        if stop_if_outside_grid
            if r_next < min(r_grid) || r_next > max(r_grid) || ...
               z_next < min(z_grid) || z_next > max(z_grid)
                r_path = r_path(1:is-1);
                z_path = z_path(1:is-1);
                return;
            end
        end

        r_now = r_next;
        th_now_deg = th_next_deg;
        z_now = z_next;

        r_path(is) = r_now;
        z_path(is) = z_now;
    end

    r_path = r_path(isfinite(r_path));
    z_path = z_path(isfinite(z_path));
end