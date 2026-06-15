%% Step_3a_buildAnnulusFromTubeWindowIntersections.m
% Build a fixed annular Q_lm map, restricted to the flux tubes that cover
% the helicon cylindrical window over a chosen axial range.
%
% PURPOSE
%   1) Read the measured/inferred annular window map from tube/window intersections
%   2) Keep only source bins whose mapped cylindrical-window hits lie in
%      the helicon window axial range z in [z_window_min_keep, z_window_max_keep]
%   3) Preserve Step_3 q_window_lm as Q_lm by default
%   4) Optionally project onto a full radial mesh; inner region is left as NaN
%   5) Save the source map (measured region only - no central fill)
%
% IMPORTANT
%   - This script is fully independent
%   - It does NOT load any previous 3a output
%   - It enforces helicon-window axial coverage using z_window_lm
%   - The default mode keeps the same heat-flux scale as Step_3:
%
%         Q_lm = q_window_lm
%
%   - Optional parallel conversion is available. For the cylindrical window,
%     the local surface normal is radial:
%
%         sin(alpha_w) = |Br| / |B|
%
%     so, only in qTransferMode = 'parallel_from_window_incidence',
%
%         q_parallel_lm = q_window_lm / sin(alpha_w)
%
% INPUT
%   1) window_bundle_qWindow_fromIntersections.mat
%   2) bfield_protoMPEX_shotSeries_1.nc
%
% OUTPUT
%   fullSource_Qparallel.mat
%
% -------------------------------------------------------------------------

clear; close all; clc;

%% ------------------------------------------------------------------------
% USER INPUTS
% -------------------------------------------------------------------------
inFile  = 'window_bundle_qWindow_fromIntersections.mat';
profileFile = 'bfield_protoMPEX_shotSeries_1.nc';
outFile = 'fullSource_Qparallel.mat';


% Fallback axial location used only if mapped window coordinates are absent
z_source = 1.75915;      % [m]

% Physical cylindrical window radius fallback
R_window_fallback = 0.0626;   % [m]

% Helicon cylindrical window axial range to keep
z_window_min_keep = 1.60;   % [m]
z_window_max_keep = 1.90;   % [m]
useWindowZMask = true;

% Require a minimum theta coverage for a measured row to be considered
minRowCoverageFraction = 0.70;

% Full source radial mesh
r_inner_full = 0.000;  % [m]
r_outer_full = 0.070;  % [m]
nL_full      = 24;

% Step_3-like default: keep the trusted q_window_lm heat-flux scale and
% radial binning. Use the optional parallel conversion/full mesh only when
% you intentionally want q_window_lm / sin(alpha_w).
qTransferMode = 'preserve_step3_q_lm';  % 'preserve_step3_q_lm' or 'parallel_from_window_incidence'
useStep3RadialMesh = true;

% Optional smoothing of final full source map
doRadialSmooth = false;
smoothPasses   = 2;

% Optional light smoothing of cylindrical incidence factor before division
doSmoothSinAlpha = false;
smoothSinAlphaPasses = 2;

% Optional smoothing of row-to-row banding in known Q_lm
doSmoothKnownRows = false;
smoothKnownRowsPasses = 2;

% Minimum allowed sin(alpha_w) = |Br|/|B| at the cylindrical window.
% Bins where the field crosses the window at a shallower angle than this
% are excluded (set to NaN) because:
%   - Br -> 0 at the helicon-coil mid-plane (perfectly axial field there)
%   - These near-tangential crossings deposit negligible power on the window
%   - Dividing q_window by a near-zero sin(alpha_w) amplifies noise
% Physically, sin_alpha_w_min ~ 0.05 means a minimum ~3 deg incidence angle.
% Raise toward 0.10 to be more conservative; lower (never below ~0.03) to
% include more bins.  The existing warning triggers at sin_alpha_w < 0.05.
sin_alpha_w_min = 0.05;

makePlots = true;

%% ------------------------------------------------------------------------
% STEP 0. LOAD EXISTING ANNULAR MAP
% -------------------------------------------------------------------------
fprintf('Reading %s ...\n', inFile);

S = load(inFile);

requiredVars = {'theta_cent_deg','theta_edges_deg','r_cent','r_edges'};
for k = 1:numel(requiredVars)
    if ~isfield(S, requiredVars{k})
        error('Missing variable %s in %s.', requiredVars{k}, inFile);
    end
end

theta_cent_deg  = double(S.theta_cent_deg(:)).';
theta_edges_deg = double(S.theta_edges_deg(:)).';
r_cent_in       = double(S.r_cent(:));
r_edges_in      = double(S.r_edges(:));

nM    = numel(theta_cent_deg);
nL_in = numel(r_cent_in);

q_window_lm_in = [];
sin_alpha_w_in = [];
valid_window_bin_mask = true(nL_in, nM);

% Optional helicon-window axial mask
if useWindowZMask
    if isfield(S, 'z_window_lm')
        z_window_lm_in = double(S.z_window_lm);

        if isfield(S, 'theta_window_deg_lm')
            theta_window_deg_lm_in = double(S.theta_window_deg_lm);
        else
            theta_window_deg_lm_in = repmat(theta_cent_deg, nL_in, 1);
        end

        if ~isequal(size(z_window_lm_in), [nL_in, nM])
            error('z_window_lm must be [%d x %d].', nL_in, nM);
        end

        valid_window_bin_mask = isfinite(z_window_lm_in) & ...
                                (z_window_lm_in >= z_window_min_keep) & ...
                                (z_window_lm_in <= z_window_max_keep);

        fprintf('Using helicon-window z mask: [%g, %g] m\n', ...
            z_window_min_keep, z_window_max_keep);
        fprintf('Valid window-connected bins = %d / %d\n', ...
            nnz(valid_window_bin_mask), numel(valid_window_bin_mask));
    else
        warning('z_window_lm not found in input file. Cannot enforce helicon window z-range.');
    end
end

%% ------------------------------------------------------------------------
% STEP 1. GET / CONSTRUCT Q_lm ON THE STEP_3 BINNING
% -------------------------------------------------------------------------
if isfield(S, 'q_window_lm')
    fprintf('Using Step_3 q_window_lm as the trusted Q_lm input.\n');
    q_window_lm_in = double(S.q_window_lm);

    if ~isequal(size(q_window_lm_in), [nL_in, nM])
        error('q_window_lm must be [%d x %d].', nL_in, nM);
    end

    % Load magnetic field
    r = double(ncread(profileFile, 'r')); r = r(:);
    z = double(ncread(profileFile, 'z')); z = z(:);

    Br = double(ncread(profileFile, 'br'));
    Bt = double(ncread(profileFile, 'bt'));
    Bz = double(ncread(profileFile, 'bz'));

    nR = numel(r);
    nZ = numel(z);

    if ~isequal(size(Br), [nR, nZ])
        Br = permute(Br, [2 1]);
        Bt = permute(Bt, [2 1]);
        Bz = permute(Bz, [2 1]);
    end

    FBr = griddedInterpolant({r, z}, Br, 'linear', 'none');
    FBt = griddedInterpolant({r, z}, Bt, 'linear', 'none');
    FBz = griddedInterpolant({r, z}, Bz, 'linear', 'none');

    % Use mapped cylindrical window hit coordinates if available
    if isfield(S, 'r_window_lm') && isfield(S, 'z_window_lm')
        fprintf('Using mapped window hit locations (r_window_lm, z_window_lm) for incidence angle.\n');

        r_eval_w = double(S.r_window_lm);
        z_eval_w = double(S.z_window_lm);

        if ~isequal(size(r_eval_w), [nL_in, nM]) || ~isequal(size(z_eval_w), [nL_in, nM])
            error('r_window_lm and z_window_lm must both be [%d x %d].', nL_in, nM);
        end
    else
        fprintf('Mapped window hit locations not found. Falling back to fixed cylindrical location.\n');
        [~, R_LM_IN] = meshgrid(theta_cent_deg, r_cent_in);
        r_eval_w = R_window_fallback * ones(size(R_LM_IN));
        z_eval_w = z_source * ones(size(R_LM_IN));
    end

    % Fill missing mapped coordinates with fallback values
    badCoord = ~isfinite(r_eval_w) | ~isfinite(z_eval_w);
    if any(badCoord(:))
        fprintf('Mapped window coordinates contain NaNs. Filling missing entries with fallback cylindrical location.\n');
        r_eval_w(badCoord) = R_window_fallback;
        z_eval_w(badCoord) = z_source;
    end

    br_w = FBr(r_eval_w, z_eval_w);
    bt_w = FBt(r_eval_w, z_eval_w);
    bz_w = FBz(r_eval_w, z_eval_w);

    Bmag_w = sqrt(br_w.^2 + bt_w.^2 + bz_w.^2);
    Bmag_w(Bmag_w < 1e-12) = 1e-12;

    % Cylindrical window: normal = e_r
    sin_alpha_w_in = abs(br_w) ./ Bmag_w;

    % Optional light radial smoothing of sin(alpha_w)
    if doSmoothSinAlpha
        fprintf('Applying light radial smoothing to sin(alpha_w) ...\n');
        for pass = 1:smoothSinAlphaPasses
            SinTmp = sin_alpha_w_in;
            for il = 2:nL_in-1
                vals = sin_alpha_w_in(il-1:il+1,:);
                good = isfinite(vals);
                vals(~good) = 0;
                num = 0.25 * vals(1,:) + 0.50 * vals(2,:) + 0.25 * vals(3,:);
                den = 0.25 * good(1,:) + 0.50 * good(2,:) + 0.25 * good(3,:);
                upd = den > 0;
                SinTmp(il,upd) = num(upd) ./ den(upd);
            end
            sin_alpha_w_in = SinTmp;
        end
    end

    finiteSin = sin_alpha_w_in(isfinite(sin_alpha_w_in));
    if ~isempty(finiteSin)
        fprintf('Window sin(alpha_w) min/median/max = %.6g / %.6g / %.6g\n', ...
            min(finiteSin), median(finiteSin), max(finiteSin));
    end

    switch lower(qTransferMode)
        case 'preserve_step3_q_lm'
            q_parallel_lm_in = q_window_lm_in;
            fprintf(['qTransferMode = preserve_step3_q_lm: preserving Step_3 ', ...
                     'q_window_lm scale; no division by sin(alpha_w).\n']);

        case 'parallel_from_window_incidence'
            % Exclude near-tangential bins before computing q_parallel.
            % Below this threshold, dividing by sin(alpha_w) can amplify
            % measurement noise into unphysical MW/m^2 spikes.
            tooShallow = isfinite(sin_alpha_w_in) & (sin_alpha_w_in < sin_alpha_w_min);
            n_excluded = sum(tooShallow(:));
            n_valid    = sum(isfinite(sin_alpha_w_in(:)));
            if n_excluded > 0
                fprintf('[sin_alpha_w] Excluding %d / %d valid bins with sin(alpha_w) < %.4g  (angle < %.1f deg)\n', ...
                    n_excluded, n_valid, sin_alpha_w_min, rad2deg(asin(sin_alpha_w_min)));
                sin_alpha_w_in(tooShallow) = NaN;   % exclude - propagates to q_parallel
            end

            sin_alpha_w_safe_in = max(sin_alpha_w_in, sin_alpha_w_min);
            q_parallel_lm_in = q_window_lm_in ./ sin_alpha_w_safe_in;
            fprintf('qTransferMode = parallel_from_window_incidence: Q_lm = q_window_lm / sin(alpha_w).\n');

        otherwise
            error('Unknown qTransferMode: %s', qTransferMode);
    end

    % Keep only bins that actually hit selected helicon-window z-range
    if useWindowZMask
        q_window_lm_in(~valid_window_bin_mask)   = NaN;
        sin_alpha_w_in(~valid_window_bin_mask)   = NaN;
        q_parallel_lm_in(~valid_window_bin_mask) = NaN;
    end

elseif isfield(S, 'q_parallel_lm')
    fprintf('q_window_lm not found. Using q_parallel_lm directly from input file.\n');
    q_parallel_lm_in = double(S.q_parallel_lm);

else
    error('Input file must contain either q_parallel_lm or q_window_lm.');
end

if ~isequal(size(q_parallel_lm_in), [nL_in, nM])
    error('Size mismatch: q_parallel_lm_in must be [%d x %d].', nL_in, nM);
end

finiteQ = q_parallel_lm_in(isfinite(q_parallel_lm_in));
if ~isempty(finiteQ)
    fprintf('[DIAG] Q_lm carried forward: min=%.3e  median=%.3e  p95=%.3e  max=%.3e  W/m^2\n', ...
        min(finiteQ), median(finiteQ), prctile(finiteQ,95), max(finiteQ));
    fprintf('[DIAG] Q_lm carried forward: min=%.1f  median=%.1f  p95=%.1f  max=%.1f  kW/m^2\n', ...
        min(finiteQ)/1e3, median(finiteQ)/1e3, prctile(finiteQ,95)/1e3, max(finiteQ)/1e3);
end

%% ------------------------------------------------------------------------
% STEP 1b. ENFORCE ROW COVERAGE BEFORE PROJECTING TO FULL SOURCE MESH
% -------------------------------------------------------------------------
rowCoverageFraction = mean(isfinite(q_parallel_lm_in), 2);
rowHasEnoughCoverage = rowCoverageFraction >= minRowCoverageFraction;

fprintf('Row coverage fraction min/median/max = %.3f / %.3f / %.3f\n', ...
    min(rowCoverageFraction), median(rowCoverageFraction), max(rowCoverageFraction));
fprintf('Rows passing minRowCoverageFraction=%.2f: %d / %d\n', ...
    minRowCoverageFraction, nnz(rowHasEnoughCoverage), nL_in);

if ~any(rowHasEnoughCoverage)
    error('No measured rows satisfy minRowCoverageFraction=%.2f after masks/incidence filtering.', ...
        minRowCoverageFraction);
end

lowCoverageRows = ~rowHasEnoughCoverage;
if any(lowCoverageRows)
    fprintf('Dropping %d low-coverage measured rows before radial interpolation.\n', ...
        nnz(lowCoverageRows));
    q_parallel_lm_in(lowCoverageRows,:) = NaN;
    if ~isempty(q_window_lm_in)
        q_window_lm_in(lowCoverageRows,:) = NaN;
    end
    if ~isempty(sin_alpha_w_in)
        sin_alpha_w_in(lowCoverageRows,:) = NaN;
    end
    if exist('valid_window_bin_mask','var')
        valid_window_bin_mask(lowCoverageRows,:) = false;
    end
end

% Optional: smooth row-to-row radial banding in q_parallel_lm_in
% while preserving theta pattern within each row. This is intentionally done
% after the coverage filter so sparse rows cannot bias neighboring row means.
if doSmoothKnownRows
    fprintf('Applying light smoothing to row means of q_parallel_lm_in ...\n');

    q_mean_row = mean(q_parallel_lm_in, 2, 'omitnan');
    q_mean_row_s = q_mean_row;

    for pass = 1:smoothKnownRowsPasses
        qtmp = q_mean_row_s;
        for il = 2:nL_in-1
            nbr = q_mean_row_s(il-1:il+1);
            if any(isfinite(nbr))
                w = [0.25; 0.50; 0.25];
                good = isfinite(nbr);
                qtmp(il) = sum(w(good).*nbr(good)) / sum(w(good));
            end
        end
        q_mean_row_s = qtmp;
    end

    q_parallel_smoothed = q_parallel_lm_in;
    for il = 1:nL_in
        row_mean = mean(q_parallel_lm_in(il,:), 'omitnan');
        if isfinite(row_mean) && row_mean > 0
            theta_shape = q_parallel_lm_in(il,:) / row_mean;
            theta_shape(~isfinite(theta_shape)) = NaN;
            q_parallel_smoothed(il,:) = q_mean_row_s(il) * theta_shape;
        end
    end

    q_parallel_lm_in = q_parallel_smoothed;
end

fprintf('Loaded/source Q_lm size = [%d x %d]\n', nL_in, nM);
fprintf('Measured radial range = [%g, %g] m\n', min(r_edges_in), max(r_edges_in));

%% ------------------------------------------------------------------------
% STEP 2. BUILD FULL RADIAL SOURCE MESH
% -------------------------------------------------------------------------
fprintf('Building full source radial mesh ...\n');

if useStep3RadialMesh
    r_edges_full = r_edges_in(:).';
    r_cent_full  = r_cent_in(:).';
    nL_full      = nL_in;
    fprintf('Using Step_3 radial mesh directly (nL = %d).\n', nL_full);
else
    r_edges_full = linspace(r_inner_full, r_outer_full, nL_full+1);
    r_cent_full  = 0.5 * (r_edges_full(1:end-1) + r_edges_full(2:end));
end

if r_outer_full < max(r_edges_in) - 1e-12
    error('r_outer_full must be >= outer radius of measured map.');
end

Q_parallel_full_lm = NaN(nL_full, nM);

%% ------------------------------------------------------------------------
% STEP 3. COPY / INTERPOLATE KNOWN OUTER REGION
% -------------------------------------------------------------------------
if useStep3RadialMesh
    fprintf('Copying Step_3 Q_lm directly onto the output mesh ...\n');
else
    fprintf('Projecting known Q_lm map onto outer region of full source mesh ...\n');
end

if useStep3RadialMesh
    Q_parallel_full_lm = q_parallel_lm_in;
else
    for im = 1:nM
        col_in = q_parallel_lm_in(:,im);
        good = isfinite(col_in);

        if nnz(good) >= 2
            Q_parallel_full_lm(:,im) = interp1(r_cent_in(good), col_in(good), r_cent_full, 'linear', NaN);
        elseif nnz(good) == 1
            Q_parallel_full_lm(:,im) = NaN(size(r_cent_full));
            [~, idxNearest] = min(abs(r_cent_full - r_cent_in(good)));
            Q_parallel_full_lm(idxNearest,im) = col_in(good);
        else
            Q_parallel_full_lm(:,im) = NaN(size(r_cent_full));
        end
    end
end

measuredMaskRows = r_cent_full >= min(r_cent_in) & r_cent_full <= max(r_cent_in);
innerMaskRows    = r_cent_full <  min(r_cent_in);

fprintf('Rows inside known radial range = %d / %d\n', nnz(measuredMaskRows), nL_full);
fprintf('Rows without measurement (left as NaN) = %d / %d\n', nnz(innerMaskRows), nL_full);

%% ------------------------------------------------------------------------
% STEP 4. OPTIONAL RADIAL SMOOTHING WITHIN MEASURED REGION ONLY
% -------------------------------------------------------------------------
if doRadialSmooth
    fprintf('Applying radial smoothing (NaN-safe, measured region only) ...\n');
    for pass = 1:smoothPasses
        Qtmp = Q_parallel_full_lm;
        for il = 2:nL_full-1
            prev = Q_parallel_full_lm(il-1,:);
            curr = Q_parallel_full_lm(il,:);
            next = Q_parallel_full_lm(il+1,:);
            prevFinite = isfinite(prev);
            currFinite = isfinite(curr);
            nextFinite = isfinite(next);
            w_prev = 0.25 * prevFinite;
            w_curr = 0.50 * currFinite;
            w_next = 0.25 * nextFinite;
            prev(~isfinite(prev)) = 0;
            curr(~isfinite(curr)) = 0;
            next(~isfinite(next)) = 0;
            wsum = w_prev + w_curr + w_next;
            blend = (w_prev .* prev + w_curr .* curr + w_next .* next) ./ wsum;
            updateMask = currFinite & wsum > 0;
            Qtmp(il, updateMask) = blend(updateMask);
        end
        Q_parallel_full_lm = Qtmp;
    end
end

%% ------------------------------------------------------------------------
% STEP 9. BUILD FULL SOURCE PATCH AREAS
% -------------------------------------------------------------------------
fprintf('Building full source patch areas ...\n');

A_full_lm = zeros(nL_full, nM);
dtheta_rad = deg2rad(theta_edges_deg(2) - theta_edges_deg(1));

for il = 1:nL_full
    for im = 1:nM
        A_full_lm(il,im) = 0.5 * (r_edges_full(il+1)^2 - r_edges_full(il)^2) * dtheta_rad;
    end
end

%% ------------------------------------------------------------------------
% STEP 10. DIAGNOSTICS
% -------------------------------------------------------------------------
fprintf('Full source map size = [%d x %d]\n', size(Q_parallel_full_lm,1), size(Q_parallel_full_lm,2));
fprintf('Min Q_full_lm = %g\n', min(Q_parallel_full_lm(:), [], 'omitnan'));
fprintf('Max Q_full_lm = %g\n', max(Q_parallel_full_lm(:), [], 'omitnan'));
fprintf('Mean Q_full_lm = %g\n', mean(Q_parallel_full_lm(:), 'omitnan'));

P_source_full = sum(Q_parallel_full_lm(:) .* A_full_lm(:), 'omitnan');
P_parallel_full = P_source_full;  % Backward-compatible alias for older scripts.
fprintf('Approx integrated Q_lm*A_source = %.6e\n', P_source_full);

%% ------------------------------------------------------------------------
% STEP 11. PLOTS
% -------------------------------------------------------------------------
if makePlots
    figure('Color','w','Position',[100 100 1100 420]);

    subplot(1,2,1);
    imagesc(theta_cent_deg, r_cent_in, q_parallel_lm_in);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('known r [m]');
    title('Known Q_{lm}');
    colorbar;

    subplot(1,2,2);
    imagesc(theta_cent_deg, r_cent_full, Q_parallel_full_lm);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('full source r [m]');
    title('Output Q_{lm}');
    colorbar;

    if ~isempty(q_window_lm_in) && ~isempty(sin_alpha_w_in)
        figure('Color','w','Position',[120 120 1400 420]);

        subplot(1,3,1);
        imagesc(theta_cent_deg, r_cent_in, q_window_lm_in);
        set(gca,'YDir','normal');
        xlabel('\theta [deg]');
        ylabel('known r [m]');
        title('Input q_{window,lm}');
        colorbar;

        subplot(1,3,2);
        imagesc(theta_cent_deg, r_cent_in, sin_alpha_w_in);
        set(gca,'YDir','normal');
        xlabel('\theta [deg]');
        ylabel('known r [m]');
        title('\sin(\alpha_w) on cylindrical window');
        colorbar;

        subplot(1,3,3);
        imagesc(theta_cent_deg, r_cent_in, q_parallel_lm_in);
        set(gca,'YDir','normal');
        xlabel('\theta [deg]');
        ylabel('known r [m]');
        title('Carried Q_{lm}');
        colorbar;
    end

    figure('Color','w');
    plot(r_cent_full, mean(Q_parallel_full_lm,2,'omitnan'), 'LineWidth', 2); hold on;
    plot(r_cent_in, mean(q_parallel_lm_in,2,'omitnan'), 'o-', 'LineWidth', 1.5);
    xlabel('r [m]');
    ylabel('<Q>_\theta [W/m^2]');
    title('Radial mean profile: full source vs known rows');
    legend('Output mesh','Known rows','Location','best');
    grid on;

    [THF_RAD, RF] = meshgrid(deg2rad(theta_cent_deg), r_cent_full);
    XF = RF .* cos(THF_RAD);
    YF = RF .* sin(THF_RAD);

    figure('Color','w');
    scatter(XF(:), YF(:), 40, Q_parallel_full_lm(:), 'filled');
    axis equal;
    xlabel('x [m]');
    ylabel('y [m]');
    title('Output Q_{lm} in x-y');
    colorbar;

    if useWindowZMask && exist('z_window_lm_in','var')
        figure('Color','w','Position',[140 140 1100 420]);

        subplot(1,2,1);
        imagesc(theta_cent_deg, r_cent_in, valid_window_bin_mask);
        set(gca,'YDir','normal');
        xlabel('\theta [deg]');
        ylabel('known r [m]');
        title('Helicon-window valid-bin mask');
        colorbar;

        subplot(1,2,2);
        imagesc(theta_cent_deg, r_cent_in, z_window_lm_in);
        set(gca,'YDir','normal');
        xlabel('\theta [deg]');
        ylabel('known r [m]');
        title('Mapped cylindrical hit z_{window,lm}');
        colorbar;
    end
end

%% ================================================================
% THETA-Z PLOT OF WINDOW HEAT FLUX (TRUE CYLINDRICAL MAP)
% ================================================================
% Use already-loaded variables - no file reload needed

if makePlots && ~isempty(q_window_lm_in) && ...
        exist('theta_window_deg_lm_in','var') && exist('z_window_lm_in','var')

    theta_scatter_plot = theta_window_deg_lm_in(:);
    z_scatter_plot     = z_window_lm_in(:);
    q_scatter_plot     = q_window_lm_in(:);

    valid_plot = isfinite(theta_scatter_plot) & isfinite(z_scatter_plot) & isfinite(q_scatter_plot);

    theta_scatter_plot = theta_scatter_plot(valid_plot);
    z_scatter_plot     = z_scatter_plot(valid_plot);
    q_scatter_plot     = q_scatter_plot(valid_plot);

    if numel(q_scatter_plot) >= 3
        theta_ext_plot = [theta_scatter_plot; theta_scatter_plot-360; theta_scatter_plot+360];
        z_ext_plot     = [z_scatter_plot; z_scatter_plot; z_scatter_plot];
        q_ext_plot     = [q_scatter_plot; q_scatter_plot; q_scatter_plot];

        theta_grid_tz = linspace(0, 360, 360);
        z_grid_tz     = linspace(min(z_ext_plot), max(z_ext_plot), 200);

        [THG_tz, ZG_tz] = meshgrid(theta_grid_tz, z_grid_tz);

        q_theta_z = griddata(theta_ext_plot, z_ext_plot, q_ext_plot, THG_tz, ZG_tz, 'linear');

        bad_tz = ~isfinite(q_theta_z);
        if any(bad_tz(:))
            q_theta_z(bad_tz) = griddata(theta_ext_plot, z_ext_plot, q_ext_plot, THG_tz(bad_tz), ZG_tz(bad_tz), 'nearest');
        end

        figure('Color','w','Position',[100 100 900 400]);

        imagesc(theta_grid_tz, z_grid_tz, q_theta_z);
        set(gca,'YDir','normal');

        xlabel('\theta [deg]');
        ylabel('window hit z [m]');
        title('q^W(\theta,z) from flux-tube mapping');

        xlim([0 360]);
        ylim([min(z_grid_tz) max(z_grid_tz)]);

        clim([0 max(q_theta_z(:),[],'omitnan')]);
        colorbar;
    end
end

%% ------------------------------------------------------------------------
% STEP 12. SAVE OUTPUT
% -------------------------------------------------------------------------
save(outFile, ...
    'Q_parallel_full_lm', ...
    'r_cent_full', ...
    'r_edges_full', ...
    'A_full_lm', ...
    'theta_cent_deg', ...
    'theta_edges_deg', ...
    'q_parallel_lm_in', ...
    'q_window_lm_in', ...
    'sin_alpha_w_in', ...
    'r_cent_in', ...
    'r_edges_in', ...
    'P_source_full', ...
    'P_parallel_full', ...
    'z_source', ...
    'z_window_min_keep', ...
    'z_window_max_keep', ...
    'useWindowZMask', ...
    'qTransferMode', ...
    'useStep3RadialMesh', ...
    'valid_window_bin_mask', ...
    'minRowCoverageFraction', ...
    'rowCoverageFraction', ...
    'rowHasEnoughCoverage');

fprintf('Saved %s\n', outFile);
