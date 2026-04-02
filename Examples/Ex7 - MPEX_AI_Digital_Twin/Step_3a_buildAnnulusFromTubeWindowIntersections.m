%% Step_3aa_buildAnnulusFromTubeWindowIntersections.m
% Build a full upstream annular source map of fixed tube-carried quantity
% Q_parallel(l,m), restricted to the flux tubes that cover the helicon
% cylindrical window over a chosen axial range.
%
% PURPOSE
%   1) Read the measured/inferred annular window map from tube/window intersections
%   2) Keep only source bins whose mapped cylindrical-window hits lie in
%      the helicon window axial range z in [z_window_min_keep, z_window_max_keep]
%   3) Convert q_window_lm -> q_parallel_lm using cylindrical-window incidence
%   4) Build the outer measured region on the full radial mesh
%   5) Fill the inner region using the SAME functional form as 3a
%   6) Use a controlled anchor level so the inner radial decay matches the
%      old-style smooth decay rather than amplified converted edge values
%   7) Save the combined full source map
%
% IMPORTANT
%   - This script is fully independent
%   - It does NOT load any previous 3a output
%   - It enforces helicon-window axial coverage using z_window_lm
%   - For the cylindrical window, the local surface normal is radial:
%
%         sin(alpha_w) = |Br| / |B|
%
%     so
%
%         q_parallel_lm = q_window_lm / sin(alpha_w)
%
% INPUT
%   1) window_bundle_qWindow_fromIntersections.mat
%   2) bfield_protoMPEX_shotSeries_1.nc
%
% OUTPUT
%   fullSource_Qparallel_3aa.mat
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
z_window_max_keep = 2.00;   % [m]
useWindowZMask = true;

% Require a minimum theta coverage for a measured row to be considered
minRowCoverageFraction = 0.70;

% Full source radial mesh
r_inner_full = 0.000;  % [m]
r_outer_full = 0.070;  % [m]
nL_full      = 24;

% Fill method for central region:
%   'flat_to_center'
%   'linear_to_center'
%   'parabolic_to_center'
%   'exp_to_center'
fillMode = 'parabolic_to_center';

% How to choose center value at r = 0
%   'same_as_first_row_mean'
%   'fraction_of_first_row_mean'
%   'user_value'
centerValueMode = 'fraction_of_first_row_mean';

centerFraction = 0.8;
centerUserValue = 1.0e4;

% Preserve angular structure
preserveThetaPattern = true;

% Initial anchor measured row guess; actual anchor may be updated later
anchorMeasuredRow = 1;

% Optional smoothing of final full source map
doRadialSmooth = true;
smoothPasses   = 2;

% Optional light smoothing of cylindrical incidence factor before division
doSmoothSinAlpha = true;
smoothSinAlphaPasses = 2;

% Optional smoothing of row-to-row banding in known q_parallel_lm_in
doSmoothKnownRows = true;
smoothKnownRowsPasses = 2;

% Control inner-fill amplitude from several edge rows
nEdgeRowsForFill = 5;

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

if anchorMeasuredRow < 1 || anchorMeasuredRow > nL_in
    error('anchorMeasuredRow must be between 1 and %d.', nL_in);
end

q_window_lm_in = [];
sin_alpha_w_in = [];
valid_window_bin_mask = true(nL_in, nM);

% Optional helicon-window axial mask
if useWindowZMask
    if isfield(S, 'z_window_lm')
        z_window_lm_in = double(S.z_window_lm);

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
% STEP 1. GET / CONSTRUCT q_parallel_lm_in
% -------------------------------------------------------------------------
if isfield(S, 'q_parallel_lm')
    fprintf('Using q_parallel_lm directly from input file.\n');
    q_parallel_lm_in = double(S.q_parallel_lm);

elseif isfield(S, 'q_window_lm')
    fprintf('q_parallel_lm not found. Converting q_window_lm to q_parallel_lm using cylindrical-window incidence.\n');

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
                num = 0.25 * vals(1,:) + 0.50 * vals(2,:) + 0.25 * vals(3,:);
                den = 0.25 * good(1,:) + 0.50 * good(2,:) + 0.25 * good(3,:);
                upd = den > 0;
                SinTmp(il,upd) = num(upd) ./ den(upd);
            end
            sin_alpha_w_in = SinTmp;
        end
    end

    sin_alpha_w_safe_in = max(sin_alpha_w_in, 1e-6);
    q_parallel_lm_in = q_window_lm_in ./ sin_alpha_w_safe_in;

    % Keep only bins that actually hit selected helicon-window z-range
    if useWindowZMask
        q_window_lm_in(~valid_window_bin_mask)   = NaN;
        sin_alpha_w_in(~valid_window_bin_mask)   = NaN;
        q_parallel_lm_in(~valid_window_bin_mask) = NaN;
    end

    % Optional: smooth row-to-row radial banding in q_parallel_lm_in
    % while preserving theta pattern within each row
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

else
    error('Input file must contain either q_parallel_lm or q_window_lm.');
end

if ~isequal(size(q_parallel_lm_in), [nL_in, nM])
    error('Size mismatch: q_parallel_lm_in must be [%d x %d].', nL_in, nM);
end

fprintf('Loaded/source q_parallel_lm_in size = [%d x %d]\n', nL_in, nM);
fprintf('Measured radial range = [%g, %g] m\n', min(r_edges_in), max(r_edges_in));

%% ------------------------------------------------------------------------
% STEP 1b. CHOOSE ANCHOR ROW FROM ACTUAL HELICON COVERAGE
% -------------------------------------------------------------------------
rowCoverageFraction = mean(isfinite(q_parallel_lm_in), 2, 'omitnan');
candidateRows = find(rowCoverageFraction >= minRowCoverageFraction);

if isempty(candidateRows)
    warning('No row satisfies minRowCoverageFraction. Falling back to first finite row.');
    candidateRows = find(any(isfinite(q_parallel_lm_in), 2));
end

if isempty(candidateRows)
    error('No valid measured rows remain after applying helicon window z mask.');
end

anchorMeasuredRow = candidateRows(1);

fprintf('Selected anchorMeasuredRow = %d\n', anchorMeasuredRow);
fprintf('Anchor row coverage fraction = %.3f\n', rowCoverageFraction(anchorMeasuredRow));

%% ------------------------------------------------------------------------
% STEP 2. BUILD FULL RADIAL SOURCE MESH
% -------------------------------------------------------------------------
fprintf('Building full source radial mesh ...\n');

r_edges_full = linspace(r_inner_full, r_outer_full, nL_full+1);
r_cent_full  = 0.5 * (r_edges_full(1:end-1) + r_edges_full(2:end));

if r_outer_full < max(r_edges_in) - 1e-12
    error('r_outer_full must be >= outer radius of measured map.');
end

Q_parallel_full_lm = NaN(nL_full, nM);

%% ------------------------------------------------------------------------
% STEP 3. COPY / INTERPOLATE KNOWN OUTER REGION
% -------------------------------------------------------------------------
fprintf('Projecting known q_parallel map onto outer region of full source mesh ...\n');

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

measuredMaskRows = r_cent_full >= min(r_cent_in) & r_cent_full <= max(r_cent_in);
innerMaskRows    = r_cent_full <  min(r_cent_in);

fprintf('Rows inside known radial range = %d / %d\n', nnz(measuredMaskRows), nL_full);
fprintf('Rows requiring central fill    = %d / %d\n', nnz(innerMaskRows), nL_full);

%% ------------------------------------------------------------------------
% STEP 4. BUILD THETA PATTERN FROM ANCHOR ROW
% -------------------------------------------------------------------------
fprintf('Building anchor theta pattern for inner fill ...\n');

q_anchor_raw = q_parallel_lm_in(anchorMeasuredRow, :);
q_anchor_raw_mean = mean(q_anchor_raw, 'omitnan');

if ~isfinite(q_anchor_raw_mean) || q_anchor_raw_mean <= 0
    error('Anchor row has invalid mean value.');
end

thetaPattern = q_anchor_raw ./ q_anchor_raw_mean;
thetaPattern(~isfinite(thetaPattern)) = 1.0;

if ~preserveThetaPattern
    thetaPattern(:) = 1.0;
end

%% ------------------------------------------------------------------------
% STEP 5. CHOOSE FILL ANCHOR LEVEL
% -------------------------------------------------------------------------
fprintf('Choosing fill anchor level ...\n');

q_mean_measured = mean(q_parallel_lm_in, 2, 'omitnan');

nEdgeRowsForFill = min(nEdgeRowsForFill, nL_in);
candidateFillRows = candidateRows(1:min(nEdgeRowsForFill, numel(candidateRows)));
q_anchor_fill_mean = mean(q_mean_measured(candidateFillRows), 'omitnan');

if ~isfinite(q_anchor_fill_mean) || q_anchor_fill_mean <= 0
    q_anchor_fill_mean = q_anchor_raw_mean;
end

fprintf('Raw anchor row mean Q_parallel         = %.6e W/m^2\n', q_anchor_raw_mean);
fprintf('Controlled fill anchor mean Q_parallel = %.6e W/m^2\n', q_anchor_fill_mean);

switch lower(centerValueMode)
    case 'same_as_first_row_mean'
        q_center_mean = q_anchor_fill_mean;

    case 'fraction_of_first_row_mean'
        q_center_mean = centerFraction * q_anchor_fill_mean;

    case 'user_value'
        q_center_mean = centerUserValue;

    otherwise
        error('Unknown centerValueMode: %s', centerValueMode);
end

fprintf('Chosen center mean for fill            = %.6e W/m^2\n', q_center_mean);

% %% ------------------------------------------------------------------------
% % STEP 6. FILL CENTRAL REGION
% % -------------------------------------------------------------------------
% fprintf('Filling central region using %s ...\n', fillMode);
% 
% r_meas_anchor = r_cent_in(anchorMeasuredRow);
% 
% for il = 1:nL_full
%     if ~innerMaskRows(il)
%         continue;
%     end
% 
%     r_now = r_cent_full(il);
% 
%     if r_meas_anchor > 0
%         xi = r_now / r_meas_anchor;
%     else
%         xi = 0;
%     end
%     xi = max(0, min(1, xi));
% 
%     switch lower(fillMode)
%         case 'flat_to_center'
%             q_mean_now = q_anchor_fill_mean;
% 
%         case 'linear_to_center'
%             q_mean_now = q_center_mean + (q_anchor_fill_mean - q_center_mean) * xi;
% 
%         case 'parabolic_to_center'
%             q_mean_now = q_center_mean + (q_anchor_fill_mean - q_center_mean) * xi^2;
% 
%         case 'exp_to_center'
%             beta = 3.0;
%             q_mean_now = q_center_mean + (q_anchor_fill_mean - q_center_mean) * ...
%                          (exp(beta*xi) - 1) / (exp(beta) - 1);
% 
%         otherwise
%             error('Unknown fillMode: %s', fillMode);
%     end
% 
%     rowVals = q_mean_now * thetaPattern;
%     rowVals(~isfinite(rowVals)) = q_mean_now;
%     rowVals(rowVals < 0) = 0;
% 
%     Q_parallel_full_lm(il,:) = rowVals;
% end
% 
% %% ------------------------------------------------------------------------
% % STEP 7. FILL ANY REMAINING NaNs
% % -------------------------------------------------------------------------
% fprintf('Filling any remaining NaNs ...\n');
% 
% for im = 1:nM
%     col = Q_parallel_full_lm(:,im);
%     valid = isfinite(col);
% 
%     if all(~valid)
%         error('Column %d has no valid values after fill.', im);
%     end
% 
%     if any(~valid)
%         validRows = find(valid);
%         for il = find(~valid).'
%             [~, idx] = min(abs(validRows - il));
%             col(il) = col(validRows(idx));
%         end
%     end
% 
%     Q_parallel_full_lm(:,im) = col;
% end
% 
% %% ------------------------------------------------------------------------
% % STEP 8. OPTIONAL RADIAL SMOOTHING OF FINAL FULL MAP
% % -------------------------------------------------------------------------
% if doRadialSmooth
%     fprintf('Applying radial smoothing ...\n');
%     for pass = 1:smoothPasses
%         Qtmp = Q_parallel_full_lm;
%         for il = 2:nL_full-1
%             Qtmp(il,:) = 0.25 * Q_parallel_full_lm(il-1,:) + ...
%                          0.50 * Q_parallel_full_lm(il,:)   + ...
%                          0.25 * Q_parallel_full_lm(il+1,:);
%         end
%         Q_parallel_full_lm = Qtmp;
%     end
% end
% 
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
fprintf('Min Q_parallel_full = %g\n', min(Q_parallel_full_lm(:), [], 'omitnan'));
fprintf('Max Q_parallel_full = %g\n', max(Q_parallel_full_lm(:), [], 'omitnan'));
fprintf('Mean Q_parallel_full = %g\n', mean(Q_parallel_full_lm(:), 'omitnan'));

P_parallel_full = nansum(Q_parallel_full_lm(:) .* A_full_lm(:));
fprintf('Approx integrated source quantity = %.6e\n', P_parallel_full);

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
    title('Known q_{\parallel,lm}');
    colorbar;

    subplot(1,2,2);
    imagesc(theta_cent_deg, r_cent_full, Q_parallel_full_lm);
    set(gca,'YDir','normal');
    xlabel('\theta [deg]');
    ylabel('full source r [m]');
    title('Filled full source Q_{\parallel,lm}');
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
        title('Converted q_{\parallel,lm}');
        colorbar;
    end

    figure('Color','w');
    plot(r_cent_full, mean(Q_parallel_full_lm,2,'omitnan'), 'LineWidth', 2); hold on;
    plot(r_cent_in, mean(q_parallel_lm_in,2,'omitnan'), 'o-', 'LineWidth', 1.5);
    xlabel('r [m]');
    ylabel('<Q_{\parallel}>_\theta [W/m^2]');
    title('Radial mean profile: full source vs known rows');
    legend('Filled full source','Known rows','Location','best');
    grid on;

    [THF_RAD, RF] = meshgrid(deg2rad(theta_cent_deg), r_cent_full);
    XF = RF .* cos(THF_RAD);
    YF = RF .* sin(THF_RAD);

    figure('Color','w');
    scatter(XF(:), YF(:), 40, Q_parallel_full_lm(:), 'filled');
    axis equal;
    xlabel('x [m]');
    ylabel('y [m]');
    title('Filled full source Q_{\parallel,lm} in x-y');
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
% THETA–Z PLOT OF WINDOW HEAT FLUX (TRUE CYLINDRICAL MAP)
% ================================================================

S = load('window_bundle_qWindow_fromIntersections.mat');

% Required variables
theta_cent_deg = double(S.theta_cent_deg(:)).';
q_window_lm    = double(S.q_window_lm);

if ~isfield(S,'z_window_lm') || ~isfield(S,'theta_window_deg_lm')
    error('Need z_window_lm and theta_window_deg_lm for proper (theta,z) plot.');
end

z_window_lm = double(S.z_window_lm);
theta_hit_lm = double(S.theta_window_deg_lm);

[nL,nM] = size(q_window_lm);

% Flatten into scattered data
theta_scatter = theta_hit_lm(:);
z_scatter     = z_window_lm(:);
q_scatter     = q_window_lm(:);

valid = isfinite(theta_scatter) & isfinite(z_scatter) & isfinite(q_scatter);

theta_scatter = theta_scatter(valid);
z_scatter     = z_scatter(valid);
q_scatter     = q_scatter(valid);

% Handle periodic theta (important!)
theta_ext = [theta_scatter; theta_scatter-360; theta_scatter+360];
z_ext     = [z_scatter; z_scatter; z_scatter];
q_ext     = [q_scatter; q_scatter; q_scatter];

% Build regular grid
theta_grid = linspace(0,360,360);
z_grid     = linspace(min(z_ext), max(z_ext), 200);

[THG,ZG] = meshgrid(theta_grid, z_grid);

% Interpolate
q_theta_z = griddata(theta_ext, z_ext, q_ext, THG, ZG, 'linear');

% Fill holes
bad = ~isfinite(q_theta_z);
if any(bad(:))
    q_theta_z(bad) = griddata(theta_ext, z_ext, q_ext, THG(bad), ZG(bad), 'nearest');
end

%% ---------------- PLOT ----------------
figure('Color','w','Position',[100 100 900 400]);

imagesc(theta_grid, z_grid, q_theta_z);
set(gca,'YDir','normal');

xlabel('\theta [deg]');
ylabel('window hit z [m]');
title('q^W(\theta,z) from flux-tube mapping');

xlim([0 360]);
ylim([min(z_grid) max(z_grid)]);

clim([0 max(q_theta_z(:),[],'omitnan')]); % start at 0
colorbar;

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
    'fillMode', ...
    'centerValueMode', ...
    'centerFraction', ...
    'centerUserValue', ...
    'anchorMeasuredRow', ...
    'preserveThetaPattern', ...
    'P_parallel_full', ...
    'z_source', ...
    'q_anchor_raw_mean', ...
    'q_anchor_fill_mean', ...
    'q_center_mean', ...
    'nEdgeRowsForFill', ...
    'z_window_min_keep', ...
    'z_window_max_keep', ...
    'useWindowZMask', ...
    'valid_window_bin_mask', ...
    'minRowCoverageFraction', ...
    'rowCoverageFraction');

fprintf('Saved %s\n', outFile);