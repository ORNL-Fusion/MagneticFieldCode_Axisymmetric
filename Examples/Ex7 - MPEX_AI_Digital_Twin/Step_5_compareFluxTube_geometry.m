% Compare flux-tube intersections on a fixed diagnostic cylinder for two B-field cases
%
% PURPOSE
%   Geometry-only comparison:
%   - Keep the same fixed source annulus (same tube identities l,m)
%   - Trace tubes in two different B-field cases
%   - Find where each tube hits the diagnostic cylinder
%   - Compare hit locations only (no heat-flux mapping)
%
% DIAGNOSTIC CYLINDER
%   r = R_cyl
%   z in [z_cyl_min, z_cyl_max]
%
% OUTPUTS
%   - z_hit(l,m) for old and new fields
%   - theta_hit(l,m) for old and new fields
%   - delta z and delta theta
%   - hit masks
%   - optional 3D trajectories with cylinder
%
% INPUT
%   1) fullSource_Qparallel.mat
%   2) reference B-field .nc
%   3) new B-field .nc
%
% OUTPUT
%   compareFluxTubeHitsOnCylinder_geometryOnly_output.nc
%
% -------------------------------------------------------------------------

clear; close all; clc;

%% ------------------------------------------------------------------------
% USER INPUTS
% -------------------------------------------------------------------------
sourceFile = 'fullSource_Qparallel.mat';

refBFile = 'bfield_protoMPEX_shotSeries_1.nc';
newBFile = 'bfield_protoMPEX_shotSeries_5.nc';

outFile = 'compareFluxTubeHitsOnCylinder_geometryOnly_output.nc';

% Diagnostic cylinder
R_cyl = 0.0626;        % [m]
z_cyl_min = 1.2;       % [m]
z_cyl_max = 3.0;       % [m]

% Tracing
ds_trace = 0.001;      % [m]
nStepsMax = 9000;
stop_if_outside_grid = true;
try_reverse_if_forward_fails = true;

% Plot controls
makePlots = true;
makeTrajectoryPlot = true;
nTrajToPlot = 40;

%% ------------------------------------------------------------------------
% STEP 0. LOAD FIXED SOURCE ANNULUS
% -------------------------------------------------------------------------
fprintf('Reading fixed source file: %s ...\n', sourceFile);

S = load(sourceFile);

requiredVars = {'Q_parallel_full_lm','r_cent_full','theta_cent_deg'};
for k = 1:numel(requiredVars)
    if ~isfield(S, requiredVars{k})
        error('Missing variable %s in %s.', requiredVars{k}, sourceFile);
    end
end

r_cent = double(S.r_cent_full(:));
theta_cent_deg = double(S.theta_cent_deg(:)).';

if isfield(S, 'z_source')
    z_source = double(S.z_source);
else
    z_source = 1.759;
end

Q_parallel_lm = double(S.Q_parallel_full_lm);
[nL, nM] = size(S.Q_parallel_full_lm);

fprintf('Loaded source annulus size = [%d x %d]\n', nL, nM);
fprintf('z_source = %.6f m\n', z_source);

[THETA_LM_DEG, R_LM] = meshgrid(theta_cent_deg, r_cent);

%% ------------------------------------------------------------------------
% STEP 1. TRACE REFERENCE CASE
% -------------------------------------------------------------------------
fprintf('Tracing reference case ...\n');

ref = trace_case_to_cylinder_geometryOnly( ...
    refBFile, R_cyl, z_cyl_min, z_cyl_max, ...
    R_LM, THETA_LM_DEG, z_source, ...
    ds_trace, nStepsMax, stop_if_outside_grid, try_reverse_if_forward_fails);

%% ------------------------------------------------------------------------
% STEP 2. TRACE NEW CASE
% -------------------------------------------------------------------------
fprintf('Tracing new case ...\n');

new = trace_case_to_cylinder_geometryOnly( ...
    newBFile, R_cyl, z_cyl_min, z_cyl_max, ...
    R_LM, THETA_LM_DEG, z_source, ...
    ds_trace, nStepsMax, stop_if_outside_grid, try_reverse_if_forward_fails);

%% ------------------------------------------------------------------------
% STEP 3. COMPARE GEOMETRY
% -------------------------------------------------------------------------
bothValid = ref.hit_in_range_lm & new.hit_in_range_lm;

dz_lm = NaN(nL, nM);
dtheta_lm = NaN(nL, nM);

dz_lm(bothValid) = new.z_hit_lm(bothValid) - ref.z_hit_lm(bothValid);

dtheta_lm(bothValid) = new.theta_hit_deg_lm(bothValid) - ref.theta_hit_deg_lm(bothValid);
dtheta_lm(bothValid) = mod(dtheta_lm(bothValid) + 180, 360) - 180;

%% ------------------------------------------------------------------------
% STEP 4. DIAGNOSTICS
% -------------------------------------------------------------------------
fprintf('\nREFERENCE CASE\n');
fprintf('Total in-range cylinder hits = %d / %d\n', nnz(ref.hit_in_range_lm), numel(ref.hit_in_range_lm));
fprintf('z_hit range = [%g, %g] m\n', ...
    min(ref.z_hit_lm(:), [], 'omitnan'), max(ref.z_hit_lm(:), [], 'omitnan'));

fprintf('\nNEW CASE\n');
fprintf('Total in-range cylinder hits = %d / %d\n', nnz(new.hit_in_range_lm), numel(new.hit_in_range_lm));
fprintf('z_hit range = [%g, %g] m\n', ...
    min(new.z_hit_lm(:), [], 'omitnan'), max(new.z_hit_lm(:), [], 'omitnan'));

fprintf('\nCOMPARISON\n');
fprintf('Valid in both cases = %d / %d\n', nnz(bothValid), numel(bothValid));
if any(bothValid(:))
    fprintf('Mean |dz| = %g m\n', mean(abs(dz_lm(bothValid)), 'omitnan'));
    fprintf('Mean |dtheta| = %g deg\n', mean(abs(dtheta_lm(bothValid)), 'omitnan'));
end

%% ------------------------------------------------------------------------
% STEP 5. PLOTS
% -------------------------------------------------------------------------
if makePlots
    % ------------------------------------------------------------
    % Figure 1: z_hit maps
    % ------------------------------------------------------------
    figure('Color','w','Position',[80 80 1250 900]);

    subplot(2,2,1);
    imagesc(theta_cent_deg, r_cent, ref.z_hit_lm);
    set(gca,'YDir','normal');
    xlabel('\theta source [deg]');
    ylabel('source annulus r [m]');
    title('Reference z_{hit} on cylinder');
    colorbar;

    subplot(2,2,2);
    imagesc(theta_cent_deg, r_cent, new.z_hit_lm);
    set(gca,'YDir','normal');
    xlabel('\theta source [deg]');
    ylabel('source annulus r [m]');
    title('New z_{hit} on cylinder');
    colorbar;

    subplot(2,2,3);
    imagesc(theta_cent_deg, r_cent, dz_lm);
    set(gca,'YDir','normal');
    xlabel('\theta source [deg]');
    ylabel('source annulus r [m]');
    title('\Delta z_{hit} = z_{new} - z_{ref}');
    colorbar;

    subplot(2,2,4);
    imagesc(theta_cent_deg, r_cent, bothValid);
    set(gca,'YDir','normal');
    xlabel('\theta source [deg]');
    ylabel('source annulus r [m]');
    title('Valid hits in both cases');
    colorbar;

    sgtitle('Flux-tube hit comparison on diagnostic cylinder');

    % ------------------------------------------------------------
    % Figure 2: theta_hit maps
    % ------------------------------------------------------------
    figure('Color','w','Position',[100 100 1250 900]);

    subplot(2,2,1);
    imagesc(theta_cent_deg, r_cent, ref.theta_hit_deg_lm);
    set(gca,'YDir','normal');
    xlabel('\theta source [deg]');
    ylabel('source annulus r [m]');
    title('Reference \theta_{hit} on cylinder');
    colorbar;

    subplot(2,2,2);
    imagesc(theta_cent_deg, r_cent, new.theta_hit_deg_lm);
    set(gca,'YDir','normal');
    xlabel('\theta source [deg]');
    ylabel('source annulus r [m]');
    title('New \theta_{hit} on cylinder');
    colorbar;

    subplot(2,2,3);
    imagesc(theta_cent_deg, r_cent, dtheta_lm);
    set(gca,'YDir','normal');
    xlabel('\theta source [deg]');
    ylabel('source annulus r [m]');
    title('\Delta \theta_{hit} [deg]');
    colorbar;

    subplot(2,2,4);
    imagesc(theta_cent_deg, r_cent, ref.hit_in_range_lm + 2*new.hit_in_range_lm);
    set(gca,'YDir','normal');
    xlabel('\theta source [deg]');
    ylabel('source annulus r [m]');
    title('Hit mask code: 0 none, 1 ref only, 2 new only, 3 both');
    colorbar;

        % ------------------------------------------------------------
    % Figure 3: theta-z scatter on cylinder
    % ------------------------------------------------------------
    figure('Color','w','Position',[120 120 1150 420]);

    Rcolor = repmat(r_cent, 1, nM);

    subplot(1,2,1);
    scatter(ref.theta_hit_deg_lm(ref.hit_in_range_lm), ...
            ref.z_hit_lm(ref.hit_in_range_lm), ...
            16, ...
            Rcolor(ref.hit_in_range_lm), ...
            'filled');
    xlabel('\theta_{hit} [deg]');
    ylabel('z_{hit} [m]');
    title('Reference hit locations on cylinder');
    xlim([0 360]);
    ylim([z_cyl_min z_cyl_max]);
    colorbar;
    grid on;

    subplot(1,2,2);
    scatter(new.theta_hit_deg_lm(new.hit_in_range_lm), ...
            new.z_hit_lm(new.hit_in_range_lm), ...
            16, ...
            Rcolor(new.hit_in_range_lm), ...
            'filled');
    xlabel('\theta_{hit} [deg]');
    ylabel('z_{hit} [m]');
    title('New hit locations on cylinder');
    xlim([0 360]);
    ylim([z_cyl_min z_cyl_max]);
    colorbar;
    grid on;

    % ------------------------------------------------------------
    % Figure 4: r-z view of cylinder hits (axisymmetric projection)
    % ------------------------------------------------------------
    figure('Color','w','Position',[140 140 1150 420]);

    subplot(1,2,1);
    hold on;
    plot([R_cyl R_cyl], [z_cyl_min z_cyl_max], 'k-', 'LineWidth', 2);
    scatter(R_cyl*ones(nnz(ref.hit_in_range_lm),1), ref.z_hit_lm(ref.hit_in_range_lm), ...
        20, ref.theta_hit_deg_lm(ref.hit_in_range_lm), 'filled');
    xlabel('r [m]');
    ylabel('z [m]');
    title('Reference hits on diagnostic cylinder (r-z view)');
    xlim([R_cyl-0.005, R_cyl+0.005]);
    ylim([z_cyl_min, z_cyl_max]);
    colorbar;
    grid on;

    subplot(1,2,2);
    hold on;
    plot([R_cyl R_cyl], [z_cyl_min z_cyl_max], 'k-', 'LineWidth', 2);
    scatter(R_cyl*ones(nnz(new.hit_in_range_lm),1), new.z_hit_lm(new.hit_in_range_lm), ...
        20, new.theta_hit_deg_lm(new.hit_in_range_lm), 'filled');
    xlabel('r [m]');
    ylabel('z [m]');
    title('New hits on diagnostic cylinder (r-z view)');
    xlim([R_cyl-0.005, R_cyl+0.005]);
    ylim([z_cyl_min, z_cyl_max]);
    colorbar;
    grid on;

    % ------------------------------------------------------------
    % Figure 5: optional 3D trajectories
    % ------------------------------------------------------------
    if makeTrajectoryPlot
        figure('Color','w','Position',[160 160 1200 500]);

        subplot(1,2,1);
        hold on;
        plot_trajectories_with_cylinder(ref.x_path, ref.y_path, ref.z_path, ...
            R_cyl, z_cyl_min, z_cyl_max, nTrajToPlot);
        title('Reference trajectories');
        xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
        axis equal; grid on; view(3);

        subplot(1,2,2);
        hold on;
        plot_trajectories_with_cylinder(new.x_path, new.y_path, new.z_path, ...
            R_cyl, z_cyl_min, z_cyl_max, nTrajToPlot);
        title('New trajectories');
        xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
        axis equal; grid on; view(3);
    end
end

%% ------------------------------------------------------------------------
% STEP 6. SAVE OUTPUT
% -------------------------------------------------------------------------
fprintf('Writing %s ...\n', outFile);

if exist(outFile, 'file')
    delete(outFile);
end

nccreate(outFile, 'theta_cent_deg', 'Dimensions', {'m', nM});
nccreate(outFile, 'r_cent', 'Dimensions', {'l', nL});

nccreate(outFile, 'z_hit_ref_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'theta_hit_ref_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'z_hit_new_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'theta_hit_new_lm', 'Dimensions', {'l', nL, 'm', nM});

nccreate(outFile, 'dz_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'dtheta_lm', 'Dimensions', {'l', nL, 'm', nM});

nccreate(outFile, 'hit_in_range_ref_lm', 'Dimensions', {'l', nL, 'm', nM});
nccreate(outFile, 'hit_in_range_new_lm', 'Dimensions', {'l', nL, 'm', nM});

ncwrite(outFile, 'theta_cent_deg', theta_cent_deg);
ncwrite(outFile, 'r_cent', r_cent);

ncwrite(outFile, 'z_hit_ref_lm', ref.z_hit_lm);
ncwrite(outFile, 'theta_hit_ref_lm', ref.theta_hit_deg_lm);
ncwrite(outFile, 'z_hit_new_lm', new.z_hit_lm);
ncwrite(outFile, 'theta_hit_new_lm', new.theta_hit_deg_lm);

ncwrite(outFile, 'dz_lm', dz_lm);
ncwrite(outFile, 'dtheta_lm', dtheta_lm);

ncwrite(outFile, 'hit_in_range_ref_lm', double(ref.hit_in_range_lm));
ncwrite(outFile, 'hit_in_range_new_lm', double(new.hit_in_range_lm));

fprintf('Done.\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================
function out = trace_case_to_cylinder_geometryOnly( ...
    bfile, R_cyl, z_cyl_min, z_cyl_max, ...
    R_LM, THETA_LM_DEG, z_source, ...
    ds_trace, nStepsMax, stop_if_outside_grid, try_reverse_if_forward_fails)

    fprintf('  Reading B-field: %s ...\n', bfile);

    r = double(ncread(bfile, 'r')); r = r(:);
    z = double(ncread(bfile, 'z')); z = z(:);

    Br = double(ncread(bfile, 'br'));
    Bt = double(ncread(bfile, 'bt'));
    Bz = double(ncread(bfile, 'bz'));

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

    [nL, nM] = size(R_LM);

    theta_hit_deg_lm = NaN(nL, nM);
    z_hit_lm = NaN(nL, nM);
    hit_in_range_lm = false(nL, nM);

    x_path = NaN(nL, nM, nStepsMax);
    y_path = NaN(nL, nM, nStepsMax);
    z_path = NaN(nL, nM, nStepsMax);

    for il = 1:nL
        for im = 1:nM
            r0 = R_LM(il,im);
            th0 = THETA_LM_DEG(il,im);

            [ok, th_hit, z_hit, xpath, ypath, zpath_] = trace_point_to_cylinder( ...
                r0, th0, z_source, ...
                FBr, FBt, FBz, r, z, R_cyl, ds_trace, nStepsMax, stop_if_outside_grid, +1);

            if ~ok && try_reverse_if_forward_fails
                [ok, th_hit, z_hit, xpath, ypath, zpath_] = trace_point_to_cylinder( ...
                    r0, th0, z_source, ...
                    FBr, FBt, FBz, r, z, R_cyl, ds_trace, nStepsMax, stop_if_outside_grid, -1);
            end

            if ok
                ns = min(numel(xpath), nStepsMax);
                x_path(il,im,1:ns) = xpath(1:ns);
                y_path(il,im,1:ns) = ypath(1:ns);
                z_path(il,im,1:ns) = zpath_(1:ns);

                if z_hit >= z_cyl_min && z_hit <= z_cyl_max
                    theta_hit_deg_lm(il,im) = mod(th_hit, 360);
                    z_hit_lm(il,im) = z_hit;
                    hit_in_range_lm(il,im) = true;
                end
            end
        end
    end

    out.theta_hit_deg_lm = theta_hit_deg_lm;
    out.z_hit_lm = z_hit_lm;
    out.hit_in_range_lm = hit_in_range_lm;
    out.x_path = x_path;
    out.y_path = y_path;
    out.z_path = z_path;
end

function [ok, th_hit, z_hit, xpath, ypath, zpath_] = trace_point_to_cylinder( ...
    r0, th0_deg, z0, ...
    FBr, FBt, FBz, r_grid, z_grid, R_cyl, ds_trace, ...
    nStepsMax, stop_if_outside_grid, trace_sign)

    ok = false;
    th_hit = NaN;
    z_hit = NaN;

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

        crossed = ((r_now - R_cyl) * (r_next - R_cyl) <= 0) && abs(r_next - r_now) > 1e-12;

        if crossed
            f = (R_cyl - r_now) / (r_next - r_now);
            if f >= 0 && f <= 1
                th_hit = mod(th_now_deg + f * (th_next_deg - th_now_deg), 360);
                z_hit = z_now + f * (z_next - z_now);

                xhit = R_cyl * cosd(th_hit);
                yhit = R_cyl * sind(th_hit);

                xpath(is) = xhit;
                ypath(is) = yhit;
                zpath_(is) = z_hit;
                ok = true;
                return;
            end
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

function plot_trajectories_with_cylinder(x_path, y_path, z_path, R_cyl, zmin_c, zmax_c, nTrajToPlot)
    [nL, nM, ~] = size(x_path);

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

    th_c = linspace(0, 2*pi, 100);
    z_c = [zmin_c zmax_c];
    [THC, ZC] = meshgrid(th_c, z_c);
    XC = R_cyl * cos(THC);
    YC = R_cyl * sin(THC);
    mesh(XC, YC, ZC, 'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0);
end