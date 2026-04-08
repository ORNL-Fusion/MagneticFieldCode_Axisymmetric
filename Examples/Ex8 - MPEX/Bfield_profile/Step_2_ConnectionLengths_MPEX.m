function Lout = Step_2_ConnectionLengths_MPEX(boundaryFile, fluxMapFile)
% Step_2_ConnectionLengths_MPEX
% Separate post-processing script to trace field lines and compute
% connection lengths for MPEX using the same vessel geometry and magnetic
% field workflow as the flux-mapping scripts.
%
% Usage
%   Lout = Step_2_ConnectionLengths_MPEX;
%   Lout = Step_2_ConnectionLengths_MPEX('MPEX_innerbound.xlsx');
%   Lout = Step_2_ConnectionLengths_MPEX('MPEX_innerbound.xlsx', 'bfield_MPEX_scenario_14.nc');
%
% Notes
% - This script traces field lines launched from the helicon window plane.
% - It computes forward, backward, and total connection length.
% - It uses the stepped vessel boundary carefully, without assuming unique z.
% - If a NetCDF field file exists, it uses that. Otherwise it rebuilds the
%   field from the coil setup and current scenario.

clearvars -except boundaryFile fluxMapFile
close all
clc

%% ------------------------- USER SETTINGS ---------------------------------
if nargin < 1 || isempty(boundaryFile)
    boundaryFile = 'MPEX_innerbound.xlsx';
end
if nargin < 2
    fluxMapFile = '';
end

confType      = 'conf_1';
scenarioIndex = 14;

% Computational domain [m]
zMin = -3.3;
zMax =  8.2;
rMax =  0.42;
Nz   = 900;
Nr   = 320;

% Launch radii normalized to the local window radius
xiLaunch = [0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.98];

% Field-line step controls
ds       = 0.002;   % m per integration step
maxSteps = 60000;   % maximum steps per trace

% Plot settings
doPlotGeometry = true;
plotMirror     = true;
%% ------------------------------------------------------------------------

%% SECTION 1: Get field either from NetCDF or from coil/current files
useNetCDF = false;
if ~isempty(fluxMapFile) && isfile(fluxMapFile)
    useNetCDF = true;
elseif isempty(fluxMapFile)
    fluxMapFile = sprintf('bfield_MPEX_scenario_%02d.nc', scenarioIndex);
    if isfile(fluxMapFile)
        useNetCDF = true;
    end
end

if useNetCDF
    fprintf('Reading magnetic field from %s\n', fluxMapFile);

    r1D  = ncread(fluxMapFile, 'r');
    z1D  = ncread(fluxMapFile, 'z');
    Br2D = double(ncread(fluxMapFile, 'br'));
    Bz2D = double(ncread(fluxMapFile, 'bz'));

    % Stored as [r,z] in file. Convert to [z,r] for consistency here.
    if isequal(size(Br2D), [numel(r1D), numel(z1D)])
        Br2D = permute(Br2D, [2 1]);
    end
    if isequal(size(Bz2D), [numel(r1D), numel(z1D)])
        Bz2D = permute(Bz2D, [2 1]);
    end

    coil = {};
else
    fprintf('Rebuilding magnetic field from coil setup and current scenario.\n');

    coilSetup = readtable('CoilSetup_MPEX.xlsx', 'Sheet', confType);
    coilSetup.z       = coilSetup.z       / 1000;
    coilSetup.dz      = coilSetup.dz      / 1000;
    coilSetup.r_inner = coilSetup.r_inner / 1000;
    coilSetup.r_outer = coilSetup.r_outer / 1000;

    currents = readmatrix('MPEX Coil Parameters - 2021-04-15.xlsx', ...
        'Sheet', 'MPEX Preliminary Design - Field', 'Range', 'K5:AH27');

    if scenarioIndex > size(currents,2)
        error('scenarioIndex exceeds the number of current scenarios in the spreadsheet.');
    end

    Iscenario = currents(:,scenarioIndex);
    if numel(Iscenario) ~= height(coilSetup)
        error('Number of currents does not match the number of coils in CoilSetup_MPEX.xlsx.');
    end

    coilCurrents = struct();
    for cc = 1:height(coilSetup)
        fieldName = sprintf('C%d', cc);
        coilCurrents.(fieldName) = Iscenario(cc);
    end

    coil = CreateCoilStructure(coilSetup, coilCurrents);

    r1D = linspace(1e-5, rMax, Nr);
    z1D = linspace(zMin,  zMax, Nz);
    [Br2D, Bz2D, ~, ~, ~, ~] = CalculateMagField(coil, z1D, r1D, 'grid');
end

%% SECTION 2: Read vessel, window, and target
[vessel, window, target, ~] = Step_0a_DrawVacuumVessel_MPEX_final(boundaryFile, false);

% Force helicon window center at z = 0 and use local inner-wall radius
window.z = 0.0;
window.r = localRadiusAtZStep(vessel.z, vessel.r, window.z);
window.label = 'helicon window center';

fprintf('Window reference: z = %.4f m, r = %.4f m\n', window.z, window.r);
fprintf('Target plane    : z = %.4f m, r = %.4f m\n', target.z, target.r);

%% SECTION 3: Launch points at the window plane
rStart = xiLaunch(:) * window.r;
zStart = window.z * ones(size(rStart));
rStart = min(rStart, 0.995 * window.r);

%% SECTION 4: Trace field lines and compute connection lengths
Lforward  = nan(size(rStart));
Lbackward = nan(size(rStart));
Ltotal    = nan(size(rStart));
hitF      = strings(size(rStart));
hitB      = strings(size(rStart));

zLineF = cell(size(rStart));
rLineF = cell(size(rStart));
zLineB = cell(size(rStart));
rLineB = cell(size(rStart));

for ii = 1:numel(rStart)
    [zF, rF, LF, hitTypeF] = localTraceConnectionLength( ...
        zStart(ii), rStart(ii), +1, z1D, r1D, Br2D, Bz2D, vessel, target, ds, maxSteps);

    [zB, rB, LB, hitTypeB] = localTraceConnectionLength( ...
        zStart(ii), rStart(ii), -1, z1D, r1D, Br2D, Bz2D, vessel, target, ds, maxSteps);

    zLineF{ii} = zF;
    rLineF{ii} = rF;
    zLineB{ii} = zB;
    rLineB{ii} = rB;

    Lforward(ii)  = LF;
    Lbackward(ii) = LB;
    Ltotal(ii)    = LF + LB;
    hitF(ii)      = string(hitTypeF);
    hitB(ii)      = string(hitTypeB);

    fprintf('Launch %2d: z0 = %.4f m, r0 = %.4f m | L+ = %.3f m (%s), L- = %.3f m (%s), Ltot = %.3f m\n', ...
        ii, zStart(ii), rStart(ii), LF, hitTypeF, LB, hitTypeB, LF+LB);
end

%% SECTION 5: Plot geometry and traced connection lengths
if doPlotGeometry
    figure('Color', 'w', 'Name', 'MPEX connection lengths');
    ax = axes;
    hold(ax, 'on')

    % ---------------- Vessel ----------------
    hVessel = plot(ax, vessel.z, vessel.r, 'k-', 'LineWidth', 1.5);
    if plotMirror
        plot(ax, vessel.z, -vessel.r, 'k-', 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end

    % ---------------- Coils ----------------
    hCoil = [];
    if ~isempty(coil)
        hCoil = plot(ax, coil{1}.zfil, coil{1}.rfil, 'r.');
        if plotMirror
            plot(ax, coil{1}.zfil, -coil{1}.rfil, 'r.', ...
                'HandleVisibility', 'off');
        end
        for jj = 2:numel(coil)
            plot(ax, coil{jj}.zfil, coil{jj}.rfil, 'r.', ...
                'HandleVisibility', 'off');
            if plotMirror
                plot(ax, coil{jj}.zfil, -coil{jj}.rfil, 'r.', ...
                    'HandleVisibility', 'off');
            end
        end
    end

    % ---------------- Window and target ----------------
    hWindow = plot(ax, window.z, window.r, 'bo', ...
        'MarkerFaceColor', 'b');
    if plotMirror
        plot(ax, window.z, -window.r, 'bo', ...
            'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
    end
    plot(ax, [window.z window.z], [-window.r window.r], 'b--', ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');

    hTarget = plot(ax, [target.z target.z], [-target.r target.r], 'k-', ...
        'LineWidth', 3);

    % ---------------- Colormap by connection length ----------------
    cmap = turbo(256);
    Lmin = min(Ltotal(isfinite(Ltotal)));
    Lmax = max(Ltotal(isfinite(Ltotal)));
    if isempty(Lmin) || isempty(Lmax) || ~isfinite(Lmin) || ~isfinite(Lmax)
        Lmin = 0;
        Lmax = 1;
    end
    if abs(Lmax - Lmin) < eps
        Lmax = Lmin + 1;
    end

    hStart = [];
    hForwardLegend = [];
    hBackwardLegend = [];

    for ii = 1:numel(rStart)

        % Normalize connection length and assign color
        Lnorm = (Ltotal(ii) - Lmin) / (Lmax - Lmin);
        Lnorm = max(0, min(1, Lnorm));
        cidx = max(1, min(size(cmap,1), 1 + round(Lnorm * (size(cmap,1)-1))));
        color_i = cmap(cidx, :);

        % Forward trace: solid, colored by total connection length
        hf = plot(ax, zLineF{ii}, rLineF{ii}, '-', ...
            'Color', color_i, 'LineWidth', 2.0);

        % Backward trace: dashed, same color
        hb = plot(ax, zLineB{ii}, rLineB{ii}, '--', ...
            'Color', color_i, 'LineWidth', 1.3);

        if plotMirror
            plot(ax, zLineF{ii}, -rLineF{ii}, '-', ...
                'Color', color_i, 'LineWidth', 2.0, ...
                'HandleVisibility', 'off');
            plot(ax, zLineB{ii}, -rLineB{ii}, '--', ...
                'Color', color_i, 'LineWidth', 1.3, ...
                'HandleVisibility', 'off');
        end

        % Launch point
        hs = plot(ax, zStart(ii), rStart(ii), 'ko', ...
            'MarkerFaceColor', 'y', 'MarkerSize', 5);

        if plotMirror
            plot(ax, zStart(ii), -rStart(ii), 'ko', ...
                'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
        end

        % Keep only one representative handle for legend
        if ii == 1
            hForwardLegend = hf;
            hBackwardLegend = hb;
            hStart = hs;
        else
            set(hf, 'HandleVisibility', 'off');
            set(hb, 'HandleVisibility', 'off');
            set(hs, 'HandleVisibility', 'off');
        end

        % Annotate with total connection length
        txt = sprintf('%.1f m', Ltotal(ii));
        text(ax, zStart(ii) + 0.03, rStart(ii), txt, ...
            'FontSize', 8, 'Color', color_i);

        if plotMirror
            text(ax, zStart(ii) + 0.03, -rStart(ii), txt, ...
                'FontSize', 8, 'Color', color_i);
        end
    end

    % ---------------- Axes formatting ----------------
    box(ax, 'on')
    grid(ax, 'on')
    xlabel(ax, 'z [m]')
    ylabel(ax, 'r [m]')
    title(ax, sprintf('MPEX connection lengths, scenario %d', scenarioIndex))
    xlim(ax, [min(z1D) max(z1D)])
    ylim(ax, [-max(r1D) max(r1D)])
    pbaspect(ax, [(max(z1D)-min(z1D)) (2*max(r1D)) 1])

    % ---------------- Colorbar ----------------
    colormap(ax, cmap)
    clim(ax, [Lmin Lmax])
    cb = colorbar(ax);
    cb.Label.String = 'Total connection length [m]';

    % ---------------- Legend ----------------
    legHandles = [hVessel hWindow hTarget hStart hForwardLegend hBackwardLegend];
    legLabels  = {'vessel', 'window', 'target', 'launch points', ...
                  'forward trace', 'backward trace'};

    if ~isempty(hCoil)
        legHandles = [hCoil(1) legHandles];
        legLabels  = [{'coil filaments'} legLabels];
    end

    legend(ax, legHandles, legLabels, 'Location', 'eastoutside');
end

%% SECTION 6: Summary plot of connection length versus launch radius
figure('Color', 'w', 'Name', 'Connection length summary');
plot(rStart, Lforward,  'o-', 'LineWidth', 1.5, 'DisplayName', 'forward')
hold on
plot(rStart, Lbackward, 's-', 'LineWidth', 1.5, 'DisplayName', 'backward')
plot(rStart, Ltotal,    'd-', 'LineWidth', 1.8, 'DisplayName', 'total')
grid on
box on
xlabel('Launch radius at window [m]')
ylabel('Connection length [m]')
title(sprintf('MPEX connection lengths, scenario %d', scenarioIndex))
legend('Location', 'best')

%% SECTION 7: Output structure
Lout = struct();
Lout.boundaryFile   = boundaryFile;
Lout.fluxMapFile    = fluxMapFile;
Lout.scenarioIndex  = scenarioIndex;
Lout.window         = window;
Lout.target         = target;
Lout.rStart         = rStart;
Lout.zStart         = zStart;
Lout.Lforward       = Lforward;
Lout.Lbackward      = Lbackward;
Lout.Ltotal         = Ltotal;
Lout.hitForward     = hitF;
Lout.hitBackward    = hitB;
Lout.zLineForward   = zLineF;
Lout.rLineForward   = rLineF;
Lout.zLineBackward  = zLineB;
Lout.rLineBackward  = rLineB;

disp('Done.')

end

%% ========================== LOCAL FUNCTIONS ==============================
function [zLine, rLine, Lconn, hitType] = localTraceConnectionLength( ...
    z0, r0, direction, z1D, r1D, Br2D, Bz2D, vessel, target, ds, maxSteps)

zLine = nan(maxSteps,1);
rLine = nan(maxSteps,1);
zLine(1) = z0;
rLine(1) = r0;

Lconn = 0;
hitType = "none";

for k = 1:maxSteps-1
    zc = zLine(k);
    rc = rLine(k);

    Br = interp2(z1D, r1D, Br2D', zc, rc, 'linear', NaN);
    Bz = interp2(z1D, r1D, Bz2D', zc, rc, 'linear', NaN);

    if ~isfinite(Br) || ~isfinite(Bz)
        hitType = "out_of_grid";
        zLine = zLine(1:k);
        rLine = rLine(1:k);
        return
    end

    Bmag = hypot(Br, Bz);
    if Bmag <= 0 || ~isfinite(Bmag)
        hitType = "zero_field";
        zLine = zLine(1:k);
        rLine = rLine(1:k);
        return
    end

    dz = direction * ds * (Bz / Bmag);
    dr = direction * ds * (Br / Bmag);

    zNext = zc + dz;
    rNext = rc + dr;

    if ((zc - target.z) * (zNext - target.z) <= 0) && abs(rNext) <= target.r
        Lconn = Lconn + hypot(zNext-zc, rNext-rc);
        zLine(k+1) = zNext;
        rLine(k+1) = rNext;
        hitType = "target";
        zLine = zLine(1:k+1);
        rLine = rLine(1:k+1);
        return
    end

    if zNext < min(z1D) || zNext > max(z1D) || rNext < 0 || rNext > max(r1D)
        Lconn = Lconn + hypot(zNext-zc, rNext-rc);
        zLine(k+1) = zNext;
        rLine(k+1) = rNext;
        hitType = "domain_exit";
        zLine = zLine(1:k+1);
        rLine = rLine(1:k+1);
        return
    end

    rw = localRadiusAtZStep(vessel.z, vessel.r, zNext);
    if ~isfinite(rw)
        Lconn = Lconn + hypot(zNext-zc, rNext-rc);
        zLine(k+1) = zNext;
        rLine(k+1) = rNext;
        hitType = "wall_unknown";
        zLine = zLine(1:k+1);
        rLine = rLine(1:k+1);
        return
    end

    if abs(rNext) >= rw
    % Approximate sub-step intersection with wall
    rAbs0 = abs(rc);
    rAbs1 = abs(rNext);
    denom = (rAbs1 - rAbs0);

    if abs(denom) > 1e-12
        alpha = (rw - rAbs0) / denom;
        alpha = max(0, min(1, alpha));
    else
        alpha = 1.0;
    end

    zHit = zc + alpha * (zNext - zc);
    rHit = rc + alpha * (rNext - rc);

    Lconn = Lconn + hypot(zHit-zc, rHit-rc);
    zLine(k+1) = zHit;
    rLine(k+1) = rHit;
    hitType = "wall";
    zLine = zLine(1:k+1);
    rLine = rLine(1:k+1);
    return
end

    zLine(k+1) = zNext;
    rLine(k+1) = rNext;
    Lconn = Lconn + hypot(zNext-zc, rNext-rc);
end

zLine = zLine(~isnan(zLine));
rLine = rLine(~isnan(rLine));
hitType = "max_steps";
end

function rw = localRadiusAtZStep(zv, rv, z0)
tol = 1e-10;
idx = find(abs(zv - z0) < tol);

if ~isempty(idx)
    rw = min(rv(idx));
    return
end

idxL = find(zv < z0, 1, 'last');
idxR = find(zv > z0, 1, 'first');

if isempty(idxL) && isempty(idxR)
    rw = NaN;
elseif isempty(idxL)
    rw = rv(idxR);
elseif isempty(idxR)
    rw = rv(idxL);
else
    rw = min(rv([idxL idxR]));
end
end
