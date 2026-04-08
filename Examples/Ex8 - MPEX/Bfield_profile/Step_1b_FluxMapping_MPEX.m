% Step_1b_FluxMapping_MPEX
% MPEX flux mapping using either:
%   (1) wall-limited reference (last uninterrupted / last wall-avoiding surface),
%   (2) window/source referenced flux tubes,
%   (3) or both.
%
% Reference choice:
%   fluxReferenceMode = 'wall'   -> only plot the last uninterrupted wall-limited contour
%   fluxReferenceMode = 'window' -> only plot source/window-referenced flux tubes
%   fluxReferenceMode = 'both'   -> plot both

clearvars
close all
clc

%% ------------------------- USER SETTINGS ---------------------------------
confType      = 'conf_1';
scenarioIndex = 14;
makeNetCDF    = true;
saveFigure    = false;

% Choose what to plot: 'wall', 'window', or 'both'
fluxReferenceMode = 'wall';

% Computational domain [m]
zMin = -3.3;
zMax =  8.2;
rMax =  0.42;
Nz   = 900;
Nr   = 320;

% Window/source referenced contour levels
xi_lines = [0.15 0.30 0.45 0.60 0.75 0.90 1.00];
%% ------------------------------------------------------------------------

showWindowFlux = any(strcmpi(fluxReferenceMode, {'window','both'}));
showWallFlux   = any(strcmpi(fluxReferenceMode, {'wall','both'}));

%% SECTION 1: Read MPEX coil setup
coilSetup = readtable('CoilSetup_MPEX.xlsx', 'Sheet', confType);
coilSetup.z       = coilSetup.z       / 1000;  % mm -> m
coilSetup.dz      = coilSetup.dz      / 1000;
coilSetup.r_inner = coilSetup.r_inner / 1000;
coilSetup.r_outer = coilSetup.r_outer / 1000;

%% SECTION 2: Read one current scenario from the design spreadsheet
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

%% SECTION 3: Create coil structure and evaluate magnetic field
coil = CreateCoilStructure(coilSetup, coilCurrents);

r1D = linspace(1e-5, rMax, Nr);
z1D = linspace(zMin,  zMax, Nz);
[Br2D, Bz2D, Atheta2D, phi2D, z2D, r2D] = CalculateMagField(coil, z1D, r1D, 'grid'); 
B2D = sqrt(Br2D.^2 + Bz2D.^2);

%% SECTION 4: Read vacuum vessel / source window / target geometry carefully
[vessel, window, target, vessel_raw] = Step_0a_DrawVacuumVessel_MPEX('MPEX_innerbound.xlsx', false);

%% SECTION 5: Wall-limited reference flux (LCFS-like boundary)
phiWall = interp2(z1D, r1D, phi2D', vessel.z, vessel.r, 'linear', NaN);
valid   = isfinite(phiWall);
phiWall(~valid) = inf;
[phiWallMin, idxWall] = min(phiWall);

limitPoint.z = vessel.z(idxWall);
limitPoint.r = vessel.r(idxWall);

xi_wall = phi2D / phiWallMin;
Cwall = contourc(z1D, r1D, xi_wall', [1 1]);
[zWallFlux, rWallFlux] = localParseContourMatrix(Cwall);

%% SECTION 6: Window-referenced flux-tube mapping
zFlux = {};
rFlux = {};
phiWindowEdge = NaN;
xi_window = nan(size(phi2D));
if showWindowFlux
    phiWindowEdge = interp2(z1D, r1D, phi2D', window.z, window.r, 'linear');
    if ~isfinite(phiWindowEdge)
        error(['Window reference point lies outside the field grid. ' ...
               'Increase domain or adjust window.z/window.r.']);
    end
    xi_window = phi2D / phiWindowEdge;

    zFlux = cell(1, numel(xi_lines));
    rFlux = cell(1, numel(xi_lines));
    for jj = 1:numel(xi_lines)
        C = contourc(z1D, r1D, xi_window', [xi_lines(jj) xi_lines(jj)]);
        [zFlux{jj}, rFlux{jj}] = localParseContourMatrix(C);
    end
end

%% SECTION 7: Plot flux mapping + B field
figure('Color', 'w', 'Name', 'MPEX flux mapping and B field');

% ---------------- Top subplot: geometry + flux tubes ----------------
ax1 = subplot(2,1,1);
hold(ax1, 'on')

% Coils
hCoil = plot(ax1, coil{1}.zfil, coil{1}.rfil, 'r.', 'DisplayName', 'coil filaments');
plot(ax1, coil{1}.zfil, -coil{1}.rfil, 'r.', 'HandleVisibility', 'off');
for ii = 2:numel(coil)
    plot(ax1, coil{ii}.zfil,  coil{ii}.rfil, 'r.', 'HandleVisibility', 'off');
    plot(ax1, coil{ii}.zfil, -coil{ii}.rfil, 'r.', 'HandleVisibility', 'off');
end

% Vessel
hVessel = plot(ax1, vessel.z,  vessel.r, 'k-', 'LineWidth', 1.5, 'DisplayName', 'vessel');
plot(ax1, vessel.z, -vessel.r, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Window marker
hWindow = plot(ax1, window.z,  window.r, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'window ref');
plot(ax1, window.z, -window.r, 'bo', 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
plot(ax1, [window.z window.z], [-window.r window.r], 'b--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Target plane
hTarget = plot(ax1, [target.z target.z], [-target.r target.r], 'k-', 'LineWidth', 3, 'DisplayName', 'target');

% Window-referenced flux tubes
hFlux = [];
gotFluxHandle = false;
if showWindowFlux
    for jj = 1:numel(zFlux)
        for kk = 1:numel(zFlux{jj})
            htmp = plot(ax1, zFlux{jj}{kk},  rFlux{jj}{kk}, 'b-', 'LineWidth', 1.0);
            plot(ax1, zFlux{jj}{kk}, -rFlux{jj}{kk}, 'b-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
            if ~gotFluxHandle
                hFlux = htmp;
                set(hFlux, 'DisplayName', 'window flux tubes');
                gotFluxHandle = true;
            else
                set(htmp, 'HandleVisibility', 'off');
            end
        end
    end
end

% Wall-limited contour = last uninterrupted wall-avoiding surface
hWall = [];
gotWallHandle = false;
if showWallFlux
    for kk = 1:numel(zWallFlux)
        htmp = plot(ax1, zWallFlux{kk},  rWallFlux{kk}, 'm-', 'LineWidth', 2.0);
        plot(ax1, zWallFlux{kk}, -rWallFlux{kk}, 'm-', 'LineWidth', 2.0, 'HandleVisibility', 'off');
        if ~gotWallHandle
            hWall = htmp;
            set(hWall, 'DisplayName', 'last uninterrupted surface');
            gotWallHandle = true;
        else
            set(htmp, 'HandleVisibility', 'off');
        end
    end
    plot(ax1, limitPoint.z,  limitPoint.r, 'mo', 'MarkerFaceColor', 'm', 'HandleVisibility', 'off');
    plot(ax1, limitPoint.z, -limitPoint.r, 'mo', 'MarkerFaceColor', 'm', 'HandleVisibility', 'off');
end

box(ax1, 'on')
grid(ax1, 'on')
xlabel(ax1, 'z [m]')
ylabel(ax1, 'r [m]')
title(ax1, sprintf('MPEX flux-tube mapping, scenario %d', scenarioIndex))
xlim(ax1, [zMin zMax])
ylim(ax1, [-rMax rMax])
pbaspect(ax1, [(zMax-zMin) (2*rMax) 1])

% Legend handles
legHandles = [hCoil(1) hVessel(1)];
legLabels  = {'coil filaments', 'vessel'};
if gotFluxHandle
    legHandles = [legHandles hFlux(1)];
    legLabels  = [legLabels {'window flux tubes'}];
end
legHandles = [legHandles hWindow(1) hTarget(1)];
legLabels  = [legLabels {'window ref', 'target'}];
if gotWallHandle
    legHandles = [legHandles hWall(1)];
    legLabels  = [legLabels {'last uninterrupted surface'}];
end
legend(ax1, legHandles, legLabels, 'Location', 'eastoutside');

% ---------------- Bottom subplot: on-axis B field ----------------
ax2 = subplot(2,1,2);
hold(ax2, 'on')
plot(ax2, z1D, B2D(:,1), 'k-', 'LineWidth', 2);

yl = [min(B2D(:,1)) max(B2D(:,1))];
if diff(yl) <= 0
    yl = [yl(1)-0.1 yl(2)+0.1];
else
    pad = 0.02 * diff(yl);
    yl = [yl(1)-pad yl(2)+pad];
end
ylim(ax2, yl)

plot(ax2, [window.z window.z], yl, 'b--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax2, [target.z target.z], yl, 'r--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
if gotWallHandle
    plot(ax2, [limitPoint.z limitPoint.z], yl, 'm--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

box(ax2, 'on')
grid(ax2, 'on')
xlabel(ax2, 'z [m]')
ylabel(ax2, '|B|(r = 0) [T]')
title(ax2, 'On-axis magnetic field magnitude')
xlim(ax2, [zMin zMax])

% Tight spacing and x alignment
linkaxes([ax1 ax2], 'x');
set(ax1, 'Position', [0.08 0.56 0.80 0.34]);
set(ax2, 'Position', [0.08 0.10 0.80 0.34]);

%%
%% SECTION 7: Plot flux mapping (single plot)
figure('Color', 'w', 'Name', 'MPEX flux mapping');

ax1 = axes;
hold(ax1, 'on')

% ---------------- Coils ----------------
hCoil = plot(ax1, coil{1}.zfil, coil{1}.rfil, 'r.', ...
    'DisplayName', 'coil filaments');
plot(ax1, coil{1}.zfil, -coil{1}.rfil, 'r.', 'HandleVisibility', 'off');

for ii = 2:numel(coil)
    plot(ax1, coil{ii}.zfil,  coil{ii}.rfil, 'r.', 'HandleVisibility', 'off');
    plot(ax1, coil{ii}.zfil, -coil{ii}.rfil, 'r.', 'HandleVisibility', 'off');
end

% ---------------- Vessel ----------------
hVessel = plot(ax1, vessel.z,  vessel.r, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'vessel');
plot(ax1, vessel.z, -vessel.r, 'k-', 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');

% ---------------- Window marker ----------------
hWindow = plot(ax1, window.z,  window.r, 'bo', ...
    'MarkerFaceColor', 'b', 'DisplayName', 'window ref');
plot(ax1, window.z, -window.r, 'bo', ...
    'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
plot(ax1, [window.z window.z], [-window.r window.r], 'b--', ...
    'LineWidth', 1.2, 'HandleVisibility', 'off');

% ---------------- Target plane ----------------
hTarget = plot(ax1, [target.z target.z], [-target.r target.r], 'k-', ...
    'LineWidth', 3, 'DisplayName', 'target');

% ---------------- Window-referenced flux tubes ----------------
hFlux = [];
gotFluxHandle = false;
if showWindowFlux
    for jj = 1:numel(zFlux)
        for kk = 1:numel(zFlux{jj})
            htmp = plot(ax1, zFlux{jj}{kk},  rFlux{jj}{kk}, 'b-', 'LineWidth', 1.0);
            plot(ax1, zFlux{jj}{kk}, -rFlux{jj}{kk}, 'b-', 'LineWidth', 1.0, ...
                'HandleVisibility', 'off');

            if ~gotFluxHandle
                hFlux = htmp;
                set(hFlux, 'DisplayName', 'window flux tubes');
                gotFluxHandle = true;
            else
                set(htmp, 'HandleVisibility', 'off');
            end
        end
    end
end

% ---------------- Wall-limited contour ----------------
hWall = [];
gotWallHandle = false;
if showWallFlux
    for kk = 1:numel(zWallFlux)
        htmp = plot(ax1, zWallFlux{kk},  rWallFlux{kk}, 'm-', 'LineWidth', 2.0);
        plot(ax1, zWallFlux{kk}, -rWallFlux{kk}, 'm-', 'LineWidth', 2.0, ...
            'HandleVisibility', 'off');

        if ~gotWallHandle
            hWall = htmp;
            set(hWall, 'DisplayName', 'last uninterrupted surface');
            gotWallHandle = true;
        else
            set(htmp, 'HandleVisibility', 'off');
        end
    end

    plot(ax1, limitPoint.z,  limitPoint.r, 'mo', 'MarkerFaceColor', 'm', ...
        'HandleVisibility', 'off');
    plot(ax1, limitPoint.z, -limitPoint.r, 'mo', 'MarkerFaceColor', 'm', ...
        'HandleVisibility', 'off');
end

% ---------------- Axes formatting ----------------
box(ax1, 'on')
grid(ax1, 'on')
xlabel(ax1, 'z [m]')
ylabel(ax1, 'r [m]')
title(ax1, sprintf('MPEX flux-tube mapping, scenario %d', scenarioIndex))
xlim(ax1, [zMin zMax])
ylim(ax1, [-rMax rMax])
pbaspect(ax1, [(zMax-zMin) (2*rMax) 1])

% ---------------- Legend ----------------
legHandles = [hCoil(1) hVessel(1)];
legLabels  = {'coil filaments', 'vessel'};

if gotFluxHandle
    legHandles = [legHandles hFlux(1)];
    legLabels  = [legLabels {'window flux tubes'}];
end

legHandles = [legHandles hWindow(1) hTarget(1)];
legLabels  = [legLabels {'window ref', 'target'}];

if gotWallHandle
    legHandles = [legHandles hWall(1)];
    legLabels  = [legLabels {'last uninterrupted surface'}];
end

legend(ax1, legHandles, legLabels, 'Location', 'eastoutside');

%% SECTION 8: Save field data for external tools
if makeNetCDF
    file = sprintf('bfield_MPEX_scenario_%02d.nc', scenarioIndex);
    if exist(file, 'file')
        delete(file);
    end

    r = r1D(:);
    z = z1D(:);
    Br = Br2D;
    Bz = Bz2D;
    Bt = zeros(size(Bz2D), 'like', Bz2D);

    nccreate(file, 'r', 'Dimensions', {'r', numel(r)}, 'Datatype', 'double');
    nccreate(file, 'z', 'Dimensions', {'z', numel(z)}, 'Datatype', 'double');
    ncwrite(file, 'r', r);
    ncwrite(file, 'z', z);

    nccreate(file, 'br', 'Dimensions', {'r', numel(r), 'z', numel(z)}, 'Datatype', 'single');
    nccreate(file, 'bt', 'Dimensions', {'r', numel(r), 'z', numel(z)}, 'Datatype', 'single');
    nccreate(file, 'bz', 'Dimensions', {'r', numel(r), 'z', numel(z)}, 'Datatype', 'single');
    nccreate(file, 'psi_wall',   'Dimensions', {'r', numel(r), 'z', numel(z)}, 'Datatype', 'single');
    nccreate(file, 'psi_window', 'Dimensions', {'r', numel(r), 'z', numel(z)}, 'Datatype', 'single');

    if isequal(size(Br), [numel(z) numel(r)]), Br = permute(Br, [2 1]); end
    if isequal(size(Bt), [numel(z) numel(r)]), Bt = permute(Bt, [2 1]); end
    if isequal(size(Bz), [numel(z) numel(r)]), Bz = permute(Bz, [2 1]); end

    psiWallOut = xi_wall;
    if isequal(size(psiWallOut), [numel(z) numel(r)]), psiWallOut = permute(psiWallOut, [2 1]); end

    psiWindowOut = xi_window;
    if isequal(size(psiWindowOut), [numel(z) numel(r)]), psiWindowOut = permute(psiWindowOut, [2 1]); end

    ncwrite(file, 'br', single(Br));
    ncwrite(file, 'bt', single(Bt));
    ncwrite(file, 'bz', single(Bz));
    ncwrite(file, 'psi_wall', single(psiWallOut));
    ncwrite(file, 'psi_window', single(psiWindowOut));
end

if saveFigure
    exportgraphics(gcf, sprintf('MPEX_FluxMapping_scenario_%02d.pdf', scenarioIndex), 'Resolution', 300)
end

disp('Done.')

%% -------------------------- LOCAL FUNCTION -------------------------------
function [zSegments, rSegments] = localParseContourMatrix(C)
zSegments = {};
rSegments = {};
if isempty(C)
    return
end

idx = 1;
seg = 0;
while idx < size(C,2)
    nPts = C(2,idx);
    seg  = seg + 1;
    cols = (idx+1):(idx+nPts);
    zSegments{seg} = C(1, cols);
    rSegments{seg} = C(2, cols);
    idx = idx + nPts + 1;
end
end
