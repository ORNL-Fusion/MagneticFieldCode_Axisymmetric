% OBJECTIVE:
% 1- Demonstrate the use of the functions: "CreateCoilStructure" and
%    "CalculateMagField" to calculate the magnetic field on Proto-MPEX
% 2- Demonstrate the use of the "CoilSetup" spreadsheet
% 3- Verify the numerical implementation against calculations from other codes
% 4- Export the magnetic field to NetCDF and CSV
%
% Written by J.F. Caneses Marin
% Updated for debugging/cleanup

clearvars
clc
close all

%% SECTION 1: Read "CoilSetup" spreadsheet
% =========================================================================
homeFolder = cd;
cd ../..
currentFolder = cd;
functionFolder = [currentFolder, '/Functions'];
addpath(genpath(functionFolder));
cd(homeFolder)

dum1 = tic;
disp('Reading "CoilSetup" spreadsheet...')
coilSetup = readtable('CoilSetup_ProtoMPEX.xlsx','Sheet','conf_A');
disp('Reading complete!')

% Small correction to the coil geometry
coilSetup.r_inner = coilSetup.r_inner - 0.5e-2;

% Assignment of currents per power supply
coilCurrents.TR1 = 600;
coilCurrents.TR2 = 160;
coilCurrents.PS1 = 4500;
coilCurrents.PS2 = 4500;

% Create coil structure
disp('Creating "coil" structure...')
coil = CreateCoilStructure(coilSetup, coilCurrents);
disp(['Complete! Elapsed time: ', num2str(toc(dum1)), ' s'])
clearvars dum1

disp(coilSetup)

%% SECTION 2: Calculate magnetic field
% =========================================================================
tic

% Use finer radial resolution to avoid blocky banding in imagesc
r1D = linspace(1e-3, 0.1, 100);
z1D = linspace(0.0, 5.0, 501);

dum1 = tic;
disp('Calculating magnetic field...')
[Br2D, Bz2D, Atheta2D, Phi2D, z2D, r2D] = CalculateMagField(coil, z1D, r1D, 'grid');
disp(['Complete! Elapsed time: ', num2str(toc(dum1)), ' s'])
clearvars dum1

% Magnetic field magnitude
B2D = sqrt(Br2D.^2 + Bz2D.^2);
toc

%% SECTION 3: Get data from Jeremy Lore's code
% =========================================================================
fileName = 'PS1_4500_PS2_4500_TR1_160_TR2_600.txt';
d = load(fileName);
z_Bjl = d(:,1);
Bjl   = d(:,2);

%% SECTION 4: Verification plots
% =========================================================================
close all

figure('Color','w')

% -------------------------------------------------------------------------
% Compare present calculation and Jeremy's code
subplot(7,1,[1:5])
hold on

% Near-axis profile: use the first radial location
n_r = 1;
B_onAxis = B2D(:, n_r);

% Make sure comparison arrays have compatible lengths
nCompare = min(numel(z1D), numel(z_Bjl));
nCompareB = min(numel(B_onAxis), numel(Bjl));
nCompareAll = min(nCompare, nCompareB);

hdum(1) = plot(z1D(1:nCompareAll), B_onAxis(1:nCompareAll), 'k', 'LineWidth', 3);
hdum(2) = plot(z_Bjl(1:nCompareAll), Bjl(1:nCompareAll), 'g.', 'LineWidth', 1);

for ii = 1:numel(coil)
    plot(coil{ii}.zfil, coil{ii}.rfil, 'r.');
end

set(gca,'FontName','times','FontSize',11)
ylabel('B [T]','Interpreter','latex','FontSize',13)
grid on
box on
xlim([z1D(1), z1D(end)])
hL = legend(hdum,'Present calc','JL');
set(hL,'Interpreter','latex','FontSize',11)
set(gca,'XTickLabel',[])

% -------------------------------------------------------------------------
% Plot error relative to max(B)
subplot(7,1,[6,7])
errorBfield = abs(B_onAxis(1:nCompareAll) - Bjl(1:nCompareAll)) ./ max(B_onAxis(1:nCompareAll));
plot(z1D(1:nCompareAll), 100*errorBfield, 'k', 'LineWidth', 3);
set(gca,'FontName','times','FontSize',11)
ylabel('|Error| (\%)','Interpreter','latex','FontSize',12)
xlabel('z [m]','Interpreter','latex','FontSize',13)
grid on
ylim([0,2])

%% SECTION 5: Optional save figure
% =========================================================================
saveFig = 0;   % set to 1 to save automatically

if saveFig
    saveas(gcf, 'Validating_MagneticField_Code', 'tiffn')
end

%% SECTION 6: Write NetCDF file
% =========================================================================
file = 'bfield_protoMPEX.nc';
if exist(file, 'file')
    delete(file);
end

r  = r1D(:);      % [nR x 1]
z  = z1D(:);      % [nZ x 1]
Br = Br2D;        % likely [nZ x nR]
Bz = Bz2D;        % likely [nZ x nR]
Bt = zeros(size(Bz2D), 'like', Bz2D);

% Coordinates
nccreate(file, 'r', 'Dimensions', {'r', numel(r)}, 'Datatype', 'double');
nccreate(file, 'z', 'Dimensions', {'z', numel(z)}, 'Datatype', 'double');
ncwrite(file, 'r', r);
ncwrite(file, 'z', z);

% Fields stored as [r x z]
nccreate(file, 'br', 'Dimensions', {'r', numel(r), 'z', numel(z)}, 'Datatype', 'single');
nccreate(file, 'bt', 'Dimensions', {'r', numel(r), 'z', numel(z)}, 'Datatype', 'single');
nccreate(file, 'bz', 'Dimensions', {'r', numel(r), 'z', numel(z)}, 'Datatype', 'single');

% Convert from [z x r] to [r x z] if needed
if isequal(size(Br), [numel(z), numel(r)])
    Br = permute(Br, [2 1]);
end
if isequal(size(Bt), [numel(z), numel(r)])
    Bt = permute(Bt, [2 1]);
end
if isequal(size(Bz), [numel(z), numel(r)])
    Bz = permute(Bz, [2 1]);
end

assert(isequal(size(Br), [numel(r), numel(z)]), 'Br size mismatch after permute');
assert(isequal(size(Bt), [numel(r), numel(z)]), 'Bt size mismatch after permute');
assert(isequal(size(Bz), [numel(r), numel(z)]), 'Bz size mismatch after permute');

ncwrite(file, 'br', single(Br));
ncwrite(file, 'bt', single(Bt));
ncwrite(file, 'bz', single(Bz));

%% SECTION 7: Read NetCDF back and plot Bz
% =========================================================================
r0  = ncread(file, 'r');
z0  = ncread(file, 'z');
br0 = ncread(file, 'br');
bt0 = ncread(file, 'bt');
bz0 = ncread(file, 'bz');

figure('Color','w')
imagesc(z0, r0, bz0)
set(gca, 'YDir', 'normal')
colorbar
xlabel('z [m]')
ylabel('r [m]')
title('B_z from NetCDF')

%% SECTION 8: Radial slice check
% =========================================================================
% Pick a representative axial index near midplane
iz = round(numel(z0)/2);

figure('Color','w')
plot(r0, bz0(:,iz), 'LineWidth', 2)
xlabel('r [m]')
ylabel('B_z [T]')
title(['Radial profile of B_z at z = ', num2str(z0(iz)), ' m'])
grid on
box on

%% SECTION 9: On-axis / near-axis field from NetCDF arrays
% =========================================================================
B0 = sqrt(br0.^2 + bz0.^2);

figure('Color','w')
plot(z0, B0(1,:), 'k', 'LineWidth', 2)
box on
grid on
set(gca, 'PlotBoxAspectRatio', [4 1 1])
set(gca, 'FontName', 'times')
xlabel('z [m]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('B$_0$ [T]', 'Interpreter', 'latex', 'FontSize', 13)
xlim([0, 4.5])
title('Near-axis magnetic field magnitude')

%% SECTION 10: Br/B0 profile at one radial index
% =========================================================================
ir = min( round(numel(r0)/2), size(br0,1) );

figure('Color','w')
plot(z0, br0(ir,:) ./ B0(ir,:), 'k', 'LineWidth', 2)
grid on
box on
xlabel('z [m]')
ylabel('B_r / B_0')
title(['B_r / B_0 at r = ', num2str(r0(ir)), ' m'])

%% SECTION 11: Export CSV
% =========================================================================
writematrix(z1D(:),  'z1D.csv')
writematrix(r1D(:),  'r1D.csv')
writematrix(Br2D.',  'Br2D.csv')
writematrix(Bz2D.',  'Bz2D.csv')
writematrix(B2D.',   'Bmag2D.csv')

disp('End of script')