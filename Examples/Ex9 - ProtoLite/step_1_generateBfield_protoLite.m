% main_program.m - Main script to compute and plot magnetic flux and axial field

clc;
clear all;

% ----- Read in coil data and assign coil parameters ----------
data1 = readtable("Coil_setup_I_locs_8inbackward.xlsx");

% Coil parameters
nR = data1.layers_r;
nZ = data1.layers_z;
R1 = data1.r_inner;
R2 = data1.r_outer;
wZ = data1.dz;
X  = data1.z;
ps = data1.ps;

z_start = -0.25;   % beginning of dump chamber
z_end   = 3.5;     % end of device

% ----- Map string current names to numeric current values ----------
mapping = containers.Map({'PS1', 'TR2', 'TR1'}, [1575, 220, 1575]);

I = zeros(size(ps));
for idx = 1:length(ps)
    I(idx) = mapping(ps{idx});
end

% ----- Define grid for flux calculation ----------
y = linspace(0.001, 0.13, 150);       % radial grid [m]
Z = linspace(z_start, z_end, 500);    % axial grid [m]

[Y, Z_matrix] = meshgrid(y, Z);
dist = sqrt((Z_matrix).^2 + (Y).^2);

% ----- Compute magnetic flux function ----------
A_phi_total = computeMagneticFlux(data1, Y, Z_matrix, I);

% ----- Compute axial magnetic field profile ----------
ZB = linspace(z_start, z_end, 500);
B_total = computeAxialField(data1, ZB, 4*pi*1e-7, I);

b_data = [ZB', B_total'];

% ----- Plot axial field profile ----------
figure;
subplot(2,1,1);

plot(ZB, B_total, 'LineWidth', 2.5);
hold on;

for j = 1:length(X)
    xline(X(j), 'r--', ['Coil ' num2str(j)]);
end

grid on;
hold off;

xlabel('Z (m)');
ylabel('Axial magnetic field');
title('Axial Magnetic Field Profile');

% ----- Plot device geometry and flux contours ----------
subplot(2,1,2);

Find_LCFS(A_phi_total, y, Z, X, data1, z_start);

hold on;

for j = 1:length(X)
    dR = R2(j) - R1(j);

    rectangle('Position', [(X(j)-wZ(j)/2), R1(j),  wZ(j), dR], ...
              'FaceColor', [1 0 0], ...
              'EdgeColor', 'r', ...
              'LineWidth', 1);

    rectangle('Position', [(X(j)-wZ(j)/2), -R2(j), wZ(j), dR], ...
              'FaceColor', [1 0 0], ...
              'EdgeColor', 'r', ...
              'LineWidth', 1);
end

hold off;

% =========================================================================
% Save B-field profile to NetCDF and plot/check output
% =========================================================================

file = 'bfield_protoLite.nc';

if exist(file, 'file')
    delete(file);
end

r  = y(:);       % radial coordinate [nR x 1]
z1D = Z(:);      % axial coordinate [nZ x 1]

% This script computes only axial B-field.
% Therefore Br and Bt are set to zero unless you compute them elsewhere.
Br = zeros(numel(r), numel(z1D), 'single');
Bt = zeros(numel(r), numel(z1D), 'single');

% Interpolate axial B profile onto the Z grid
Bz_line = interp1(ZB(:), B_total(:), z1D, 'linear', 'extrap');

% Replicate axial B-field across radius
Bz = repmat(Bz_line(:).', numel(r), 1);   % [nR x nZ]

% --- declare coordinates
nccreate(file, 'r', ...
    'Dimensions', {'r', numel(r)}, ...
    'Datatype', 'double');

nccreate(file, 'z', ...
    'Dimensions', {'z', numel(z1D)}, ...
    'Datatype', 'double');

ncwrite(file, 'r', r);
ncwrite(file, 'z', z1D);

% --- declare fields in [r x z] order
nccreate(file, 'br', ...
    'Dimensions', {'r', numel(r), 'z', numel(z1D)}, ...
    'Datatype', 'single');

nccreate(file, 'bt', ...
    'Dimensions', {'r', numel(r), 'z', numel(z1D)}, ...
    'Datatype', 'single');

nccreate(file, 'bz', ...
    'Dimensions', {'r', numel(r), 'z', numel(z1D)}, ...
    'Datatype', 'single');

% --- sanity checks
assert(isequal(size(Br), [numel(r), numel(z1D)]), ...
    'Br size mismatch');

assert(isequal(size(Bt), [numel(r), numel(z1D)]), ...
    'Bt size mismatch');

assert(isequal(size(Bz), [numel(r), numel(z1D)]), ...
    'Bz size mismatch');

% --- write fields
ncwrite(file, 'br', single(Br));
ncwrite(file, 'bt', single(Bt));
ncwrite(file, 'bz', single(Bz));

% --- read/plot check
r0  = ncread(file, 'r');
z0  = ncread(file, 'z');
bz0 = ncread(file, 'bz');   % [nR x nZ]

figure;
imagesc(z0, r0, bz0);
set(gca, 'YDir', 'normal');
colorbar;

xlabel('z [m]');
ylabel('r [m]');
title('B_z from NetCDF');

% --- save CSV files
writematrix(z1D, 'z1D.csv');
writematrix(r,   'r1D.csv');
writematrix(Br,  'Br2D.csv');
writematrix(Bz,  'Bz2D.csv');

disp('Saved B-field profile to bfield_protoMPEX.nc');
disp('End of script');