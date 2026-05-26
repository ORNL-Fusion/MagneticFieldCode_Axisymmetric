% main_program.m - Compute magnetic flux, Br, Bz, and write ProtoLite B-field

clc;
clear all;
close all;

% ----- Read in coil data ----------
data1 = readtable("Coil_setup_I_locs_8inbackward.xlsx");

R1 = data1.r_inner;
R2 = data1.r_outer;
wZ = data1.dz;
X  = data1.z;
ps = data1.ps;

z_start = -0.25;
z_end   = 3.5;

% ----- Map current names ----------
mapping = containers.Map({'PS1', 'TR2', 'TR1'}, [1575, 220, 1575]);

I = zeros(size(ps));
for idx = 1:length(ps)
    I(idx) = mapping(ps{idx});
end

% ----- Define grid ----------
r1D = linspace(0.001, 0.13, 150);
z1D = linspace(z_start, z_end, 500);

% meshgrid gives arrays as [nZ x nR]
[R, Z] = meshgrid(r1D, z1D);

% ----- Compute magnetic vector potential A_phi ----------
A_phi_total = computeMagneticFlux(data1, R, Z, I);   % [nZ x nR]

% =========================================================================
% Compute psi, Br, and Bz from A_phi
% =========================================================================

psi_zr = 2*pi .* R .* A_phi_total;   % [nZ x nR]

dr = r1D(2) - r1D(1);
dz = z1D(2) - z1D(1);

% For matrix [z x r]:
% first output = derivative along z
% second output = derivative along r
[dpsidz, dpsidr] = gradient(psi_zr, dz, dr);

rSafe = R;
rSafe(rSafe < 1e-8) = 1e-8;

Br_zr = -dpsidz ./ (2*pi .* rSafe);
Bz_zr =  dpsidr ./ (2*pi .* rSafe);
Bt_zr = zeros(size(Bz_zr));

% Convert to [r x z] for NetCDF
Br = Br_zr.';
Bt = Bt_zr.';
Bz = Bz_zr.';
psi = psi_zr.';

r = r1D(:);
z = z1D(:);

fprintf('max(abs(Br)) before write = %.6e T\n', max(abs(Br(:))));
fprintf('max(abs(Bz)) before write = %.6e T\n', max(abs(Bz(:))));

% ----- Optional axial field comparison ----------
ZB = z1D;
B_total = computeAxialField(data1, ZB, 4*pi*1e-7, I);

figure;
subplot(2,1,1);
plot(ZB, B_total, 'k-', 'LineWidth', 2.5);
hold on;
plot(z1D, Bz(1,:), 'r--', 'LineWidth', 2.0);

for j = 1:length(X)
    xline(X(j), 'b--', ['Coil ' num2str(j)]);
end

grid on;
box on;
xlabel('z [m]');
ylabel('B_z [T]');
title('Axial magnetic field comparison');
legend('computeAxialField', 'B_z from \psi near axis', 'Location', 'best');

% ----- Plot device geometry and flux contours ----------
subplot(2,1,2);
Find_LCFS(A_phi_total, r1D, z1D, X, data1, z_start);

hold on;

for j = 1:length(X)
    dR = R2(j) - R1(j);

    rectangle('Position', [(X(j)-wZ(j)/2), R1(j), wZ(j), dR], ...
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
% Plot Br, Bz, psi
% =========================================================================

figure;
imagesc(z, r, Bz);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
colorbar;
xlabel('z [m]');
ylabel('r [m]');
title('ProtoLite B_z [T]');

figure;
imagesc(z, r, Br);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
colorbar;
xlabel('z [m]');
ylabel('r [m]');
title('ProtoLite B_r [T]');

figure;
imagesc(z, r, psi);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
colorbar;
xlabel('z [m]');
ylabel('r [m]');
title('ProtoLite \psi [Wb]');

% =========================================================================
% Save B-field profile to NetCDF
% =========================================================================

file = 'bfield_protoLite.nc';

if exist(file, 'file')
    delete(file);
end

assert(isequal(size(Br), [numel(r), numel(z)]), 'Br size mismatch');
assert(isequal(size(Bt), [numel(r), numel(z)]), 'Bt size mismatch');
assert(isequal(size(Bz), [numel(r), numel(z)]), 'Bz size mismatch');
assert(isequal(size(psi), [numel(r), numel(z)]), 'psi size mismatch');

nccreate(file, 'r', ...
    'Dimensions', {'r', numel(r)}, ...
    'Datatype', 'double');

nccreate(file, 'z', ...
    'Dimensions', {'z', numel(z)}, ...
    'Datatype', 'double');

ncwrite(file, 'r', r);
ncwrite(file, 'z', z);

nccreate(file, 'br', ...
    'Dimensions', {'r', numel(r), 'z', numel(z)}, ...
    'Datatype', 'single');

nccreate(file, 'bt', ...
    'Dimensions', {'r', numel(r), 'z', numel(z)}, ...
    'Datatype', 'single');

nccreate(file, 'bz', ...
    'Dimensions', {'r', numel(r), 'z', numel(z)}, ...
    'Datatype', 'single');

nccreate(file, 'psi', ...
    'Dimensions', {'r', numel(r), 'z', numel(z)}, ...
    'Datatype', 'double');

ncwrite(file, 'br', single(Br));
ncwrite(file, 'bt', single(Bt));
ncwrite(file, 'bz', single(Bz));
ncwrite(file, 'psi', psi);

% =========================================================================
% Read-back check
% =========================================================================

br0 = ncread(file, 'br');
bz0 = ncread(file, 'bz');

fprintf('max(abs(Br)) after write = %.6e T\n', max(abs(br0(:))));
fprintf('max(abs(Bz)) after write = %.6e T\n', max(abs(bz0(:))));

figure;
imagesc(z, r, br0);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('z [m]');
ylabel('r [m]');
title('Read-back B_r from NetCDF');

figure;
imagesc(z, r, bz0);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('z [m]');
ylabel('r [m]');
title('Read-back B_z from NetCDF');

% =========================================================================
% Save CSV files
% =========================================================================

writematrix(z,   'z1D.csv');
writematrix(r,   'r1D.csv');
writematrix(Br,  'Br2D.csv');
writematrix(Bz,  'Bz2D.csv');
writematrix(psi, 'psi2D.csv');

disp('Saved B-field profile to bfield_protoLite.nc');
disp('End of script');