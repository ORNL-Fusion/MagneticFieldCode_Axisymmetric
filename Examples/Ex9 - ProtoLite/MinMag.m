% Minimize Mag %

clc;
clear all;
close all;

tag = 'cur_design';

base_currents = [1575, 220, 1575];
range_currents = [[1 2000];[1 280];[1 2000]];
n_samples = 5000;

% ----- Read in coil data ----------
data1 = readtable("Coil_setup_I_locs_8inbackward.xlsx");

R1 = data1.r_inner;
R2 = data1.r_outer;
wZ = data1.dz;
X  = data1.z;
ps = data1.ps;

ds.z_start = -0.25;
ds.z_end = 3.5;
ds.z_n = 500;
ds.r_start = 0.001;
ds.r_end = 0.13;
ds.r_n = 150;
ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, base_currents);



parms=lhsu(range_currents(:,1),range_currents(:,2),n_samples);

% ----- Map current names ----------


I = zeros(size(ps));
for idx = 1:length(ps)
    I(idx) = ds.mapping(ps{idx});
end

% ----- Define grid ----------
r1D = linspace(ds.r_start, ds.r_end, ds.r_n);
z1D = linspace(ds.z_start, ds.z_end, ds.z_n);

r = r1D(:);
z = z1D(:);

ZB = z1D;

dr = r1D(2) - r1D(1);
dz = z1D(2) - z1D(1);

% meshgrid gives arrays as [nZ x nR]
[R, Z] = meshgrid(r1D, z1D);

% ----- Compute magnetic vector potential A_phi ----------
A_phi_total = computeMagneticFlux(data1, R, Z, I);   % [nZ x nR]

% =========================================================================
% Compute psi, Br, and Bz from A_phi
% =========================================================================

psi_zr = 2*pi .* R .* A_phi_total;   % [nZ x nR]

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

fprintf('max(abs(Br)) before write = %.6e T\n', max(abs(Br(:))));
fprintf('max(abs(Bz)) before write = %.6e T\n', max(abs(Bz(:))));

% ----- Optional axial field comparison ----------

B_total = computeAxialField(data1, ZB, 4*pi*1e-7, I);

figure;
subplot(2,1,1);
plot(ZB, B_total, 'k-', 'LineWidth', 2.5);
hold on;
plot(z1D, Bz(1,:), 'r--', 'LineWidth', 2.0);

for j = 1:length(X)
    xline(X(j), 'b--', ['Coil ' num2str(j) ' (' ps{j} ')']);
end

grid on;
box on;
xlabel('z [m]');
ylabel('B_z [T]');
title(strcat('Axial magnetic field comparison(','PS1=',num2str(base_currents(1)),' TR2=',num2str(base_currents(2)), ' TR1=',num2str(base_currents(3)),')'));
legend('computeAxialField', 'B_z from \psi near axis', 'Location', 'best');

subplot(2,1,2);
[C_min, max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
ds.C_min = C_min;
Find_LCFS_ML(A_phi_total, r1D, z1D, X, data1, ds);

us_z = find(abs(ZB-1)<.01);
us_zc = find(abs(ZB-1.795)<.01);
zc_val = B_total(us_zc);
max_org = max(B_total(us_z));


opt_var = zeros(n_samples,6);
B_var = zeros(n_samples,size(B_total,2));
parms(1,:) = base_currents;
h = waitbar(0, 'Sampling');
for j=1:n_samples
    waitbar(j/n_samples, h)
    % ----- Map current names ----------
    ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, parms(j,:));

    I = zeros(size(ps));
    for idx = 1:length(ps)
        I(idx) = ds.mapping(ps{idx});
    end

    % ----- Compute magnetic vector potential A_phi ----------
    A_phi_total = computeMagneticFlux(data1, R, Z, I);   % [nZ x nR]

    [C_min,C_max,max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
    opt_var(j,:) = [C_min;C_max;max_p_hel(:);max_p_lim(:)]';

    psi_zr = 2*pi .* R .* A_phi_total;   % [nZ x nR]
    [dpsidz, dpsidr] = gradient(psi_zr, dz, dr);
    rSafe = R;
    rSafe(rSafe < 1e-8) = 1e-8;
    Br_zr = -dpsidz ./ (2*pi .* rSafe);
    Bz_zr =  dpsidr ./ (2*pi .* rSafe);
    Bt_zr = zeros(size(Bz_zr));
    Br = Br_zr.';
    Bt = Bt_zr.';
    Bz = Bz_zr.';
    psi = psi_zr.';
    B_total = computeAxialField(data1, ZB, 4*pi*1e-7, I);
    B_var(j,:) = B_total;
end

save(strcat('data_',tag,'.mat'),'parms','opt_var','B_var');


close(h)
rind = find(opt_var(:,6)>.06 & max(abs(B_var(:,us_zc)-repmat(zc_val,[n_samples 1])),[],2)<.01);


[min_v,bind] = max(min(B_var(rind,us_z),[],2));
us_ind = rind(bind);

ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, parms(us_ind,:));

I = zeros(size(ps));
for idx = 1:length(ps)
    I(idx) = ds.mapping(ps{idx});
end

% ----- Compute magnetic vector potential A_phi ----------
A_phi_total = computeMagneticFlux(data1, R, Z, I);   % [nZ x nR]

[C_min,C_max,max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
psi_zr = 2*pi .* R .* A_phi_total;   % [nZ x nR]
[dpsidz, dpsidr] = gradient(psi_zr, dz, dr);
rSafe = R;
rSafe(rSafe < 1e-8) = 1e-8;
Br_zr = -dpsidz ./ (2*pi .* rSafe);
Bz_zr =  dpsidr ./ (2*pi .* rSafe);
Bt_zr = zeros(size(Bz_zr));
Br = Br_zr.';
Bt = Bt_zr.';
Bz = Bz_zr.';
psi = psi_zr.';
B_total = computeAxialField(data1, ZB, 4*pi*1e-7, I);

figure;
subplot(2,1,1);
plot(ZB, B_total, 'k-', 'LineWidth', 2.5);
hold on;
plot(z1D, Bz(1,:), 'r--', 'LineWidth', 2.0);

for j = 1:length(X)
    xline(X(j), 'b--', ['Coil ' num2str(j) ' (' ps{j} ')']);
end

grid on;
box on;
xlabel('z [m]');
ylabel('B_z [T]');
title(strcat('Axial magnetic field comparison(','PS1=',num2str(parms(us_ind,1)),' TR2=',num2str(parms(us_ind,2)), ' TR1=',num2str(parms(us_ind,3)),')'));
legend('computeAxialField', 'B_z from \psi near axis', 'Location', 'best');

subplot(2,1,2);
[C_min, max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
ds.C_min = C_min;
Find_LCFS_ML(A_phi_total, r1D, z1D, X, data1, ds);


[min_v,bind] = min(min(B_var(rind,us_z),[],2));
us_ind = rind(bind);


ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, parms(us_ind,:));

I = zeros(size(ps));
for idx = 1:length(ps)
    I(idx) = ds.mapping(ps{idx});
end

% ----- Compute magnetic vector potential A_phi ----------
A_phi_total = computeMagneticFlux(data1, R, Z, I);   % [nZ x nR]

[C_min,C_max,max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
psi_zr = 2*pi .* R .* A_phi_total;   % [nZ x nR]
[dpsidz, dpsidr] = gradient(psi_zr, dz, dr);
rSafe = R;
rSafe(rSafe < 1e-8) = 1e-8;
Br_zr = -dpsidz ./ (2*pi .* rSafe);
Bz_zr =  dpsidr ./ (2*pi .* rSafe);
Bt_zr = zeros(size(Bz_zr));
Br = Br_zr.';
Bt = Bt_zr.';
Bz = Bz_zr.';
psi = psi_zr.';
B_total = computeAxialField(data1, ZB, 4*pi*1e-7, I);

figure;
subplot(2,1,1);
plot(ZB, B_total, 'k-', 'LineWidth', 2.5);
hold on;
plot(z1D, Bz(1,:), 'r--', 'LineWidth', 2.0);

for j = 1:length(X)
    xline(X(j), 'b--', ['Coil ' num2str(j) ' (' ps{j} ')']);
end

grid on;
box on;
xlabel('z [m]');
ylabel('B_z [T]');
title(strcat('Axial magnetic field comparison(','PS1=',num2str(parms(us_ind,1)),' TR2=',num2str(parms(us_ind,2)), ' TR1=',num2str(parms(us_ind,3)),')'));
legend('computeAxialField', 'B_z from \psi near axis', 'Location', 'best');

subplot(2,1,2);
[C_min, max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
ds.C_min = C_min;
Find_LCFS_ML(A_phi_total, r1D, z1D, X, data1, ds);



rmax = .6;
[min_v,bind] = sort(min(B_var(rind,us_z),[],2), 'descend');
n_frames = length(find(min_v>max_org));
h_fig = figure('Visible','on');
h_fig.WindowState = 'maximized';
set(gca, 'NextPlot', 'replaceChildren');
v = VideoWriter(strcat('movie_max_',tag,'.mat'), 'MPEG-4');
v.FrameRate = 2;
open(v);

for jf = 1:n_frames
    us_ind = rind(bind(jf));

    ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, parms(us_ind,:));

    I = zeros(size(ps));
    for idx = 1:length(ps)
        I(idx) = ds.mapping(ps{idx});
    end

    % ----- Compute magnetic vector potential A_phi ----------
    A_phi_total = computeMagneticFlux(data1, R, Z, I);   % [nZ x nR]

    [C_min,C_max,max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
    psi_zr = 2*pi .* R .* A_phi_total;   % [nZ x nR]
    [dpsidz, dpsidr] = gradient(psi_zr, dz, dr);
    rSafe = R;
    rSafe(rSafe < 1e-8) = 1e-8;
    Br_zr = -dpsidz ./ (2*pi .* rSafe);
    Bz_zr =  dpsidr ./ (2*pi .* rSafe);
    Bt_zr = zeros(size(Bz_zr));
    Br = Br_zr.';
    Bt = Bt_zr.';
    Bz = Bz_zr.';
    psi = psi_zr.';
    B_total = computeAxialField(data1, ZB, 4*pi*1e-7, I);

    figure(h_fig);
    clf
    subplot(2,1,1);
    plot(ZB, B_total, 'k-', 'LineWidth', 2.5);
    hold on;
    plot(z1D, Bz(1,:), 'r--', 'LineWidth', 2.0);

    for j = 1:length(X)
        xline(X(j), 'b--', ['Coil ' num2str(j) ' (' ps{j} ')']);
    end

    grid on;
    box on;
    xlabel('z [m]');
    ylabel('B_z [T]');
    title(strcat('Axial magnetic field comparison(','PS1=',num2str(parms(us_ind,1)),' TR2=',num2str(parms(us_ind,2)), ' TR1=',num2str(parms(us_ind,3)),', max B_T =',num2str(min(B_var(us_ind,us_z))),')'));
    legend('computeAxialField', 'B_z from \psi near axis', 'Location', 'best');
     xlim([ds.z_start ds.z_end]);
    ylim([0 rmax]);

    subplot(2,1,2);
    [C_min, max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
    ds.C_min = C_min;
    Find_LCFS_ML(A_phi_total, r1D, z1D, X, data1, ds);
    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v);



[min_v,bind] = sort(min(B_var(rind,us_z),[],2));
h_fig = figure('Visible','on');
h_fig.WindowState = 'maximized';
set(gca, 'NextPlot', 'replaceChildren');
v = VideoWriter(strcat('movie_min_',tag,'.mat'), 'MPEG-4');
v.FrameRate = 2;
open(v);

for jf = 1:n_frames
    us_ind = rind(bind(jf));

    ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, parms(us_ind,:));

    I = zeros(size(ps));
    for idx = 1:length(ps)
        I(idx) = ds.mapping(ps{idx});
    end

    % ----- Compute magnetic vector potential A_phi ----------
    A_phi_total = computeMagneticFlux(data1, R, Z, I);   % [nZ x nR]

    [C_min,C_max,max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
    psi_zr = 2*pi .* R .* A_phi_total;   % [nZ x nR]
    [dpsidz, dpsidr] = gradient(psi_zr, dz, dr);
    rSafe = R;
    rSafe(rSafe < 1e-8) = 1e-8;
    Br_zr = -dpsidz ./ (2*pi .* rSafe);
    Bz_zr =  dpsidr ./ (2*pi .* rSafe);
    Bt_zr = zeros(size(Bz_zr));
    Br = Br_zr.';
    Bt = Bt_zr.';
    Bz = Bz_zr.';
    psi = psi_zr.';
    B_total = computeAxialField(data1, ZB, 4*pi*1e-7, I);

    figure(h_fig);
    clf
    subplot(2,1,1);
    plot(ZB, B_total, 'k-', 'LineWidth', 2.5);
    hold on;
    plot(z1D, Bz(1,:), 'r--', 'LineWidth', 2.0);

    for j = 1:length(X)
        xline(X(j), 'b--', ['Coil ' num2str(j) ' (' ps{j} ')']);
    end

    grid on;
    box on;
    xlabel('z [m]');
    ylabel('B_z [T]');
    title(strcat('Axial magnetic field comparison(','PS1=',num2str(parms(us_ind,1)),' TR2=',num2str(parms(us_ind,2)), ' TR1=',num2str(parms(us_ind,3)),', max B_T =',num2str(min(B_var(us_ind,us_z))),')'));
    legend('computeAxialField', 'B_z from \psi near axis', 'Location', 'best');
     xlim([ds.z_start ds.z_end]);
    ylim([0 rmax]);

    subplot(2,1,2);
    [C_min, max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
    ds.C_min = C_min;
    Find_LCFS_ML(A_phi_total, r1D, z1D, X, data1, ds);
    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v);



