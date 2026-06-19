% Minimize Mag %

clc;
clear all;
close all;

tag = 'table_design';

PS1_base_current = 1600;
m_targ = .05;
range_currents = [[0 280];[0 2000]];

TS1_dis = [1000 1300 1600 1900];
PS1_dis = [1000 1300 1600 1900];
lu_tab_v = linspace(1000,2000,11);

headers = {'PS1', 'TR1', 'TR2'};
n_lab = length(lu_tab_v);
z_lab = zeros(n_lab,1);
lu_tab = table(z_lab+PS1_base_current,lu_tab_v(:),z_lab,'VariableNames', headers);

n_samples = 3;
n_interp = 20;
z_point = 1;
% do_test_mid = 0;


% ----- Read in coil data ----------
data1 = readtable("Coil_setup_I_locs_8inbackward.xlsx");

R1 = data1.r_inner;
R2 = data1.r_outer;
wZ = data1.dz;
X  = data1.z;
ps = data1.ps;

ds.z_start = -0.25;
ds.z_end = 3.5;
ds.z_n = 376;
ds.r_start = 0.001;
ds.r_end = 0.13;
ds.r_n = 100;

TR1_cur = linspace(range_currents(2,1), range_currents(2,2),n_samples);
TR2_cur = linspace(range_currents(1,1), range_currents(1,2),n_samples);

TR1_curS = linspace(range_currents(2,1), range_currents(2,2),n_samples*n_interp);
TR2_curS = linspace(range_currents(1,1), range_currents(1,2),n_samples*n_interp);

[X,Y] = meshgrid(TR1_cur,TR2_cur);
[Xq,Yq] = meshgrid(TR1_curS,TR2_curS);


% ----- Define grid ----------
r1D = linspace(ds.r_start, ds.r_end, ds.r_n);
z1D = linspace(ds.z_start, ds.z_end, ds.z_n);
[minv z1_ind] = min(abs(z1D-z_point));

r = r1D(:);
z = z1D(:);

ZB = z1D;

dr = r1D(2) - r1D(1);
dz = z1D(2) - z1D(1);

% meshgrid gives arrays as [nZ x nR]
[R, Z] = meshgrid(r1D, z1D);


BT_z1 = zeros(n_samples,n_samples);
h = waitbar(0, 'Sampling');
for j1=1:n_samples
    for j2=1:n_samples
        waitbar((j2+(j1-1)*n_samples)/(n_samples*n_samples), h)
        cparm = [PS1_base_current TR2_cur(j2) TR1_cur(j1)];
        % ----- Map current names ----------
        ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, cparm);

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
        BT_z1(j1,j2) = B_total(z1_ind);
    end
end
close(h)
figure,
imagesc(TR1_cur,TR2_cur,BT_z1)
xlabel('TR1 Current (A)')
ylabel('TR2 Current (A)')
title(strcat('Total Magenetic Field at z=',num2str(z_point),'m'))
c = colorbar;             % Create the colorbar
title(c,'Tesla');  % Add the label
c.Label.FontSize = 12; 

Surr_BT_z1 = interp2(X,Y,BT_z1,Xq,Yq);

figure,
imagesc(TR1_curS,TR2_curS,Surr_BT_z1)
xlabel('TR1 Current (A)')
ylabel('TR2 Current (A)')
title(strcat('Surrogate of Total Magenetic Field at z=',num2str(z_point),'m'))
c = colorbar;             % Create the colorbar
title(c,'Tesla');  % Add the label
c.Label.FontSize = 12; 


BT_z1_true = zeros(n_samples,n_samples);
BT_z1_est = zeros(n_samples,n_samples);

C_matrix = contourc(TR2_curS,TR1_curS,Surr_BT_z1,[m_targ m_targ]);
x_data = C_matrix(1, 2:end); 
y_data = C_matrix(2, 2:end);

p = polyfit(y_data,x_data, 1);

BT_z1p1 = zeros(n_samples,n_samples);
h = waitbar(0, 'Sampling for PS1');
for j1=1:n_samples
    for j2=1:n_samples
        waitbar((j2+(j1-1)*n_samples)/(n_samples*n_samples), h)
        cparm = [PS1_base_current-1 TR2_cur(j2) TR1_cur(j1)];
        % ----- Map current names ----------
        ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, cparm);

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
        BT_z1p1(j1,j2) = B_total(z1_ind);
    end
end
close(h)

Surr_BT_z1p1 = interp2(X,Y,BT_z1p1,Xq,Yq);
C_matrixp1 = contourc(TR2_curS,TR1_curS,Surr_BT_z1p1,[m_targ m_targ]);
x_datap1 = C_matrixp1(1, 2:end); 
y_datap1 = C_matrixp1(2, 2:end);

pp1 = polyfit(y_datap1,x_datap1, 1);

m = (p(2)-pp1(2));
b = p(2)-m*(PS1_base_current);

h = waitbar(0, 'Generating Figure');
figure,
sgtitle(strcat('TS2 = ',num2str(p(1)),'\cdot','TS1',num2str(m),'\cdot','PS1+',num2str(b)));
subplot(length(TS1_dis),length(PS1_dis),1)
ax = zeros(length(PS1_dis)*length(TS1_dis),1);
maxy = 0;
for jf = 1:length(TS1_dis)
    for jp = 1:length(PS1_dis)
        waitbar((jp+(jf-1)*length(TS1_dis))/(length(TS1_dis)*length(PS1_dis)), h)
        cparm = [PS1_dis(jp) p(1)*TS1_dis(jf)+m*PS1_dis(jp)+b TS1_dis(jf)];

        ds.mapping = containers.Map({'PS1', 'TR2', 'TR1'}, cparm);

        I = zeros(size(ps));
        for idx = 1:length(ps)
            I(idx) = ds.mapping(ps{idx});
        end

        % ----- Compute magnetic vector potential A_phi ----------
        A_phi_total = computeMagneticFlux(data1, R, Z, I);   % [nZ x nR]

        [C_min,C_max,max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, r1D, z1D, X, data1, ds.z_start);
        opt_var(j1,j2,:) = [C_min;C_max;max_p_hel(:);max_p_lim(:)]';

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
        ax((jf-1)*length(PS1_dis)+jp) = subplot(length(TS1_dis),length(PS1_dis),(jf-1)*length(PS1_dis)+jp);
        plot(ZB, B_total, 'k-', 'LineWidth', 2.5);
        title(strcat('B total (','PS1=',num2str(cparm(1)),' TR2=',num2str(cparm(2)), ' TR1=',num2str(cparm(3)),', max B_T =',num2str(B_total(z1_ind)),')'))
        maxy = max(maxy,max(B_total));
    end
end
close(h)

linkaxes(ax, 'xy'); 
xlim([ds.z_start, ds.z_end]);
ylim([0, 1.1*maxy]);

for j = 1:n_lab
    lu_tab.TR2(j) = p(1)*lu_tab.TR1(j)+m*lu_tab.PS1(j)+b;
end
disp(lu_tab)


