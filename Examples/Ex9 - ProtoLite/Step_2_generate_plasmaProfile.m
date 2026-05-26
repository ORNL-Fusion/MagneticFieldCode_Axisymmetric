%% ProtoLite: Compute ne, te using psi formulation and write to one NetCDF

close all; clear; clc;

% -------------------------------------------------------------------------
% READ EXISTING PROTOLITE B-FIELD FILE
% -------------------------------------------------------------------------
fileB = 'bfield_protoLite.nc';

r  = ncread(fileB,'r');      % [m]
z  = ncread(fileB,'z');      % [m]
Br = ncread(fileB,'br');     % [r x z]
Bt = ncread(fileB,'bt');     % [r x z]
Bz = ncread(fileB,'bz');     % [r x z]

if ~isequal(size(Bz), [numel(r) numel(z)])
    Br = permute(Br,[2 1]);
    Bt = permute(Bt,[2 1]);
    Bz = permute(Bz,[2 1]);
end

nR = numel(r);
nZ = numel(z);

zL = z(1);
zR = z(end);

% -------------------------------------------------------------------------
% PLOT INPUT Bz FIELD
% -------------------------------------------------------------------------
figure;
imagesc(z, r, Bz);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite input $B_z$ [T]','Interpreter','latex');
colorbar;
axis image tight;
xlim([zL zR]);

% -------------------------------------------------------------------------
% CALCULATE PSI
% -------------------------------------------------------------------------
disp('Calculating ProtoLite psi...');

dr = r(2) - r(1);
dz = z(2) - z(1);

psi = zeros(nR, nZ);

for ii = 1:nR
    for jj = 1:nZ
        psi(ii,jj) = 2*pi * sum(Bz(1:ii,jj) .* r(1:ii)) * dr;
    end
end

% Reference location for psi normalization
z0 = min(max(1.0, zL), zR);
r0 = min(max(0.06, r(1)), r(end));

psi_norm_val = interp2(z, r, psi, z0, r0, 'linear');

if ~isfinite(psi_norm_val) || abs(psi_norm_val) < 1e-30
    error('Bad psi normalization value. Check z0/r0 and B-field file.');
end

psiN = psi ./ psi_norm_val;

figure;
imagesc(z, r, psiN);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite normalized $\psi_N$','Interpreter','latex');
colorbar;
xlim([zL zR]);

% -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES
% -------------------------------------------------------------------------
disp('Constructing ProtoLite Te and ne profiles...');

qe = 1.602176634e-19;

Te0       = 2.0;     % eV
Te_edge   = 6.0;     % eV
Te_floor  = 0.2;     % eV
alpha     = 1.75;
lambdaPsi = 0.15;

Te = zeros(size(psiN));

inside = psiN <= 1;
outside = psiN > 1;

Te(inside) = Te_edge + (Te0 - Te_edge) .* (1 - psiN(inside)).^alpha;
Te(outside) = Te_floor + (Te_edge - Te_floor) .* ...
    exp(-(psiN(outside) - 1) ./ lambdaPsi);

% Density model
gasPuff_Gps = 1.0e20;
PRF_kW      = 180.0;
Eion_eV     = 1.0e3;

ne_min = 1.0e16;
ne_cap = 1.0e20;
ne_max = min(max(ne_cap, 5.0e18), 5.0e20);

Preq_kW = gasPuff_Gps * Eion_eV * qe / 1e3;
powerFactor = min(PRF_kW / max(Preq_kW,1e-12), 1.0);

ne_ref_for_power = 1.0e20;
Preq_ne_kW = Preq_kW * (ne_max / ne_ref_for_power);
powerMargin = PRF_kW / max(Preq_ne_kW,1e-12);

Ne = ne_min .* ones(size(psiN));
Ne(inside) = (ne_max - ne_min) .* (1 - psiN(inside)).^1.75 + ne_min;

disp(sprintf(['ProtoLite ne_max = %.3e m^-3\n', ...
    'Required power for ne_max: %.2f kW, available PRF = %.1f kW'], ...
    ne_max, Preq_ne_kW, PRF_kW));

% -------------------------------------------------------------------------
% PARALLEL VELOCITY PROFILE
% ProtoLite: helicon window centered at z = 1.0 m
% Stagnation point is imposed at z = 1.0 m
% -------------------------------------------------------------------------
disp('Constructing ProtoLite parallel velocity profile with stagnation at z = 1.0 m...');

mp = 1.67262192369e-27;
Cs = sqrt(qe .* Te ./ mp);

zStag = 1.0;       % helicon-window center / stagnation point [m]

% Clamp to simulation domain
zStag = min(max(zStag, zL), zR);

% Endpoints for Mach profile
zLeft  = zL;
zRight = zR;

Mz = zeros(1,nZ);

% Left side: M = -1 at left boundary, rising to 0 at zStag
idxL = z < zStag;
if zStag > zLeft
    Mz(idxL) = -1 .* (zStag - z(idxL)) ./ (zStag - zLeft);
else
    Mz(idxL) = 0;
end

% Right side: M = 0 at zStag, rising to +1 at right boundary
idxR = z > zStag;
if zRight > zStag
    Mz(idxR) = +1 .* (z(idxR) - zStag) ./ (zRight - zStag);
else
    Mz(idxR) = 0;
end

% Exactly at/nearest stagnation point: M = 0
[~, iStag] = min(abs(z - zStag));
Mz(iStag) = 0;

% Clip for numerical safety
Mz = max(min(Mz,1),-1);

M = repmat(Mz,nR,1);
Vpar = M .* Cs;
Mach = Vpar ./ Cs;
% -------------------------------------------------------------------------
% PARALLEL FLUXES
% -------------------------------------------------------------------------
disp('Constructing ProtoLite parallel particle and heat fluxes...');

GammaPar = Ne .* Vpar;
Qpar = (5/2) * qe .* Te .* GammaPar;

Bmag = sqrt(Br.^2 + Bt.^2 + Bz.^2);
bvec = Bmag(:);
bvec = bvec(isfinite(bvec) & bvec > 0);

if isempty(bvec)
    b_floor = 1e-12;
else
    b_floor = max(prctile(bvec,2),1e-12);
end

Bmag_safe = max(Bmag,b_floor);

bhat_r = Br ./ Bmag_safe;
bhat_t = Bt ./ Bmag_safe;
bhat_z = Bz ./ Bmag_safe;

GammaPar_r = GammaPar .* bhat_r;
GammaPar_t = GammaPar .* bhat_t;
GammaPar_z = GammaPar .* bhat_z;

Qpar_r = Qpar .* bhat_r;
Qpar_t = Qpar .* bhat_t;
Qpar_z = Qpar .* bhat_z;

% -------------------------------------------------------------------------
% PERPENDICULAR PARTICLE FLUX
% -------------------------------------------------------------------------
disp('Computing ProtoLite perpendicular particle flux...');

Dperp = 0.45;

[dndr, dndz] = gradient(Ne, dr, dz);

GammaPerp_r = -Dperp .* dndr;
GammaPerp_z = -Dperp .* dndz;

rMat = repmat(r(:),1,nZ);

dpsidr =  2*pi .* rMat .* Bz;
dpsidz = -2*pi .* rMat .* Br;

gradPsiMag = sqrt(dpsidr.^2 + dpsidz.^2);

gvec = gradPsiMag(:);
gvec = gvec(isfinite(gvec) & gvec > 0);

if isempty(gvec)
    eps_floor = 1e-12;
else
    eps_floor = max(prctile(gvec,2),1e-12);
end

gradPsiMag_safe = max(gradPsiMag,eps_floor);

nr = dpsidr ./ gradPsiMag_safe;
nz = dpsidz ./ gradPsiMag_safe;

mask_bad = gradPsiMag < eps_floor;
nr(mask_bad) = 0;
nz(mask_bad) = 0;

GammaPerp = GammaPerp_r .* nr + GammaPerp_z .* nz;
GammaPerp(~isfinite(GammaPerp)) = 0;

% -------------------------------------------------------------------------
% PERPENDICULAR HEAT FLUX
% -------------------------------------------------------------------------
disp('Computing ProtoLite perpendicular heat flux...');

chiPerp = 1.0;

[dTdr, dTdz] = gradient(Te, dr, dz);

Qperp_r = -chiPerp .* Ne .* qe .* dTdr;
Qperp_z = -chiPerp .* Ne .* qe .* dTdz;

Qperp = Qperp_r .* nr + Qperp_z .* nz;
Qperp(~isfinite(Qperp)) = 0;

% -------------------------------------------------------------------------
% PSI-NORMAL PARALLEL LEAKAGE CHECK
% -------------------------------------------------------------------------
bDotNpsi = bhat_r .* nr + bhat_z .* nz;

bDotN_tol = 5e-3;
bDotNpsi_f = bDotNpsi;
bDotNpsi_f(abs(bDotNpsi_f) < bDotN_tol) = 0;

GammaPar_psi = GammaPar .* bDotNpsi_f;
GammaPar_psi(~isfinite(GammaPar_psi)) = 0;

Qpar_psi = Qpar .* bDotNpsi_f;
Qpar_psi(~isfinite(Qpar_psi)) = 0;

% -------------------------------------------------------------------------
% TOTAL FLUXES
% -------------------------------------------------------------------------
GammaTot_r = GammaPar_r + GammaPerp_r;
GammaTot_t = GammaPar_t;
GammaTot_z = GammaPar_z + GammaPerp_z;
GammaTot_psi = GammaPerp;

Qtot_r = Qpar_r + Qperp_r;
Qtot_t = Qpar_t;
Qtot_z = Qpar_z + Qperp_z;
Qtot_psi = Qperp;

% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------
figure; imagesc(z,r,Te);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $T_e$ [eV]','Interpreter','latex');
colorbar; xlim([zL zR]);

figure; imagesc(z,r,Ne);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $n_e$ [m$^{-3}$]','Interpreter','latex');
colorbar; xlim([zL zR]);

figure; imagesc(z,r,Vpar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $v_\parallel$ [m/s]','Interpreter','latex');
colorbar; xlim([zL zR]);

figure;

imagesc(z, r, Mach);

set(gca,'YDir','normal','FontName','Times','FontSize',18);

xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');

title('ProtoLite $M_\parallel = v_\parallel/c_s$','Interpreter','latex');

colorbar;
caxis([-1 1]);
xlim([zL zR]);

hold on;

% -------------------------------------------------------------------------
% 1D lineout at r = 0.06 m
% -------------------------------------------------------------------------
rLine = 0.0;

[~, iRline] = min(abs(r - rLine));

% Show lineout location on 2D map
plot([z(1) z(end)], ...
     [r(iRline) r(iRline)], ...
     'w--', 'LineWidth', 2);

% Secondary axis for 1D lineout
yyaxis right

plot(z, Mach(iRline,:), ...
     'k-', 'LineWidth', 3);

ylabel(sprintf('$M_{\\parallel}(z)$ at $r=%.3f$ m', r(iRline)), ...
       'Interpreter','latex');

ylim([-1.05 1.05]);

yyaxis left

hold off;

figure; imagesc(z,r,GammaPar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $\Gamma_\parallel$ [m$^{-2}$ s$^{-1}$]','Interpreter','latex');
colorbar; xlim([zL zR]);

figure; imagesc(z,r,Qpar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $q_\parallel$ [W m$^{-2}$]','Interpreter','latex');
colorbar; xlim([zL zR]);

figure; imagesc(z,r,GammaPerp);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $\Gamma_\perp$ across $\psi$ [m$^{-2}$ s$^{-1}$]','Interpreter','latex');
colorbar; xlim([zL zR]);

figure; imagesc(z,r,Qperp);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $q_\perp$ across $\psi$ [W m$^{-2}$]','Interpreter','latex');
colorbar; xlim([zL zR]);

figure; imagesc(z,r,Qtot_psi);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('ProtoLite $q_{\mathrm{tot},\psi}$ [W m$^{-2}$]','Interpreter','latex');
colorbar; xlim([zL zR]);

% -------------------------------------------------------------------------
% WRITE TO NetCDF
% -------------------------------------------------------------------------
disp('Writing ProtoLite profiles to NetCDF...');

outFile = 'protoLite_profiles.nc';

if exist(outFile,'file')
    delete(outFile);
end

cmpr = {'Format','netcdf4','DeflateLevel',4,'Shuffle',true};

nccreate(outFile,'r','Dimensions',{'r',nR},'Datatype','double',cmpr{:});
nccreate(outFile,'z','Dimensions',{'z',nZ},'Datatype','double',cmpr{:});

ncwrite(outFile,'r',r);
ncwrite(outFile,'z',z);

vars = {
    'br',Br,'single';
    'bt',Bt,'single';
    'bz',Bz,'single';
    'psi',psi,'double';
    'psiN',psiN,'double';
    'ne',Ne,'double';
    'te',Te,'double';
    'vpar',Vpar,'double';
    'gamma_par',GammaPar,'double';
    'qpar',Qpar,'double';
    'gamma_par_r',GammaPar_r,'double';
    'gamma_par_t',GammaPar_t,'double';
    'gamma_par_z',GammaPar_z,'double';
    'gamma_par_psi',GammaPar_psi,'double';
    'bdotnpsi',bDotNpsi,'double';
    'qpar_r',Qpar_r,'double';
    'qpar_t',Qpar_t,'double';
    'qpar_z',Qpar_z,'double';
    'qpar_psi',Qpar_psi,'double';
    'gamma_perp',GammaPerp,'double';
    'gamma_perp_r',GammaPerp_r,'double';
    'gamma_perp_z',GammaPerp_z,'double';
    'qperp',Qperp,'double';
    'qperp_r',Qperp_r,'double';
    'qperp_z',Qperp_z,'double';
    'gamma_tot_r',GammaTot_r,'double';
    'gamma_tot_t',GammaTot_t,'double';
    'gamma_tot_z',GammaTot_z,'double';
    'gamma_tot_psi',GammaTot_psi,'double';
    'qtot_r',Qtot_r,'double';
    'qtot_t',Qtot_t,'double';
    'qtot_z',Qtot_z,'double';
    'qtot_psi',Qtot_psi,'double';
};

for ii = 1:size(vars,1)
    name = vars{ii,1};
    data = vars{ii,2};
    dtype = vars{ii,3};

    nccreate(outFile,name, ...
        'Dimensions',{'r',nR,'z',nZ}, ...
        'Datatype',dtype,cmpr{:});

    ncwrite(outFile,name,data);
end

% Absolute-value companions
fluxVarNames = {
    'gamma_par','qpar','gamma_par_r','gamma_par_t','gamma_par_z','gamma_par_psi', ...
    'qpar_r','qpar_t','qpar_z','qpar_psi', ...
    'gamma_perp','gamma_perp_r','gamma_perp_z', ...
    'qperp','qperp_r','qperp_z', ...
    'gamma_tot_r','gamma_tot_t','gamma_tot_z','gamma_tot_psi', ...
    'qtot_r','qtot_t','qtot_z','qtot_psi'
};

fluxData = {
    GammaPar,Qpar,GammaPar_r,GammaPar_t,GammaPar_z,GammaPar_psi, ...
    Qpar_r,Qpar_t,Qpar_z,Qpar_psi, ...
    GammaPerp,GammaPerp_r,GammaPerp_z, ...
    Qperp,Qperp_r,Qperp_z, ...
    GammaTot_r,GammaTot_t,GammaTot_z,GammaTot_psi, ...
    Qtot_r,Qtot_t,Qtot_z,Qtot_psi
};

for ii = 1:numel(fluxVarNames)
    absName = [fluxVarNames{ii}, '_abs'];

    nccreate(outFile,absName, ...
        'Dimensions',{'r',nR,'z',nZ}, ...
        'Datatype','double',cmpr{:});

    ncwrite(outFile,absName,abs(fluxData{ii}));
end

% -------------------------------------------------------------------------
% METADATA
% -------------------------------------------------------------------------
ncwriteatt(outFile,'/','title','ProtoLite plasma profiles with B-field, psi, psiN, ne, Te, velocity, particle fluxes, and heat fluxes');
ncwriteatt(outFile,'/','layout','All 2D variables are [r x z]');
ncwriteatt(outFile,'/','source_bfield_file',fileB);
ncwriteatt(outFile,'/','D_perp_m2_per_s',Dperp);
ncwriteatt(outFile,'/','chi_perp_m2_per_s',chiPerp);
ncwriteatt(outFile,'/','psi_normal_method','n_hat from analytic grad(psi): dpsi/dr=2*pi*r*Bz, dpsi/dz=-2*pi*r*Br');
ncwriteatt(outFile,'/','psi_norm_reference_z_m',z0);
ncwriteatt(outFile,'/','psi_norm_reference_r_m',r0);
ncwriteatt(outFile,'/','psi_norm_value_Wb',psi_norm_val);
ncwriteatt(outFile,'/','bDotNpsi_filter_tol',bDotN_tol);
ncwriteatt(outFile,'/','gasPuff_particles_per_s',gasPuff_Gps);
ncwriteatt(outFile,'/','PRF_kW',PRF_kW);
ncwriteatt(outFile,'/','Preq_kW_from_ionization_cost',Preq_kW);
ncwriteatt(outFile,'/','powerFactor_PRF_over_Preq_capped',powerFactor);
ncwriteatt(outFile,'/','Eion_eV_per_pair',Eion_eV);
ncwriteatt(outFile,'/','ne_ref_for_power_scaling_m3',ne_ref_for_power);
ncwriteatt(outFile,'/','Preq_kW_for_target_ne_max',Preq_ne_kW);
ncwriteatt(outFile,'/','power_margin_PRF_over_Preq_ne',powerMargin);
ncwriteatt(outFile,'/','ne_model_ne_cap_m3',ne_cap);
ncwriteatt(outFile,'/','ne_model_type','ProtoLite helicon fixed ne target with power check');

ncwriteatt(outFile,'r','units','m');
ncwriteatt(outFile,'z','units','m');
ncwriteatt(outFile,'br','units','T');
ncwriteatt(outFile,'bt','units','T');
ncwriteatt(outFile,'bz','units','T');
ncwriteatt(outFile,'psi','units','Wb');
ncwriteatt(outFile,'psiN','units','1');
ncwriteatt(outFile,'ne','units','m^-3');
ncwriteatt(outFile,'te','units','eV');
ncwriteatt(outFile,'vpar','units','m s^-1');

disp(['Wrote ', outFile, ' successfully.']);
disp('End of ProtoLite profile generation script.');