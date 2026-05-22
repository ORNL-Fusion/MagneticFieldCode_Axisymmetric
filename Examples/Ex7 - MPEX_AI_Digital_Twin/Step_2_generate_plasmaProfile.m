%% ProtoMPEX: Compute ne, te using original psi formulation and write to one NetCDF
% + Adds perpendicular (cross-field) particle flux using D_perp = 0.45 m^2/s
%   - Vector components: gamma_perp_r, gamma_perp_z = -D_perp * grad(ne)
%   - Scalar across-psi flux: gamma_perp = gamma_perp_vec · n_hat(psi)
%
% + Adds perpendicular (cross-field) HEAT flux with chi_perp = 1.0 m^2/s
%   - Vector components: qperp_r, qperp_z = -chi_perp * ne * e * grad(Te[eV])
%   - Scalar across-psi flux: qperp = qperp_vec · n_hat(psi)
%
% KEY CHANGE to remove striping:
%   - Compute grad(psi) analytically from B: dpsi/dr = 2*pi*r*Bz, dpsi/dz = -2*pi*r*Br
%     instead of using gradient(psi).
%
% + Writes psi and psiN to the NetCDF file.

close all; clear; clc;

% -------------------------------------------------------------------------
% READ EXISTING B-FIELD FILE (assumed [r × z])
% -------------------------------------------------------------------------
fileB = 'bfield_protoMPEX.nc';
r  = ncread(fileB,'r');                % [m]
z  = ncread(fileB,'z');                % [m]
Br = ncread(fileB,'br');               % [r × z]
Bt = ncread(fileB,'bt');               % [r × z]
Bz = ncread(fileB,'bz');               % [r × z]

% If orientation is swapped, fix once
if ~isequal(size(Bz), [numel(r) numel(z)])
    Br = permute(Br,[2 1]);
    Bt = permute(Bt,[2 1]);
    Bz = permute(Bz,[2 1]);
end

nR = numel(r);
nZ = numel(z);

% -------------------------------------------------------------------------
% PLOT INPUT Bz FIELD
% -------------------------------------------------------------------------
figure; imagesc(z, r, Bz);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('Input $B_z$ [T]','Interpreter','latex');
colorbar; axis image tight; xlim([0.5 4.2]);

% -------------------------------------------------------------------------
% CALCULATE PSI USING ORIGINAL METHOD (YOUR LOOP)
% -------------------------------------------------------------------------
disp('Calculating psi using original summation method...');
dx = r(2) - r(1);
psi = zeros(nR, nZ);
for ii = 1:nR
    for jj = 1:nZ
        psi(ii,jj) = 2*pi * sum(Bz(1:ii,jj) .* r(1:ii)) * dx;
    end
end

% Normalize psi using helicon / stagnation location (z=1.745, r=0.0625)
z0 = 1.745;         % [m]
r0 = 0.0625;        % [m]

psi_norm_val = interp2(z, r, psi, z0, r0, 'linear');
psiN = psi ./ psi_norm_val;

figure; imagesc(z, r, psiN);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('Normalized $\psi_N$','Interpreter','latex');
colorbar;

% -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES FROM psiN
% -------------------------------------------------------------------------
disp('Constructing Te and ne profiles...');
qe  = 1.602176634e-19;      % [C]

% --------- Te: profile definition ----------
Te0       = 2;      % eV at core (psiN=0)
Te_edge   = 8;      % eV at LCFS (psiN=1)
Te_floor  = 0.2;    % eV far outside
alpha     = 1.75;   % core shaping
lambdaPsi = 0.15;   % SOL falloff width (smaller = steeper)

Te = zeros(size(psiN));

% Core (psiN <= 1): Te(0)=Te0, Te(1)=Te_edge
Te(psiN <= 1) = Te_edge + (Te0 - Te_edge) .* (1 - psiN(psiN <= 1)).^alpha;

% SOL (psiN > 1): exponential starting from Te_edge at psiN=1
Te(psiN > 1)  = Te_floor + (Te_edge - Te_floor) .* exp(-(psiN(psiN > 1) - 1) ./ lambdaPsi);

% Density model (helicon-only): fixed ne_max target + required-power check
gasPuff_Gps = 1.0e20;   % gas puff/source rate [particles/s]
PRF_kW      = 180.0;    % coupled RF power [kW]

Eion_eV     = 1.0e3;    % ionization cost [eV/pair]
ne_cap      = 1E20;     % target max density [m^-3]

% RF power required to sustain source rate G at the ionization cost
Preq_kW = gasPuff_Gps * Eion_eV * qe / 1e3;
powerFactor = min(PRF_kW / max(Preq_kW,1e-12), 1.0);  % saturates at 1

% User-fixed density target
ne_min = 1.0E16;  % m^-3
ne_max = ne_cap;  % user-fixed target
ne_max = min(max(ne_max, 5.0E18), 5.0E20);  % safety bounds

% Required RF power to sustain chosen ne_max (assumption: linear with ne)
ne_ref_for_power = 1.0E20;   % [m^-3] reference density for Preq scaling
Preq_ne_kW = Preq_kW * (ne_max / ne_ref_for_power);
powerMargin = PRF_kW / max(Preq_ne_kW,1e-12);
disp(sprintf(['Helicon ne_max = %.3e m^-3 (fixed source; G=%.2e 1/s)\n', ...
    'Required power for ne_{max}: %.2f kW (available PRF=%.1f kW)'], ...
    ne_max, gasPuff_Gps, Preq_ne_kW, PRF_kW));

Ne = (psiN < 1) .* (ne_max - ne_min) .* (1 - psiN).^1.75 + ne_min;

% -------------------------------------------------------------------------
% BUILD PARALLEL VELOCITY PROFILE FROM Te
% Mach is piecewise in z
% -------------------------------------------------------------------------
disp('Constructing parallel velocity profile from Te...');

mp  = 1.67262192369e-27;    % [kg] proton mass (assume H+)

Cs = sqrt(qe .* Te ./ mp);  % [nR × nZ] ion sound speed

% Piecewise Mach profile in z
zL = z(1);
zR = z(end);

% Control points (m) and jump levels
zDropL = min(max(0.5, zL), zR);  % 0 -> -1 step location
zJump  = min(max(1.8, zL), zR);  % small jump location
zDropR = min(max(4.2, zL), zR);  % +1 -> 0 step location
Mpre   = -0.10;                  % value just before zJump
Mpost  = +0.10;                  % value just after zJump

Mz = zeros(1, nZ);

% Region A: left plateau (M=0)
idxA = z < zDropL;
Mz(idxA) = 0;

% Region B: ramp from -1 at zDropL to Mpre at zJump
idxB = (z >= zDropL) & (z < zJump);
if zJump > zDropL
    Mz(idxB) = -1 + (Mpre + 1) .* (z(idxB) - zDropL) ./ (zJump - zDropL);
else
    Mz(idxB) = -1;
end

% Region C: ramp from Mpost at zJump to +1 at zDropR
idxC = (z >= zJump) & (z < zDropR);
if zDropR > zJump
    Mz(idxC) = Mpost + (1 - Mpost) .* (z(idxC) - zJump) ./ (zDropR - zJump);
else
    Mz(idxC) = 1;
end

% Region D: right plateau (M=0)
idxD = z >= zDropR;
Mz(idxD) = 0;

% Clip for numerical safety
Mz = max(min(Mz, 1.0), -1.0);

% Expand to [nR × nZ] and compute parallel velocity
M = repmat(Mz, nR, 1);
Vpar = M .* Cs;

% For plotting/checks
Mach = Vpar ./ Cs;   % == M

% -------------------------------------------------------------------------
% ANALYTICAL PARTICLE AND HEAT FLUXES (PARALLEL)
% -------------------------------------------------------------------------
disp('Constructing analytical particle and heat fluxes...');
GammaPar = Ne .* Vpar;
Qpar     = (5/2) * qe .* Te .* GammaPar;

% Unit vector along magnetic field (for directional flux components)
Bmag = sqrt(Br.^2 + Bt.^2 + Bz.^2);
bvec = Bmag(:); 
bvec = bvec(isfinite(bvec) & bvec>0);
if isempty(bvec)
    b_floor = 1e-12;
else
    b_floor = max(prctile(bvec,2), 1e-12);
end
Bmag_safe = max(Bmag, b_floor);

bhat_r = Br ./ Bmag_safe;
bhat_t = Bt ./ Bmag_safe;
bhat_z = Bz ./ Bmag_safe;

% Parallel particle flux vector components
GammaPar_r = GammaPar .* bhat_r;
GammaPar_t = GammaPar .* bhat_t;
GammaPar_z = GammaPar .* bhat_z;

% -------------------------------------------------------------------------
% PERPENDICULAR PARTICLE FLUX (DIFFUSIVE) WITH D_perp = 0.45
% + Stripe-free across-psi projection using analytic grad(psi) from B
% -------------------------------------------------------------------------
disp('Computing perpendicular particle flux with D_perp = 0.45 m^2/s ...');

Dperp = 0.45;

dr = r(2) - r(1);
dz = z(2) - z(1);

% Diffusive particle flux vector components
[dndr, dndz] = gradient(Ne, dr, dz);
GammaPerp_r = -Dperp .* dndr;   % [m^-2 s^-1]
GammaPerp_z = -Dperp .* dndz;   % [m^-2 s^-1]

% Compute grad(psi) from B-field (axisymmetric identity)
rMat = repmat(r(:), 1, nZ);

dpsidr = 2*pi .* rMat .* Bz;      % ∂ψ/∂r
dpsidz = -2*pi .* rMat .* Br;     % ∂ψ/∂z

gradPsiMag = sqrt(dpsidr.^2 + dpsidz.^2);

% Prevent blow-up where gradPsi is tiny
gvec = gradPsiMag(:); 
gvec = gvec(isfinite(gvec) & gvec>0);
if isempty(gvec)
    eps_floor = 1e-12;
else
    eps_floor = prctile(gvec, 2);
    eps_floor = max(eps_floor, 1e-12);
end
gradPsiMag_safe = max(gradPsiMag, eps_floor);

nr = dpsidr ./ gradPsiMag_safe;
nz = dpsidz ./ gradPsiMag_safe;

% Mask truly tiny gradients
mask_bad = gradPsiMag < eps_floor;
nr(mask_bad) = 0; 
nz(mask_bad) = 0;

% Scalar particle flux across psi surfaces
GammaPerp = GammaPerp_r .* nr + GammaPerp_z .* nz;
GammaPerp(~isfinite(GammaPerp)) = 0;

% -------------------------------------------------------------------------
% PERPENDICULAR HEAT FLUX (CONDUCTIVE) WITH chi_perp = 1.0
% + Across-psi projection using same n_hat(psi) (nr, nz)
% -------------------------------------------------------------------------
disp('Computing perpendicular heat flux with chi_perp = 1.0 m^2/s ...');

chiPerp = 1.0;  % [m^2/s]

% Temperature gradient (Te is in eV, convert gradient to J/m by multiplying qe)
[dTdr, dTdz] = gradient(Te, dr, dz);     % [eV/m]

% Conductive heat flux vector components:
% q_perp_vec = -chi_perp * ne * grad(Te*qe)
Qperp_r = -chiPerp .* Ne .* (qe .* dTdr);    % [W/m^2]
Qperp_z = -chiPerp .* Ne .* (qe .* dTdz);    % [W/m^2]

% Scalar heat flux across psi surfaces
Qperp = Qperp_r .* nr + Qperp_z .* nz;       % [W/m^2]
Qperp(~isfinite(Qperp)) = 0;

% Parallel heat flux vector components
Qpar_r = Qpar .* bhat_r;
Qpar_t = Qpar .* bhat_t;
Qpar_z = Qpar .* bhat_z;

% Psi-normal projections of parallel fluxes are ideally ~0.
% Non-zero values here are numerical leakage from b_hat·n_hat_psi ~= 0.
bDotNpsi = bhat_r .* nr + bhat_z .* nz;
bDotN_tol = 5e-3;
bDotNpsi_f = bDotNpsi;
bDotNpsi_f(abs(bDotNpsi_f) < bDotN_tol) = 0;

GammaPar_psi = GammaPar .* bDotNpsi_f;
GammaPar_psi(~isfinite(GammaPar_psi)) = 0;

Qpar_psi = Qpar .* bDotNpsi_f;
Qpar_psi(~isfinite(Qpar_psi)) = 0;

% Total particle and heat flux vectors/scalars
GammaTot_r = GammaPar_r + GammaPerp_r;
GammaTot_t = GammaPar_t;  % perpendicular model has no toroidal component
GammaTot_z = GammaPar_z + GammaPerp_z;

% Physics-facing across-psi total uses perpendicular contribution.
% Parallel cross-psi transport is zero in ideal axisymmetry.
GammaTot_psi = GammaPerp;
GammaTot_psi(~isfinite(GammaTot_psi)) = 0;

Qtot_r = Qpar_r + Qperp_r;
Qtot_t = Qpar_t;          % perpendicular model has no toroidal component
Qtot_z = Qpar_z + Qperp_z;
Qtot_psi = Qperp;
Qtot_psi(~isfinite(Qtot_psi)) = 0;

% -------------------------------------------------------------------------
% PLOT PROFILES
% -------------------------------------------------------------------------
figure; imagesc(z, r, Te);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$T_e$ [eV]','Interpreter','latex');
colorbar;

figure; imagesc(z, r, Ne);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$n_e$ [m$^{-3}$]','Interpreter','latex');
colorbar;

figure; imagesc(z, r, Vpar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$v_\parallel$ [m/s]','Interpreter','latex');
colorbar;

figure; imagesc(z, r, Mach);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$M_\parallel = v_\parallel / c_s$','Interpreter','latex');
colorbar;
caxis([-1 1]);
hold on;

% Superpose Mach lineout at r = 0 (nearest grid point)
[~, iR0] = min(abs(r - 0.0));
plot([z(1) z(end)], [r(iR0) r(iR0)], 'w--', 'LineWidth', 1.5);

yyaxis right;
plot(z, Mach(iR0,:), 'k-', 'LineWidth', 2.0);
ylabel(sprintf('Lineout $M_{\\parallel}(z)$ at $r=%.4f$ m', r(iR0)), 'Interpreter','latex');
ylim([-1.05 1.05]);
yyaxis left;

figure; imagesc(z, r, GammaPar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$\Gamma_\parallel$ [m$^{-2}$ s$^{-1}$]','Interpreter','latex');
colorbar;

figure; imagesc(z, r, Qpar);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$q_\parallel$ [W m$^{-2}$]','Interpreter','latex');
colorbar;

figure; imagesc(z, r, GammaPerp);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$\Gamma_\perp$ across $\psi$ [m$^{-2}$ s$^{-1}$]','Interpreter','latex');
colorbar;

figure; imagesc(z, r, Qperp);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$q_\perp$ across $\psi$ [W m$^{-2}$] ($\chi_\perp=1$)','Interpreter','latex');
colorbar;

Qtot_map = Qpar + Qperp;
figure; imagesc(z, r, Qtot_map);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$q_{\mathrm{tot}} = q_{\parallel} + q_{\perp}$ [W m$^{-2}$]','Interpreter','latex');
colorbar;

figure; imagesc(z, r, Qtot_psi);
set(gca,'YDir','normal','FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('$r$ [m]','Interpreter','latex');
title('$q_{\mathrm{tot},\psi}$ (across $\psi$) [W m$^{-2}$]','Interpreter','latex');
colorbar;

% -------------------------------------------------------------------------
% HEAT-FLUX PROFILES
%   1) Lineout vs z at r = 0.06 m
%   2) Radial profile averaged over 1.6 < z < 1.9
% -------------------------------------------------------------------------
[~, iR06] = min(abs(r - 0.06));
zMask = (z > 1.6) & (z < 1.9);
if ~any(zMask)
    error('No z points satisfy 1.6 < z < 1.9. Please check z-grid.');
end

qpar_z_r06  = abs(Qpar(iR06, :));
qperp_z_r06 = abs(Qperp(iR06, :));
qtot_z_r06  = abs(Qpar(iR06, :) + Qperp(iR06, :));
qpsi_z_r06  = abs(Qtot_psi(iR06, :));

qpar_r_win  = mean(abs(Qpar(:, zMask)), 2);
qperp_r_win = mean(abs(Qperp(:, zMask)), 2);
qtot_r_win  = mean(abs(Qpar(:, zMask) + Qperp(:, zMask)), 2);
qpsi_r_win  = mean(abs(Qtot_psi(:, zMask)), 2);

% Positive floor for log-scale plotting
qFloor = 1e-12;
qpar_z_r06  = max(qpar_z_r06,  qFloor);
qperp_z_r06 = max(qperp_z_r06, qFloor);
qtot_z_r06  = max(qtot_z_r06,  qFloor);
qpsi_z_r06  = max(qpsi_z_r06,  qFloor);
qpar_r_win  = max(qpar_r_win,  qFloor);
qperp_r_win = max(qperp_r_win, qFloor);
qtot_r_win  = max(qtot_r_win,  qFloor);
qpsi_r_win  = max(qpsi_r_win,  qFloor);

figure;
hPar1 = semilogy(z, qpar_z_r06, 'k-', 'LineWidth', 2.0); hold on;
hPerp1 = semilogy(z, qperp_z_r06, 'r-', 'LineWidth', 2.0);
hTot1  = semilogy(z, qtot_z_r06, 'b-', 'LineWidth', 2.0);
hPsi1  = semilogy(z, qpsi_z_r06, 'm--', 'LineWidth', 1.8);
xline(1.6, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(1.9, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
set(gca,'FontName','Times','FontSize',18);
xlabel('$z$ [m]','Interpreter','latex');
ylabel('Heat flux [W m$^{-2}$]','Interpreter','latex');
ylim([1e0 1e7]);
title(sprintf('Absolute heat-flux lineout at $r=%.4f$ m', r(iR06)), 'Interpreter','latex');
leg1 = {'$|q_\parallel|$','$|q_\perp|$','$|q_{\mathrm{tot}}|$','$|q_{\mathrm{tot},\psi}|$'};
legend([hPar1 hPerp1 hTot1 hPsi1], leg1, 'Interpreter','latex','Location','best');
grid on; box on;

figure;
hPar2 = semilogy(r, qpar_r_win, 'k-', 'LineWidth', 2.0); hold on;
hPerp2 = semilogy(r, qperp_r_win, 'r-', 'LineWidth', 2.0);
hTot2  = semilogy(r, qtot_r_win, 'b-', 'LineWidth', 2.0);
hPsi2  = semilogy(r, qpsi_r_win, 'm--', 'LineWidth', 1.8);
xline(0.06, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
set(gca,'FontName','Times','FontSize',18);
xlabel('$r$ [m]','Interpreter','latex');
ylabel('Heat flux [W m$^{-2}$]','Interpreter','latex');
ylim([1e0 1e7]);
title('Absolute radial heat-flux profile averaged over $1.6<z<1.9$', 'Interpreter','latex');
leg2 = {'$|q_\parallel|$','$|q_\perp|$','$|q_{\mathrm{tot}}|$','$|q_{\mathrm{tot},\psi}|$'};
legend([hPar2 hPerp2 hTot2 hPsi2], leg2, 'Interpreter','latex','Location','best');
grid on; box on;

% -------------------------------------------------------------------------
% WRITE TO ONE NetCDF FILE
% -------------------------------------------------------------------------
disp('Writing combined ProtoMPEX profiles to NetCDF...');

outFile = 'protoMPEX_profiles.nc';
if exist(outFile,'file'); delete(outFile); end

cmpr = {'Format','netcdf4','DeflateLevel',4,'Shuffle',true};

% Dimensions
nccreate(outFile,'r','Dimensions',{'r',nR},'Datatype','double',cmpr{:});
nccreate(outFile,'z','Dimensions',{'z',nZ},'Datatype','double',cmpr{:});
ncwrite(outFile,'r',r);
ncwrite(outFile,'z',z);

% Variables (all [r × z])
nccreate(outFile,'br','Dimensions',{'r',nR,'z',nZ},'Datatype','single',cmpr{:});
nccreate(outFile,'bt','Dimensions',{'r',nR,'z',nZ},'Datatype','single',cmpr{:});
nccreate(outFile,'bz','Dimensions',{'r',nR,'z',nZ},'Datatype','single',cmpr{:});
nccreate(outFile,'psi','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'psiN','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'ne','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'te','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'vpar','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_par','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qpar','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_par_r','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_par_t','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_par_z','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_par_psi','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'bdotnpsi','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qpar_r','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qpar_t','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qpar_z','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qpar_psi','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});

% Perpendicular particle fluxes
nccreate(outFile,'gamma_perp','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_perp_r','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_perp_z','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});

% Perpendicular heat fluxes
nccreate(outFile,'qperp','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qperp_r','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qperp_z','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});

% Total particle and heat fluxes
nccreate(outFile,'gamma_tot_r','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_tot_t','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_tot_z','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'gamma_tot_psi','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qtot_r','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qtot_t','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qtot_z','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
nccreate(outFile,'qtot_psi','Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});

% Absolute-value companions for all flux variables (suffix: _abs)
fluxVarNames = {'gamma_par','qpar','gamma_par_r','gamma_par_t','gamma_par_z','gamma_par_psi', ...
                'qpar_r','qpar_t','qpar_z','qpar_psi', ...
                'gamma_perp','gamma_perp_r','gamma_perp_z', ...
                'qperp','qperp_r','qperp_z', ...
                'gamma_tot_r','gamma_tot_t','gamma_tot_z','gamma_tot_psi', ...
                'qtot_r','qtot_t','qtot_z','qtot_psi'};
for iFlux = 1:numel(fluxVarNames)
    nccreate(outFile,[fluxVarNames{iFlux}, '_abs'], ...
        'Dimensions',{'r',nR,'z',nZ},'Datatype','double',cmpr{:});
end

% Write data
ncwrite(outFile,'br',Br);
ncwrite(outFile,'bt',Bt);
ncwrite(outFile,'bz',Bz);
ncwrite(outFile,'psi',psi);
ncwrite(outFile,'psiN',psiN);
ncwrite(outFile,'ne',Ne);
ncwrite(outFile,'te',Te);
ncwrite(outFile,'vpar',Vpar);
ncwrite(outFile,'gamma_par',GammaPar);
ncwrite(outFile,'qpar',Qpar);
ncwrite(outFile,'gamma_par_r',GammaPar_r);
ncwrite(outFile,'gamma_par_t',GammaPar_t);
ncwrite(outFile,'gamma_par_z',GammaPar_z);
ncwrite(outFile,'gamma_par_psi',GammaPar_psi);
ncwrite(outFile,'bdotnpsi',bDotNpsi);
ncwrite(outFile,'qpar_r',Qpar_r);
ncwrite(outFile,'qpar_t',Qpar_t);
ncwrite(outFile,'qpar_z',Qpar_z);
ncwrite(outFile,'qpar_psi',Qpar_psi);

ncwrite(outFile,'gamma_perp',GammaPerp);
ncwrite(outFile,'gamma_perp_r',GammaPerp_r);
ncwrite(outFile,'gamma_perp_z',GammaPerp_z);

ncwrite(outFile,'qperp',Qperp);
ncwrite(outFile,'qperp_r',Qperp_r);
ncwrite(outFile,'qperp_z',Qperp_z);

ncwrite(outFile,'gamma_tot_r',GammaTot_r);
ncwrite(outFile,'gamma_tot_t',GammaTot_t);
ncwrite(outFile,'gamma_tot_z',GammaTot_z);
ncwrite(outFile,'gamma_tot_psi',GammaTot_psi);
ncwrite(outFile,'qtot_r',Qtot_r);
ncwrite(outFile,'qtot_t',Qtot_t);
ncwrite(outFile,'qtot_z',Qtot_z);
ncwrite(outFile,'qtot_psi',Qtot_psi);

% Write absolute-value companions
fluxData = {GammaPar,Qpar,GammaPar_r,GammaPar_t,GammaPar_z,GammaPar_psi, ...
            Qpar_r,Qpar_t,Qpar_z,Qpar_psi, ...
            GammaPerp,GammaPerp_r,GammaPerp_z, ...
            Qperp,Qperp_r,Qperp_z, ...
            GammaTot_r,GammaTot_t,GammaTot_z,GammaTot_psi, ...
            Qtot_r,Qtot_t,Qtot_z,Qtot_psi};
for iFlux = 1:numel(fluxVarNames)
    ncwrite(outFile,[fluxVarNames{iFlux}, '_abs'], abs(fluxData{iFlux}));
end

% Metadata
ncwriteatt(outFile,'/','title','ProtoMPEX profiles with ne, Te, psi, psiN, B-field, v_parallel, analytical fluxes, and perpendicular particle + heat fluxes');
ncwriteatt(outFile,'/','layout','All 2D variables are [r × z]');
ncwriteatt(outFile,'/','D_perp_m2_per_s',Dperp);
ncwriteatt(outFile,'/','chi_perp_m2_per_s',chiPerp);
ncwriteatt(outFile,'/','psi_normal_method','n_hat from analytic grad(psi): dpsi/dr=2*pi*r*Bz, dpsi/dz=-2*pi*r*Br');
ncwriteatt(outFile,'/','psi_norm_reference_z_m',z0);
ncwriteatt(outFile,'/','psi_norm_reference_r_m',r0);
ncwriteatt(outFile,'/','psi_norm_value_Wb',psi_norm_val);
ncwriteatt(outFile,'/','bDotNpsi_filter_tol',bDotN_tol);
ncwriteatt(outFile,'/','flux_convention','Signed flux variables are primary; companion *_abs variables store magnitudes.');
ncwriteatt(outFile,'/','gasPuff_particles_per_s',gasPuff_Gps);
ncwriteatt(outFile,'/','PRF_kW',PRF_kW);
ncwriteatt(outFile,'/','Preq_kW_from_ionization_cost',Preq_kW);
ncwriteatt(outFile,'/','powerFactor_PRF_over_Preq_capped',powerFactor);
ncwriteatt(outFile,'/','Eion_eV_per_pair',Eion_eV);
ncwriteatt(outFile,'/','ne_ref_for_power_scaling_m3',ne_ref_for_power);
ncwriteatt(outFile,'/','Preq_kW_for_target_ne_max',Preq_ne_kW);
ncwriteatt(outFile,'/','power_margin_PRF_over_Preq_ne',powerMargin);
ncwriteatt(outFile,'/','ne_model_ne_cap_m3',ne_cap);
ncwriteatt(outFile,'/','ne_model_type','helicon_fixed_ne_target_with_power_check');

ncwriteatt(outFile,'r','units','m');
ncwriteatt(outFile,'z','units','m');
ncwriteatt(outFile,'br','units','tesla');
ncwriteatt(outFile,'bt','units','tesla');
ncwriteatt(outFile,'bz','units','tesla');

ncwriteatt(outFile,'psi','units','Wb');
ncwriteatt(outFile,'psi','description','Poloidal flux computed from original summation method: psi(r,z)=2*pi*sum(Bz(1:r,z).*r)*dr');

ncwriteatt(outFile,'psiN','units','1');
ncwriteatt(outFile,'psiN','description','Normalized poloidal flux psi/psi_ref');

ncwriteatt(outFile,'ne','units','m^-3');
ncwriteatt(outFile,'te','units','eV');
ncwriteatt(outFile,'vpar','units','m s^-1');
ncwriteatt(outFile,'gamma_par','units','m^-2 s^-1');
ncwriteatt(outFile,'qpar','units','W m^-2');
ncwriteatt(outFile,'gamma_par_r','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_par_t','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_par_z','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_par_psi','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_par_psi','description','Parallel particle flux projected to psi-normal; non-zero is numerical leakage due to b_hat·n_hat_psi != 0');

ncwriteatt(outFile,'bdotnpsi','units','1');
ncwriteatt(outFile,'bdotnpsi','description','Dot product b_hat·n_hat_psi; should be approximately zero');

ncwriteatt(outFile,'qpar_r','units','W m^-2');
ncwriteatt(outFile,'qpar_t','units','W m^-2');
ncwriteatt(outFile,'qpar_z','units','W m^-2');
ncwriteatt(outFile,'qpar_psi','units','W m^-2');
ncwriteatt(outFile,'qpar_psi','description','Parallel heat flux projected to psi-normal; non-zero is numerical leakage due to b_hat·n_hat_psi != 0');

ncwriteatt(outFile,'gamma_perp','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_perp','description','Scalar cross-field particle flux across psi: (-D_perp*grad ne)·n_hat(psi), with n_hat from B identities');
ncwriteatt(outFile,'gamma_perp_r','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_perp_z','units','m^-2 s^-1');

ncwriteatt(outFile,'qperp','units','W m^-2');
ncwriteatt(outFile,'qperp','description','Scalar cross-field heat flux across psi: (-chi_perp*ne*e*grad Te)·n_hat(psi), with n_hat from B identities');
ncwriteatt(outFile,'qperp_r','units','W m^-2');
ncwriteatt(outFile,'qperp_z','units','W m^-2');

ncwriteatt(outFile,'gamma_tot_r','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_tot_t','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_tot_z','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_tot_psi','units','m^-2 s^-1');
ncwriteatt(outFile,'gamma_tot_psi','description','Across-psi total particle flux; set to gamma_perp (parallel cross-psi is ideally zero in axisymmetry)');

ncwriteatt(outFile,'qtot_r','units','W m^-2');
ncwriteatt(outFile,'qtot_t','units','W m^-2');
ncwriteatt(outFile,'qtot_z','units','W m^-2');
ncwriteatt(outFile,'qtot_psi','units','W m^-2');
ncwriteatt(outFile,'qtot_psi','description','Across-psi total heat flux; set to qperp (parallel cross-psi is ideally zero in axisymmetry)');

% Metadata for absolute-value companions
for iFlux = 1:numel(fluxVarNames)
    absName = [fluxVarNames{iFlux}, '_abs'];
    ncwriteatt(outFile,absName,'units',ncreadatt(outFile,fluxVarNames{iFlux},'units'));
    ncwriteatt(outFile,absName,'description',['Absolute value of ', fluxVarNames{iFlux}]);
end

disp(['✅ Wrote ', outFile, ' successfully.']);