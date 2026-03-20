function [r_h,z_h,r_t,z_t,Q_tot_h,Q_tot_t,r_bz0,z_bz0,bz0,Qtot_map] = ObjectiveOracle(parms,path_CoilSetup,do_vis,h_fig,parms_tune)


% Magnetic configuration of interest:
confType = 'conf_G';

% =========================================================================
% Read "CoilSetup" spreadsheet:
dum1 = tic;
disp('Reading "coilSetup" spreadsheet...')
coilSetup = readtable('CoilSetup_ProtoMPEX.xlsx','Sheet',confType);
disp('Reading complete!')
toc(dum1);


z_Dump = 0.5;
z_Target = 4.2;
r1D = linspace(1e-6,0.1  ,400 );
z1D = linspace(z_Dump,z_Target,500);

coilCurrents{1}.TR1 = parms(1);
coilCurrents{1}.TR2 = parms(2);
coilCurrents{1}.PS1 = parms(3);
coilCurrents{1}.PS2 = parms(4);
coilCurrents{1}.PS3 = parms(5);

dum1 = tic;
disp(['Calculating magnetic field for ',num2str(numel(coilCurrents)),' cases ...'])
for ii = 1:numel(coilCurrents)
    % Create "coil" structure":
    [coil] = CreateCoilStructure(coilSetup,coilCurrents{ii});
    % Calculate magnetic field and vector potential:
    [Br2D,Bz2D,~,phi2D{ii},z2D,r2D] = CalculateMagField(coil,z1D,r1D,'grid');
    % Magnetic field magnitude:
    B2D{ii} = sqrt(Br2D.*Br2D  + Bz2D.*Bz2D);
    disp(['Case ',num2str(ii),' complete!'])
end
disp(['Complete! Elapsed time: ',num2str(toc(dum1)),' s'])
clearvars dum*

%% Draw Proto-MPEX vacuum vessel:
if ~strcmpi(confType,'conf_G')
    error('Change "confType" to config_G')
end
DrawVacuumVessel_Conf_G

%% SECTION 3: Calculate reference flux:
% =========================================================================
% Define reference flux based on phi at vacuum vessel boundary
for ii = 1:numel(coilCurrents)
    % Select region of interest:
    rng_jj = find(vessel_1.z > 0.5 & vessel_1.z < 3.7);
    % Initialize variable:
    phiBoundary = ones(size(vessel_1.z));
    % Interpolate phi along vaccum vessel contour
    zq = vessel_1.z(rng_jj);
    rq = vessel_1.r(rng_jj);
    a = interp2(z1D,r1D,phi2D{ii}',zq,rq);
    phiBoundary(rng_jj) = a;
    % Find location of minimim phi along contour
    [~,jj] = min(phiBoundary);
    % Physical location of limit:
    rlimit(ii) = vessel_1.r(jj);
    zlimit(ii) = vessel_1.z(jj);
    % Extract reference magnetic flux at the limiting location:
    nr = find(r1D > rlimit(ii),1,'first');
    nz = find(z1D > zlimit(ii),1);
    phi0 = interp2(z1D,r1D,phi2D{ii}',zlimit(ii),rlimit(ii));
    % Flux coordinate
    xi{ii} = phi2D{ii}/phi0;
end
%% SECTION 4: MAGNETIC FIELD LINES AND PLASMA EDGE
% =========================================================================
% Magnetic field field lines up to the plasma edge
for ii = 1:numel(coilCurrents)
    % Define the number of flux lines to plot:
    xi_lines = linspace(1e-2,1,1);
    % Calculate the flux line trajectory r(z) in physical space:
    for jj = 1:numel(xi_lines)
        
        C = contour(z2D,r2D,xi{ii},[1,1]*xi_lines(jj));
        %C = contourc(r2D(1,1:end)',z2D(1:end,1),xi{ii},[1,1]*xi_lines(jj));
        z_fluxline{ii}{jj} = C(1,2:end);
        r_fluxline{ii}{jj} = C(2,2:end);
    end
end
if do_vis == 1
    close(gcf)
end
clearvars ii


%% SECTION 5: Plot data

% Magnetic flux mapping:
% =========================================================================
if do_vis == 1
    figure(h_fig)
    subplot(4,2,3:4)
    hold on
    % Magnetic coils:
    for ii = 1:numel(coil)
        plot(coil{ii}.zfil,+coil{ii}.rfil,'r.');
        plot(coil{ii}.zfil,-coil{ii}.rfil,'r.');
    end
    % Flux lines:
    lineColor = {'k','r','g','bl','m','c'};
    for ii = 1:numel(coilCurrents)
        for jj = 1:1:numel(xi_lines)
            dum1 = plot((z_fluxline{ii}{jj}),+r_fluxline{ii}{jj},lineColor{ii});
            dum2 = plot((z_fluxline{ii}{jj}),-r_fluxline{ii}{jj},lineColor{ii});
            if jj == numel(xi_lines)
                set(dum1,'LineStyle','-','LineWidth',2)
                set(dum2,'LineStyle','-','LineWidth',2)
            end
        end
    end

    % Vacuuum vessel
    plot(vessel_1.z,+vessel_1.r,'r','LineWidth',1)
    plot(vessel_1.z,-vessel_1.r,'r','LineWidth',1)
    plot(vessel_0.z,+vessel_0.r,'k-','LineWidth',1)
    plot(vessel_0.z,-vessel_0.r,'k-','LineWidth',1)

    % Limiting location:
    for ii = 1:numel(coilCurrents)
        hdum1 = plot(zlimit(ii),+rlimit(ii));
        hdum2 = plot(zlimit(ii),-rlimit(ii));
        set(hdum1,'Marker','o','MarkerFaceColor',lineColor{ii},'Color',lineColor{ii})
        set(hdum2,'Marker','o','MarkerFaceColor',lineColor{ii},'Color',lineColor{ii})
    end

    % Target:
    hT = line(z_Target*[1,1],0.045*[-1,+1]);
    set(hT,'color','k','LineWidth',4)

    % Formatting:
    set(gca,'FontName','times')
    zoomType = 3;
    switch zoomType
        case 1
            set(gca,'PlotBoxAspectRatio',[2 1.5 1])
            xlim([0.25,4.5])
            ylim(0.15*[0,+1])
        case 2
            set(gca,'PlotBoxAspectRatio',[2.3 1.5 1])
            xlim([1,2.5])
            ylim(0.1*[-1,+1])
        case 3
            set(gca,'PlotBoxAspectRatio',[4 1 1])
            xlim([0,4.5])
            ylim(0.35*[-1,+1])
    end
    box on
    % xlabel('z [m]','Interpreter','Latex','FontSize',13)
    ylabel('r [m]','Interpreter','Latex','FontSize',13)
    grid on
end

% Magnetic field profiles:
% =========================================================================
% figureName = 'ProtoMPEX_BzProfileVarious';
% figure('color','w','Tag',figureName)

r  = r1D(:);                 % nR×1
z  = z1D(:);                 % nZ×1
Br = Br2D;                   % likely nZ×nR
Bz = Bz2D;                   % likely nZ×nR
Bt = zeros(size(Bz2D),'like',Bz2D);

r_bz0 = r;
z_bz0 = z;
bz0 = Bz;
if do_vis == 1
    subplot(4,2,5:6)
    for ii = 1:numel(coilCurrents)
        hBz(ii) = plot(z1D,B2D{ii}(:,1),lineColor{ii},'LineWidth',2);
    end
    box on 
    grid on
    set(gca,'PlotBoxAspectRatio',[4 1 1])
    set(gca,'FontName','times')
    xlabel('z [m]','Interpreter','Latex','FontSize',13)
    ylabel('B$_0$ [T]','Interpreter','Latex','FontSize',13)
    xlim([0,4.5])


    subplot(4,2,7:8)
    imagesc(z_bz0, r_bz0, bz0');       % NO transpose
    set(gca,'YDir','normal');
end








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
% if do_vis == 1
%     figure; imagesc(z, r, Bz);
%     set(gca,'YDir','normal','FontName','Times','FontSize',18);
%     xlabel('$z$ [m]','Interpreter','latex');
%     ylabel('$r$ [m]','Interpreter','latex');
%     title('Input $B_z$ [T]','Interpreter','latex');
%     colorbar; axis image tight; xlim([0.5 4.2]);
% end

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

% if do_vis == 1
%     figure; imagesc(z, r, psiN);
%     set(gca,'YDir','normal','FontName','Times','FontSize',18);
%     xlabel('$z$ [m]','Interpreter','latex');
%     ylabel('$r$ [m]','Interpreter','latex');
%     title('Normalized $\psi_N$','Interpreter','latex');
%     colorbar;
% end

% -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES FROM psiN
% -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% CONSTRUCT Te AND ne PROFILES FROM psiN
% -------------------------------------------------------------------------
disp('Constructing Te and ne profiles...');
qe  = 1.602176634e-19;      % [C]
% te_min = 8;   % eV
% te_max = 2;   % eV
% Te = (psiN < 1) .* (te_max - te_min) .* (1 - psiN).^1.75 + te_min;
% Te(psiN > 1) = te_max;   % flatten beyond LCFS

% --------- Te: max in core, lower at LCFS, then exponential falloff ----------
Te0       = parms_tune(4);      % eV at core (psiN=0)  <-- MAX.
Te_edge   = parms_tune(3);      % eV at LCFS (psiN=1) 
Te_floor  = 0.2;    % eV far outside
alpha     = 1.75;   % core shaping
lambdaPsi = 0.15;   % SOL falloff width (smaller = steeper)

Te = zeros(size(psiN));

% Core (psiN <= 1): Te(0)=Te0, Te(1)=Te_edge
Te(psiN <= 1) = Te_edge + (Te0 - Te_edge) .* (1 - psiN(psiN <= 1)).^alpha;

% SOL (psiN > 1): exponential starting from Te_edge at psiN=1
Te(psiN > 1)  = Te_floor + (Te_edge - Te_floor) .* exp(-(psiN(psiN > 1) - 1) ./ lambdaPsi);


% Density model (helicon-only): fixed ne_max target + required-power check
% Set these two knobs per case:
gasPuff_Gps = 1.0e20;   % gas puff/source rate [particles/s]
PRF_kW      = 180.0;     % coupled RF power [kW]

Eion_eV     = 1.0e3;    % ionization cost [eV/pair]
ne_cap      = parms_tune(1);   % target max density [m^-3] % 

% RF power required to sustain source rate G at the ionization cost
Preq_kW = gasPuff_Gps * Eion_eV * qe / 1e3;
powerFactor = min(PRF_kW / max(Preq_kW,1e-12), 1.0);  % saturates at 1

% User-fixed density target
ne_min = parms_tune(2);  % m^-3     
ne_max = ne_cap;  % user-fixed target


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
% Mach is linear in z: M=-1 at left boundary, +1 at right boundary
% -------------------------------------------------------------------------
disp('Constructing parallel velocity profile from Te...');

mp  = 1.67262192369e-27;    % [kg] proton mass (assume H+)

Cs = sqrt(qe .* Te ./ mp);  % [nR × nZ] ion sound speed

% Piecewise Mach profile in z (to match requested stepped/ramped shape)
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
bvec = Bmag(:); bvec = bvec(isfinite(bvec) & bvec>0);
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

% === CRITICAL CHANGE: compute grad(psi) from B-field (axisymmetric identity)
rMat = repmat(r(:), 1, nZ);

dpsidr = 2*pi .* rMat .* Bz;      % ∂ψ/∂r
dpsidz = -2*pi .* rMat .* Br;     % ∂ψ/∂z

gradPsiMag = sqrt(dpsidr.^2 + dpsidz.^2);

% Prevent blow-up where gradPsi is tiny
gvec = gradPsiMag(:); gvec = gvec(isfinite(gvec) & gvec>0);
if isempty(gvec)
    eps_floor = 1e-12;
else
    eps_floor = prctile(gvec, 2);
    eps_floor = max(eps_floor, 1e-12);
end
gradPsiMag_safe = max(gradPsiMag, eps_floor);

nr = dpsidr ./ gradPsiMag_safe;
nz = dpsidz ./ gradPsiMag_safe;

% Mask truly tiny gradients (optional)
mask_bad = gradPsiMag < eps_floor;
nr(mask_bad) = 0; nz(mask_bad) = 0;

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
% HEAT-FLUX PROFILES
%   1) Lineout vs z at r = 0.06 m
%   2) Radial profile averaged over 1.6 < z < 1.9
% -------------------------------------------------------------------------
[~, iR06] = min(abs(r - 0.06));
zMask = (z > 1.6) & (z < 1.9);

[~, iZ45] = min(abs(z - 4.1926));
rMask = (r < 0.06);


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


Qtot_map = Qpar + Qperp;

if do_vis == 1
    subplot(4,2,2)
    imagesc(z, r, Qtot_map);
    set(gca,'YDir','normal','FontName','Times','FontSize',18);
end

r_h = r(iR06);
z_h = z(zMask);
Q_tot_h = abs(Qtot_map(iR06,zMask));
r_t = r(rMask);
z_t = z(iZ45);
Q_tot_t = Qtot_map(rMask, iZ45);

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



% if do_vis == 1
%     figure;
%     hPar1 = semilogy(z, qpar_z_r06, 'k-', 'LineWidth', 2.0); hold on;
%     hPerp1 = semilogy(z, qperp_z_r06, 'r-', 'LineWidth', 2.0);
%     hTot1  = semilogy(z, qtot_z_r06, 'b-', 'LineWidth', 2.0);
%     hPsi1  = semilogy(z, qpsi_z_r06, 'm--', 'LineWidth', 1.8);
%     xline(1.6, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
%     xline(1.9, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
%     set(gca,'FontName','Times','FontSize',18);
%     xlabel('$z$ [m]','Interpreter','latex');
%     ylabel('Heat flux [W m$^{-2}$]','Interpreter','latex');
%     ylim([1e0 1e7]);
%     title(sprintf('Absolute heat-flux lineout at $r=%.4f$ m', r(iR06)), 'Interpreter','latex');
%     leg1 = {'$|q_\parallel|$','$|q_\perp|$','$|q_{\mathrm{tot}}|$','$|q_{\mathrm{tot},\psi}|$'};
%     legend([hPar1 hPerp1 hTot1 hPsi1], leg1, 'Interpreter','latex','Location','best');
%     grid on; box on;
% 
%     figure;
%     hPar2 = semilogy(r, qpar_r_win, 'k-', 'LineWidth', 2.0); hold on;
%     hPerp2 = semilogy(r, qperp_r_win, 'r-', 'LineWidth', 2.0);
%     hTot2  = semilogy(r, qtot_r_win, 'b-', 'LineWidth', 2.0);
%     hPsi2  = semilogy(r, qpsi_r_win, 'm--', 'LineWidth', 1.8);
%     xline(0.06, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
%     set(gca,'FontName','Times','FontSize',18);
%     xlabel('$r$ [m]','Interpreter','latex');
%     ylabel('Heat flux [W m$^{-2}$]','Interpreter','latex');
%     ylim([1e0 1e7]);
%     title('Absolute radial heat-flux profile averaged over $1.6<z<1.9$', 'Interpreter','latex');
%     leg2 = {'$|q_\parallel|$','$|q_\perp|$','$|q_{\mathrm{tot}}|$','$|q_{\mathrm{tot},\psi}|$'};
%     legend([hPar2 hPerp2 hTot2 hPsi2], leg2, 'Interpreter','latex','Location','best');
%     grid on; box on;
% end





