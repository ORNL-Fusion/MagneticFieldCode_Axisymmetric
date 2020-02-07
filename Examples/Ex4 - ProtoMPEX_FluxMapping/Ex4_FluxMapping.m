% OBJECTIVE:
% 1- Demonstrate the use of the functions in this repo to calculate the
% magnetic field along ProtoMPEX and its associated flux mapping
% 2- Demonstrate the construction of the vacuum boundary

% Written by J.F. Caneses Marin
% Created on 2020-02-06

% evalContour Branch

%% SECTION 1: Read "CoilSetup" spreadsheet
clearvars
clc
close all

% =========================================================================
% Define the search path:
homeFolder = cd;
cd ..\..
currentFolder = cd;
functionFolder = [currentFolder,'\Functions'];
p = addpath(genpath(functionFolder));
cd(homeFolder)
 
% =========================================================================
% Magnetic configuration of interest:
confType = 'conf_G';

% Assignment of currents per power supply:
% =========================================================================
limitType = 3;
switch limitType
    case 1 % MPEX-like limiter
        coilCurrents.TR1 = 530;
        coilCurrents.TR2 = 2300;
        coilCurrents.PS1 = 5000;
        coilCurrents.PS2 = 4000;
        coilCurrents.PS3 = 430;
    case 2 % Window limiter
        coilCurrents.TR1 = 530;
        coilCurrents.TR2 = 2300;
        coilCurrents.PS1 = 6500;
        coilCurrents.PS2 = 4000;
        coilCurrents.PS3 = 220;
    case 3 % Upstream vaccum vessel limiter
        coilCurrents.TR1 = 0;
        coilCurrents.TR2 = 2300;
        coilCurrents.PS1 = 5000;
        coilCurrents.PS2 = 4000;
        coilCurrents.PS3 = 350;
end

% =========================================================================
% Read "CoilSetup" spreadsheet:
dum1 = tic;
disp('Reading "coilSetup" spreadsheet...')
coilSetup = readtable('CoilSetup_ProtoMPEX.xlsx','Sheet',confType);
disp('Reading complete!')

% =========================================================================
% Create "coil" structure"
disp('Creating "coil" structure...')
[coil] = CreateCoilStructure(coilSetup,coilCurrents);
disp(['Complete! Elapsed time: ',num2str(toc(dum1)),' s'])
clearvars dum*

% =========================================================================
% Display "coilSetup" on CLI:
coilSetup

% #########################################################################
% IMPORTANT QUANTITIES FROM THIS SECTION:
% #########################################################################
% - coilSetup: structure that contains all the geometrical info of the
% coils
% - coilCurrents: strucuture that contains the power supply currents
% - coil: object that describes the specific coil setup and is the input
% for the magnetic field calculator function

%% SECTION 2: Calculate magnetic field
% =========================================================================
% Define the area to evaluate the fields at:
tic
z_Dump = 0.5;
z_Target = 4.2;
r1D = linspace(1e-3,0.1  ,30 );
z1D = linspace(z_Dump,z_Target,300);

% =========================================================================
% Calculate the magnetic field and magnetic vector potential:
dum1 = tic;
disp('Calculating magnetic field...')
[Br2D,Bz2D,~,phi2D,z2D,r2D] = CalculateMagField(coil,z1D,r1D,'grid');
disp(['Complete! Elapsed time: ',num2str(toc(dum1)),' s'])
clearvars dum*

% =========================================================================
% magnetic field magnitude:
B2D = sqrt(Br2D.*Br2D  + Bz2D.*Bz2D);
toc

%% Draw Proto-MPEX vacuum vessel:
if ~strcmpi(confType,'conf_G')
    error('Change "confType" to config_G')
end
DrawVacuumVessel_Conf_G

%% SECTION 3: Calculate reference flux:
% =========================================================================
% Define reference flux
limitType = 2;
switch limitType
    case 1 % Based on helicon window
        rlimit = min(heliconWindow.r);
        nr = find(r1D > rlimit,1);     
        rng_z = find(z1D>1 & z1D<2.5);
        [~,a] = min(Bz2D(rng_z,nr));
        nz = rng_z(a);
        zlimit = z1D(nz);
        phi0 = phi2D(nz,nr);
    case 2 % Based on phi at vacuum vessel boundary
        % Select region of interest:
        rng_ii = find(vessel_wLim.z > 0.5 & vessel_wLim.z < 3.7);
        % Initialize variable:
        phiBoundary = ones(size(vessel_wLim.z));
        % Interpolate phi along vaccum vessel contour
        zq = vessel_wLim.z(rng_ii);
        rq = vessel_wLim.r(rng_ii);
        a = interp2(z1D,r1D,phi2D',zq,rq);
        phiBoundary(rng_ii) = a;
        % find location of minimim phi along contour
        [~,ii] = min(phiBoundary);
        % Limit physical location:
        rlimit = vessel_wLim.r(ii);
        zlimit = vessel_wLim.z(ii);
        % Extract reference magnetic flux at the limiting location:
        nr = find(r1D > rlimit,1,'first');
        nz = find(z1D > zlimit,1);
        phi0 = interp2(z1D,r1D,phi2D',zlimit,rlimit);
end

% Flux coordinate
xi = phi2D/phi0;

%% SECTION 4: MAGNETIC FIELD LINES AND PLASMA EDGE
% =========================================================================
% Magnetic field field lines up to the plasma edge
xi_lines = linspace(1e-2,1,10);
for ii = 1:numel(xi_lines)
    C = contour(z2D,r2D,xi,[1,1]*xi_lines(ii));
    z_fluxline{ii} = C(1,2:end);
    r_fluxline{ii} = C(2,2:end);
end
close(gcf)
clearvars ii

%% SECTION 5: Plot data
figure('color','w')
hold on
% Magnetic coils:
for ii = 1:numel(coil)
    plot(coil{ii}.zfil,+coil{ii}.rfil,'r.');
    plot(coil{ii}.zfil,-coil{ii}.rfil,'r.');
end
% Flux lines:
for ii = 1:1:numel(xi_lines)
    dum1 = plot((z_fluxline{ii}),+r_fluxline{ii},'k:');
    dum2 = plot((z_fluxline{ii}),-r_fluxline{ii},'k:');
    if ii == numel(xi_lines)
        set(dum1,'LineStyle','-','LineWidth',1)
        set(dum2,'LineStyle','-','LineWidth',1)
    end
end

% Vacuuum vessel
plot(vessel_wLim.z,+vessel_wLim.r,'r','LineWidth',2)
plot(vessel_wLim.z,-vessel_wLim.r,'r','LineWidth',2)
plot(vessel.z,+vessel.r,'k-','LineWidth',2)
plot(vessel.z,-vessel.r,'k-','LineWidth',2)

% Target
plot(zlimit,rlimit,'ro')
hT = line(z_Target*[1,1],0.08*[-1,+1]);
set(hT,'color','k','LineWidth',3)

% Formatting
% set(gca,'PlotBoxAspectRatio',[1 1 1])
set(gca,'FontName','times')
xlim([0.25,4.5])
ylim(0.35*[-1,+1])
box on
xlabel('z [m]','Interpreter','Latex','FontSize',13)
ylabel('r [m]','Interpreter','Latex','FontSize',13)

%% SECTION 6: Save figure
% =========================================================================
% Saving figure:

InputStructure.prompt = {['Would you like to save figure? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
saveFig = GetUserInput(InputStructure);

if saveFig
    figureName = 'ProtoMPEX_FluxMapping';
    saveas(gcf,figureName,'tiffn')
end

% =========================================================================
disp('End of script')