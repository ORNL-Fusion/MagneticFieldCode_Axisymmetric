% OBJECTIVE:
% 1- Build on Ex4_FluxMapping to allow simultaneous plotting of various
% magnetic field configurations (flux mapping and field profile)

% Written by J.F. Caneses Marin
% Created on 2020-02-06

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

% =========================================================================
% Read "CoilSetup" spreadsheet:
dum1 = tic;
disp('Reading "coilSetup" spreadsheet...')
coilSetup = readtable('CoilSetup_ProtoMPEX.xlsx','Sheet',confType);
disp('Reading complete!')
toc(dum1);

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
z_Dump = 0.5;
z_Target = 4.2;
r1D = linspace(1e-3,0.1  ,40 );
z1D = linspace(z_Dump,z_Target,300);

% Assignment of currents per power supply:
% =========================================================================
% MPEX-like limiter: case 1
coilCurrents{1}.TR1 = 530;
coilCurrents{1}.TR2 = 3000;
coilCurrents{1}.PS1 = 3500;
coilCurrents{1}.PS2 = 3500;
coilCurrents{1}.PS3 = 330;

% MPEX-like limiter: case 2
% coilCurrents{2}.TR1 = 530;
% coilCurrents{2}.TR2 = 2300;
% coilCurrents{2}.PS1 = 3500;
% coilCurrents{2}.PS2 = 4000;
% coilCurrents{2}.PS3 = 430;

% Window limiter: case 3
% coilCurrents{3}.TR1 = 530;
% coilCurrents{3}.TR2 = 2300;
% coilCurrents{3}.PS1 = 6800;
% coilCurrents{3}.PS2 = 4000;
% coilCurrents{3}.PS3 = 160;

% Window limiter: case 4
% coilCurrents{4}.TR1 = 530;
% coilCurrents{4}.TR2 = 2300;
% coilCurrents{4}.PS1 = 6000;
% coilCurrents{4}.PS2 = 4000;
% coilCurrents{4}.PS3 = 220;

% Upstream vaccum vessel limiter: case 5
% coilCurrents{5}.TR1 = 0;
% coilCurrents{5}.TR2 = 2300;
% coilCurrents{5}.PS1 = 5000;
% coilCurrents{5}.PS2 = 4000;
% coilCurrents{5}.PS3 = 350;

% Calculate the magnetic field and magnetic vector potential:
% =========================================================================
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
        z_fluxline{ii}{jj} = C(1,2:end);
        r_fluxline{ii}{jj} = C(2,2:end);
    end
end
close(gcf)
clearvars ii

%% SECTION 5: Plot data

% Magnetic flux mapping:
% =========================================================================
figureName = 'ProtoMPEX_FluxMappingVarious';
figure('color','w','Tag',figureName)
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
        set(gca,'PlotBoxAspectRatio',[2.3 1.5 1])
        xlim([0.25,4.5])
        ylim(0.35*[-1,+1])
end
box on
xlabel('z [m]','Interpreter','Latex','FontSize',13)
ylabel('r [m]','Interpreter','Latex','FontSize',13)
grid on

% Magnetic field profiles:
% =========================================================================
figureName = 'ProtoMPEX_BzProfileVarious';
figure('color','w','Tag',figureName)
hold on
for ii = 1:numel(coilCurrents)
    hBz(ii) = plot(z1D,B2D{ii}(:,1),lineColor{ii},'LineWidth',2);
end
box on 
grid on
set(gca,'PlotBoxAspectRatio',[2 1 1])
set(gca,'FontName','times')
xlabel('z [m]','Interpreter','Latex','FontSize',13)
ylabel('B$_0$ [T]','Interpreter','Latex','FontSize',13)
xlim([0.25,4.5])

% Target:
hT = line(z_Target*[1,1],[0,1]);
set(hT,'color','k','LineWidth',4)

%% SECTION 6: Save figure
% =========================================================================
% Saving figure:

InputStructure.prompt = {['Would you like to save figure? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
saveFig = GetUserInput(InputStructure);

if saveFig
    figureName = 'ProtoMPEX_FluxMappingVarious';
    hdum1 = findobj('Tag',figureName);
    saveas(hdum1,figureName,'tiffn')
    
    figureName = 'ProtoMPEX_BzProfileVarious';
    hdum1 = findobj('Tag',figureName);
    saveas(hdum1,figureName,'tiffn')
end

% =========================================================================
disp('End of script')