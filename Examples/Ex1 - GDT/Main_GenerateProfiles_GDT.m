% Based on a given magnetic coil geometry generate 2D profiles of the
% magnetic field, magnetic flux and plasma density in a linear plasma configuration

% Created 2019_12_09, JF Caneses

% Scripts/function needed to run this code:
% DefineCoilTypes.m
% GenerateCoilSetup.m
% AddComponent.m
% CalculateMagField.m
% CreateCoilStructure.m
% bfield_circular_coil_analytic.m

%% SECTION 1: Create "coil" structure
% =========================================================================
clearvars
clc
close all

% =========================================================================
% Define the search path:
homeFolder = cd;
% Go back two folders
cd ..\.. 
currentFolder = cd;
functionFolder = [currentFolder,'\Functions'];
p = addpath(genpath(functionFolder));
cd(homeFolder)

% =========================================================================
% Create "coilSetup" structure:
% The geometry of the coils can be defines either using an excell
% spreadsheet called "coilSetup_GDT" or from a script called
% "GenerateCoilSetup"

% The main output is to produce a structure called "coilSetup" that
% contains all the geometrical information of all the coils

coilSetupType = 2;
switch coilSetupType
    case 1
        % Read "CoilSetup" spreadsheet:
        dum1 = tic;
        disp('Reading "coilSetup" spreadsheet...')
        coilSetup = readtable('CoilSetup_GDT.xlsx','Sheet','conf_3');
    case 2
        dum1 = tic;
        disp('Running "GenerateCoilSetup" script...')
        GenerateCoilSetup
end

% =========================================================================
% Assignment of currents per power supply:
coilCurrents.expander_1  = +0e3; % Expander coils
coilCurrents.expander_2  = -0.5e3; % Expander coils
coilCurrents.mirror      = +15e3; % Mirror coils
coilCurrents.centralCell = +5e3; % Central cell coils
coilCurrents.transition  = +10e3; % Transition coils

% =========================================================================
% Create "coil" structure: 
% Based the "coilSetup" and "coilCurrents", we produce a structure called
% "coil" for each physical coil. The "coil" structure contains the position
% of each current filament loop, power supply current.

disp('Creating "coil" structure...')
[coil] = CreateCoilStructure(coilSetup,coilCurrents);
disp(['Complete! Elapsed time: ',num2str(toc(dum1)),' s'])
clearvars dum*

if 1
    % =========================================================================
    % Plot Magnetic coils and current filaments:
    figure('color','w','Tag','coilSetup')
    hold on
    for ii = 1:numel(coil)
        plot(coil{ii}.zfil,+coil{ii}.rfil,'r.');
        plot(coil{ii}.zfil,-coil{ii}.rfil,'r.');
    end
    % Formatting:
    set(gca,'FontName','times','FontSize',11)
    xlabel('z [m]','Interpreter','latex','FontSize',13)
    ylabel('r [m]','Interpreter','latex','FontSize',13)
    box on
    xlim([coil{1}.z,coil{end}.z])
    ylim([-2.0     ,2.0        ])
    grid on
end

%% SECTION 2: Calculate magnetic field
% =========================================================================
% Define the area to evaluate the fields at:
tic
r1D = linspace(1e-3,0.5,200);
z1D = linspace(-10 ,10 ,300);

% =========================================================================
% Calculate the magnetic field and magnetic vector potential:
dum1 = tic;
disp('Calculating magnetic field...')
[Br2D,Bz2D,Atheta2D,Phi2D,z2D,r2D] = CalculateMagField(coil,z1D,r1D);
disp(['Complete! Elapsed time: ',num2str(toc(dum1)),' s'])
clearvars dum*

% =========================================================================
% magnetic field magnitude:
B2D = sqrt(Br2D.*Br2D  + Bz2D.*Bz2D);
toc
%% SECTION 4: REFERENCE FLUX
% =========================================================================
% Define reference flux
limitType = 2;
switch limitType
    case 11 % Define throat at the mirror coil
        rlimit = 0.6*coil{2}.r1;
        nr = find(r1D > rlimit,1);
        zlimit = coil{2}.z;
        nz = find(z1D > zlimit ,1);
        Phi0 = Phi2D(nz,nr);
    case 2 % Define size at the midplane
        rlimit = 6/100;
        nr = find(r1D > rlimit,1);
        zlimit = 0;
        nz = find(z1D > zlimit ,1);
        Phi0 = Phi2D(nz,nr);
end

% Flux coordinate
xi = Phi2D/Phi0;

% Plasma edge:
xi_lines = [0.01:0.1:1];
for ii = 1:numel(xi_lines)
    C = contour(z2D,r2D,xi,[1,1]*xi_lines(ii));
    z_fluxline{ii} = C(1,2:end);
    r_fluxline{ii} = C(2,2:end);
end
clearvars ii

% Truncate Xi when > 1
dum1 = find(xi>1);
xi(dum1) = [-1];

%% SECTION 5: DEFINE PLASMA DENSITY and Te 2D PROFILE
% =========================================================================
% Define the radial plasma density and Te profile:
a = 2;
b = 1.75;
% Peak values:
ne_peak = 2e20;
Te_peak = 600;
% Edge values:
ne_edge = 1e17;
Te_edge = 30;

% Radial shape function:
% =========================================================================
f_r = real(( 1 - (xi.^a) ).^b); 

% =========================================================================
% Axial shape function:
axialType = 3;
switch axialType
    case 1
        Lz = 2*(z1D(end)-z1D(1));
        f_z = cos(2*pi*(z2D/Lz));
    case 2
        Lz = 4;
        dum1 = z2D/Lz;
        f_z = exp(-0.5*(dum1).^2);
    case 3 % hyperbolic tan
        Lz = 1.8;
        zmirror1 = coil{1    +1}.z;
        zmirror2 = coil{end  -1}.z;
        f_z  = 0.5*(tanh((z2D-zmirror1)/Lz) - tanh((z2D+zmirror1)/Lz));  
end
clearvars dum1
% =========================================================================
% 2D plasma profiles:
ne2D = ne_peak*f_r.*f_z  + ne_edge;
Te2D = Te_peak*f_r.*f_z  + Te_edge;

if 1
    figure('color','w','Tag','ne profile')
    subplot(2,1,1)
    hold on
    % Plasma profile at the midplate
    [~,ndum] = find(z1D>0,1);
    dum1(1) = plot(r1D,f_r(ndum,:),'k-','LineWidth',2);
    plot(-r1D,f_r(ndum,:),'k-','LineWidth',2);
    % throat profile
    [~,ndum] = find(z1D>coil{2}.z,1);
    dum1(2) = plot(r1D,f_r(ndum,:),'r-','LineWidth',2);
    plot(-r1D,f_r(ndum,:),'r-','LineWidth',2);
    set(gca,'FontName','times','FontSize',11)
    xlim([-0.15,0.15])
    box on 
    grid on
    ylabel('$f_r(r)$','Interpreter','latex','FontSize',14)
    xlabel('r [m]','Interpreter','latex','FontSize',14)
    legend(dum1,'Central cell','Mirror throat')
    
    subplot(2,1,2)
    hold on
    plot(z1D,f_z(:,1),'k-','LineWidth',2);
    % Magnetic coils:
    for ii = 1:numel(coil)
        plot(coil{ii}.zfil,coil{ii}.rfil,'r.');
    end
    set(gca,'FontName','times','FontSize',11)
    box on 
    grid on
    ylabel('$f_z(z)$','Interpreter','latex','FontSize',14)
    xlabel('z [m]','Interpreter','latex','FontSize',14)
end

%% SECTION 6: Generate Vacuum boundary
% =========================================================================
% Define boundary:
vessel.z = [0  ,5.5,6.5,7.0,8.5,10.0,10.0];
vessel.r = [0.3,0.3,0.1,0.1,0.6,0.6 ,0.0 ];

% =========================================================================
% Mirror vessel and add:
vessel_mirror.z = -fliplr(vessel.z);
vessel_mirror.r = +fliplr(vessel.r);
vessel = AddComponent(vessel,vessel_mirror);

%% SECTION 7: Cyclotron resonance contours
% RF frequency:
f_RF = [11]*1e6;
% Ion mass:
m_p = 1.6726e-27;
m_i = 2*m_p;
% fundamental cyclotron resonance: 
e_c = 1.6020e-19;
B_res = 2*pi*f_RF*m_i/e_c;

% Define contour lines:
dum1 = contour(z2D,r2D,B2D,[1,1]*B_res);
z_resLayer = dum1(1,2:end);
r_resLayer = dum1(2,2:end);

%% SECTION 8: Plot data
% close all
% =========================================================================
% Plot magnetic field profile:
figure('color','w','Tag','GDT B profile')
hold on
% Bz profile:
n_r = find(r1D<0.01,1);
B_onAxis = B2D(:,n_r);
hdum(1) = plot((z1D),B_onAxis,'k','LineWidth',3);
% Magnetic coils:
for ii = 1:numel(coil)
    plot(coil{ii}.zfil,coil{ii}.rfil,'r.');
end
% Formatting:
set(gca,'FontName','times','FontSize',11)
ylabel('B [T]','Interpreter','latex','FontSize',13)
grid on
box on
xlim([z1D(1),z1D(end)])
hL = legend(hdum,'$B_0$');
set(hL,'Interpreter','latex','FontSize',13)
clearvars hdum

% =========================================================================
% Plot flux lines and warm plasma density profile:
figure('color','w','Tag','GDT flux tubes')
hold on
% ne contours:
contourf(z2D,+r2D,ne2D,50,'LineStyle','none')
contourf(z2D,-r2D,ne2D,50,'LineStyle','none')
colormap(flipud(hot))
colorbar('TickLabelInterpreter','latex','FontSize',12)
% Magnetic coils:
for ii = 1:numel(coil)
    plot(coil{ii}.zfil,+coil{ii}.rfil,'r.');
    plot(coil{ii}.zfil,-coil{ii}.rfil,'r.');
end
% Flux lines:
for ii = 1:numel(xi_lines)
    dum1 = plot((z_fluxline{ii}),+r_fluxline{ii},'k:');
    dum2 = plot((z_fluxline{ii}),-r_fluxline{ii},'k:');
    if ii == numel(xi_lines)
        set(dum1,'LineStyle','-','LineWidth',1)
        set(dum2,'LineStyle','-','LineWidth',1)
    end
end
% Formatting:
set(gca,'FontName','times','FontSize',11)
xlabel('z [m]','Interpreter','latex','FontSize',13)
ylabel('r [m]','Interpreter','latex','FontSize',13)
box on
xlim([z1D(1),z1D(end)])
ylim([-1.0  ,     1.0])
% Add vacuum boundary:
plot(vessel.z,+vessel.r,'k','LineWidth',3);
plot(vessel.z,-vessel.r,'k','LineWidth',3);

% 1 Tesla resonance
plot(z_resLayer,+r_resLayer,'r.')
plot(z_resLayer,-r_resLayer,'r.')
if 0
    set(gca,'PlotBoxAspectRatio',[1.5 1 1])
    grid on
    axis('image')
end

%% SECTION 9: Save figure
% =========================================================================
% Saving figure:
InputStructure.prompt = {['Would you like to save Figure? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
saveFig = GetUserInput(InputStructure);

if saveFig
    disp('Saving figures...')
    hdum = findobj('Tag','GDT B profile');
    saveas(hdum,hdum.Tag,'tiffn')
    
    hdum = findobj('Tag','GDT flux tubes');
    saveas(hdum,hdum.Tag,'tiffn')
    
   hdum = findobj('Tag','ne profile');
   saveas(hdum,hdum.Tag,'tiffn')
   disp('figures saved succesfully!')
end

%% SECTION 10: Save data
% =========================================================================
% Saving data:
% r2D: 
% z2D:
% Br2D: radial component of magnetic field in 2D
% Bz2D: axial component of magnetic field in 2D
% Atheta2D: azimuthal component of the magnetic vector potential in 2D
% Phi2D: Magnetic flux in 2D calculated as 2*pi*Atheta2D.*r2D
% zlimit: axial location where the plasma edge is defined
% rlimit: radial location where plasma edge is defined
% xi: magnetic flux coordinate [0,1] based on zlimit, rlimit and Phi2D
% f_r: Radial density profile defined as ( 1 - (xi.^a) ).^b
% f_z: axial density profile
% ne2D: plasma density 2D profile defined as ne2D = ne_peak*f_r.*f_z  + ne_edge;
% Te2D: electron temperature 2D profile defined as Te2D = Te_peak*f_r.*f_z  + Te_edge;
% coil: structure for the magnetic coil geomtry and currents
% vessel: structure that contains coordinates that define the vacuum
% boundary.

InputStructure.prompt = {['Would you like to save data? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
svdt = GetUserInput(InputStructure);

if svdt
    disp('Saving data...')
    variableNames = {'r2D','z2D','Br2D','Bz2D','Atheta2D','Phi2D',...
        'zlimit','rlimit','xi','f_r','f_z','ne_peak','ne_edge','ne2D',...
        'Te_peak','Te_edge','Te2D','coil','vessel'};
    fileName = ['GDT_profileData_',date,'.mat'];
    save(fileName,variableNames{:});
    disp('Data saved succesfully!')
end

% =========================================================================
disp('End of script')