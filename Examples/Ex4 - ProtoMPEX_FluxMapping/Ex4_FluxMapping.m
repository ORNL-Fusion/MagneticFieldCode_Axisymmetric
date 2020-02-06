% OBJECTIVE:
% 1- Demonstrate the use of the functions in this repo to calculate the
% magnetic field along ProtoMPEX and its associated flux mapping
% 2- Demonstrate the construction of the vacuum boundary

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
protoMpexConf = 'conf_G';

% Assignment of currents per power supply:
% =========================================================================
coilCurrents.TR1 = 530;
coilCurrents.TR2 = 2300;
coilCurrents.PS1 = 5000;
coilCurrents.PS2 = 4000;
coilCurrents.PS3 = 430;

% =========================================================================
% Read "CoilSetup" spreadsheet:
dum1 = tic;
disp('Reading "coilSetup" spreadsheet...')
coilSetup = readtable('CoilSetup_ProtoMPEX.xlsx','Sheet',protoMpexConf);
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
z_Target = 4.0344;
r1D = linspace(1e-3,0.1  ,60 );
z1D = linspace(z_Dump,4.5,501);

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

%% Draw Proto-MPEX vacuum vessel:
% Inch to meter conversion factor
in2m = 0.0254;

% Dump plate geometry:
rDump = 15.75*in2m/2;
zDump = 0.5;

% Central chamber:
zCC1 = coil{6}.z + 0.5*coil{6}.dz;
zCC2 = coil{7}.z - 0.5*coil{7}.dz;
rCC  = 24*in2m/2;

% Vaccum vessel:
rVac1 = 5.834*in2m/2; % Vacuum vessel rad that is on either side of ALN window
rVac2a = 19.25*in2m/2; % Wide sections dowmstream of central chamber to end
rVac2b = 4.272*2*in2m/2;  % Narrow sections are just inside of coil
zVac1  = zCC1;
zVac2a = zCC2;
zVac2b = coil{end}.z + 0.5*coil{end}.dz;

% Helicon window:
L =  11.8*in2m; % L is 11.8" from Meitner
r1 = 4.95*in2m/2; % diameter is 5.1 inches --> update 4/21/16.  Nominal is ~4.95;
r2 = rVac1;
z1 = 0.5*(coil{3}.z + coil{4}.z) - L/2;
z2 = 0.5*(coil{3}.z + coil{4}.z) + L/2;
heliconWindow.r = [r2,r1,r1,r2];
heliconWindow.z = [z1,z1,z2,z2];

% MPEX-like limiter:
limiterLength = 30e-2;
limiterWidth  = 3e-3;
z1 = heliconWindow.z(end) + 1e-2;
z2 = z1 + limiterLength;
r1 = 2.5*in2m - limiterWidth;
r2 = rVac1;

limiter.r = [r2,r1,r1,r2];
limiter.z = [z1,z1,z2,z2];

% Skimmer:
r1 = (7/100)/2;
r2 = rVac1;
z1 = coil{5}.z + 0.0702 - 0.5e-2;
z2 = coil{5}.z + 0.0702 + 0.5e-2;
skimmer.r = [r2,r1,r1,r2];
skimmer.z = [z1,z1,z2,z2];

% ECH heating region:
r1 = rVac2b;
r2 = rVac2a;
z1 = coil{8}.z + 0.5*coil{8}.dz;
z2 = coil{9}.z - 0.5*coil{9}.dz;
echSection.r = [r1,r2,r2,r1];
echSection.z = [z1,z1,z2,z2];

% ICH Sleeve: -------------------------------------------------------------
r1 = rVac2b;
r2 = 0.08/2; % ID = 80 mm
z1 = coil{09}.z - 0.5*coil{09}.dz;
z2 = coil{10}.z + 0.5*coil{10}.dz;
ichSleeve.r = [r1,r2,r2,r1];
ichSleeve.z = [z1,z1,z2,z2];

% Space between coil 10-11
r1 = rVac2b;
r2 = rVac2a;
z1 = coil{10}.z + 0.5*coil{10}.dz;
z2 = coil{11}.z - 0.5*coil{11}.dz;
space1.r = [r1,r2,r2,r1];
space1.z = [z1,z1,z2,z2];

% Space between coil 11-12
r1 = rVac2b;
r2 = rVac2a;
z1 = coil{11}.z + 0.5*coil{11}.dz;
z2 = coil{12}.z - 0.5*coil{12}.dz;
space2.r = [r1,r2,r2,r1];
space2.z = [z1,z1,z2,z2];

% Space between coil 12-13
r1 = rVac2b;
r2 = rVac2a;
z1 = coil{12}.z + 0.5*coil{12}.dz;
z2 = coil{13}.z - 0.5*coil{13}.dz;
space3.r = [r1,r2,r2,r1];
space3.z = [z1,z1,z2,z2];

vessel.r = [0    ,rDump,rDump    ,rVac1     ,rVac1,rCC  ,rCC   ,rVac2b,rVac2b];
vessel.z = [zDump,zDump,zDump+0.2,zDump+0.2 ,zVac1,zVac1,zVac2a,zVac2a,zVac2b];

vessel = AddComponent(vessel,echSection);
vessel = AddComponent(vessel,ichSleeve);
vessel = AddComponent(vessel,space1);
vessel = AddComponent(vessel,space2);
vessel = AddComponent(vessel,space3);
vessel = AddComponent(vessel,heliconWindow);
vessel = AddComponent(vessel,skimmer);
vessel = AddComponent(vessel,limiter);

vessel = SegmentBoundary(vessel,0.02);

figure;
plot(vessel.z,vessel.r,'k.-')
xlim([0,5])
ylim([0,1])

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
        Phi0 = Phi2D(nz,nr);
    case 2 % Based on minimum flux at boundary
        rng_ii = find(vessel.z > 0.5 & vessel.z < 3.7);
        Phi_min = ones(size(vessel.z));
        for ii = rng_ii
            [~,~,~,Phi_min(ii),~,~] = CalculateMagField(coil,vessel.z(ii),vessel.r(ii));
        end
        [~,ii] = min(Phi_min);
        rlimit = vessel.r(ii);
        zlimit = vessel.z(ii);
        nr = find(r1D > rlimit,1);
        nz = find(z1D > zlimit,1);
        Phi0 = Phi2D(nz,nr);
end

% Flux coordinate
xi = Phi2D/Phi0;

%% SECTION 4: MAGNETIC FIELD LINES AND PLASMA EDGE
% =========================================================================
% Magnetic field field lines up to the plasma edge
xi_lines = linspace(0.01,1,20);
for ii = 1:numel(xi_lines)
    C = contour(z2D,r2D,xi,[1,1]*xi_lines(ii));
    z_fluxline{ii} = C(1,2:end);
    r_fluxline{ii} = C(2,2:end);
end
clearvars ii

% Truncate Xi when > 1
% dum1 = find(xi>1);
% xi(dum1) = [-1];

figure
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

plot(vessel.z,vessel.r,'k-')
plot(zlimit,rlimit,'ro')

xlim([0,5  ])
ylim([0,0.15])

return

%% SECTION 5: Save figure
% =========================================================================
% Saving figure:

InputStructure.prompt = {['Would you like to save figure? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
saveFig = GetUserInput(InputStructure);

if saveFig
    saveas(gcf,'Validating Magnetic field code','tiffn')
end

% =========================================================================
disp('End of script')