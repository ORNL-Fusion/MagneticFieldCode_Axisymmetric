% OBJECTIVE:
% 1- Demonstrate the use of the functions: "CreateCoilStructure" and
% "CalculateMagneticField" to calculate the magnetic field on Proto-MPEX
% along the 
% 2- Demonstrate the use of the "CoilSetup" spreadsheet
% 3- Verify the numerical implementation against calculations from other
% codes.

% COIL GEOMETRY INFORMATION:
% The coil arrangement and geometrical information is recorDed in a
% spreadsheet called "CoilSetup_ProtoMPEX".

% VERIFICATION EXERCISE:
% The magnetic field in Proto-MPEX is calculated and compared with
% previosly calculated magnetic field based on Jeremy Lore's code.

% Written by J.F. Caneses Marin
% Created on 2019-12-19

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
% Read "CoilSetup" spreadsheet:
dum1 = tic;
disp('Reading "coilSetup" spreadsheet...')
coilSetup = readtable('CoilSetup_ProtoMPEX.xlsx','Sheet','conf_A');
disp('Reading complete!')

% =========================================================================
% Small correction to the coil geometry:
coilSetup.r_inner =  coilSetup.r_inner - 0.5e-2;
 
% =========================================================================
% Assignment of currents per power supply:
coilCurrents.TR1 = 600;
coilCurrents.TR2 = 160;
coilCurrents.PS1 = 4500;
coilCurrents.PS2 = 4500;

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
r1D = linspace(1e-3,0.3 ,10 );
z1D = linspace(0   ,+5.0,501);

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

%% SECTION 3: Get data from Jeremy's code
% =========================================================================
% Get precalculated magnetic field from Jeremy's code:
fileName = 'PS1_4500_PS2_4500_TR1_160_TR2_600.txt  ';
d = load(fileName);
z_Bjl = d(:,1);
Bjl   = d(:,2);

%% SECTION 4: verification, plot data
close all
% =========================================================================
% Compare present calculation and Jeremy's code:
figure('color','w')
subplot(7,1,[1:5])
hold on
n_r = find(r1D<0.01,1);
B_onAxis = B2D(:,n_r);
hdum(1) = plot((z1D),B_onAxis,'k','LineWidth',3);
hdum(2) = plot(z_Bjl,Bjl,'g.','LineWidth',1);
for ii = 1:numel(coil)
    plot(coil{ii}.zfil,coil{ii}.rfil,'r.');
end
set(gca,'FontName','times','FontSize',11)
ylabel('B [T]','Interpreter','latex','FontSize',13)
grid on
box on
xlim([z1D(1),z1D(end)])
hL = legend(hdum,'Present calc','JL');
set(hL,'Interpreter','latex','FontSize',11)
set(gca,'XTickLabel',[])
clearvars hdum
% -------------------------------------------------------------------------
% Plotting error relative to max(B):
subplot(7,1,[6,7])
errorBfield = abs(B_onAxis-Bjl)./max(B_onAxis);
hdum(1) = plot(z1D,100*errorBfield,'k','LineWidth',3);
set(gca,'FontName','times','FontSize',11)
ylabel('|Error| %','Interpreter','latex','FontSize',12)
xlabel('z [m]','Interpreter','latex','FontSize',13)
grid on
ylim([0,2])
clearvars hdum

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