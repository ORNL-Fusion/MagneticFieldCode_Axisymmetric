% Magnetic field flux mapping in Proto-MPEX

% Written by J.F. Caneses Marin
% Created on 2020-06-28

%% SECTION 1: Read "CoilSetup" spreadsheet
% clearvars
% clc
% close all

% #########################################################################
%                       INPUT FROM USER:
% #########################################################################

% Assignment of currents per power supply:
% =========================================================================

shotType = 1;

switch shotType
    case 1     
    ii = 1;
    coilCurrents{ii}.TR1 = 530;
    coilCurrents{ii}.TR2 = 3000;
    coilCurrents{ii}.PS1 = 3500;
    coilCurrents{ii}.PS2 = 3500;
    coilCurrents{ii}.PS3 = 650;
    
    % Vacuum vessel details:
    windowType = 'SiN';
    limiterType = 'AlN-limiter';
    skimmerType = '7.0 cm';
    
    % Magnetic configuration of interest:
    confType = 'conf_G';
    
    case 2
   
    case 3 
    
    case 4
   
    case 5

end
        
% #########################################################################

% #########################################################################
 
% =========================================================================


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
r1D = linspace(1e-3,0.13,40 );
z1D = linspace(z_Dump-1,z_Target+1,300);

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

DrawVacuumVessel_Conf_G_VariousWindows

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


%% SECTION 6: CYCLOTRON RESONANCE REGIONS
% Flag to draw resonance layers:
drawResLayer = 0;

% RF frequency:
f_RF = [28]*1e9;

% Particle mass:
m_e = 9.1094e-31;

% Harmonics:
n_harmonic = [1,2,3,4,5];

% nth harmonic cyclotron resonance: 
e_c = 1.6020e-19;
B_res = (2*pi*f_RF*m_e/e_c)./n_harmonic;

% Define contour lines:
for ii = 1:numel(coilCurrents)
    for nn = 1:numel(n_harmonic)
        dum1 = contour(z2D,r2D,B2D{ii},[1,1]*B_res(nn));
        z_resLayer{nn}{ii} = dum1(1,2:end);
        r_resLayer{nn}{ii} = dum1(2,2:end);
    end
end
close all

%% SECTION 7: PLOT DATA

% Magnetic flux mapping:
% =========================================================================
figureName = 'Step_0_FluxMapping';
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

% Cyclotron resonances:
color = {'r','g','m','c','k'};
if drawResLayer
    for ii = 1:numel(coilCurrents)
        for nn = 1:numel(n_harmonic)
            plot(z_resLayer{nn}{ii},+r_resLayer{nn}{ii},color{nn},'Marker','.','LineStyle','none')
            plot(z_resLayer{nn}{ii},-r_resLayer{nn}{ii},color{nn},'Marker','.','LineStyle','none')
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
    case 4
        axis image
        xlim([0.9,2.3])
        ylim(0.25*[-1,+1])
end
box on
xlabel('z [m]','Interpreter','Latex','FontSize',13)
ylabel('r [m]','Interpreter','Latex','FontSize',13)
grid on

% Magnetic field profiles:
% =========================================================================
figureName = 'Step_0_MagneticFieldProfile';
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
xlim([0,5])

% Cyclotron resonance layer:
if drawResLayer
    for nn = 1:numel(n_harmonic)
        hR = line([0,5],[1,1]*B_res(nn));
        set(hR,'color',color{nn},'LineWidth',1)
    end
end

% Target:
hT = line(z_Target*[1,1],[0,0.3]);
set(hT,'color','k','LineWidth',4)

% Dump:
hT = line(z_Dump*[1,1],[0,0.3]);
set(hT,'color','k','LineWidth',4)

%% SECTION 8: SAVE PICTURE
% =========================================================================
% Saving figure:

InputStructure.prompt = {['Would you like to save figure? Yes [1], No [0]']};
InputStructure.option.WindowStyle = 'normal';
saveFig = GetUserInput(InputStructure);

if saveFig
    figureName = 'Step_0_FluxMapping';
    hdum1 = findobj('Tag',figureName);
    saveas(hdum1,figureName,'tiffn')
    
    figureName = 'Step_0_MagneticFieldProfile';
    hdum1 = findobj('Tag',figureName);
    saveas(hdum1,figureName,'tiffn')
end

% =========================================================================
disp('End of script')