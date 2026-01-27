% Magnetic field in Proto-MPEX

% Written by J.F. Caneses Marin
% Created on 2021-10-11

%% SECTION 1: Read "CoilSetup" spreadsheet
clearvars
clc
close all

saveFig = 1;
saveData = 1;

% #########################################################################
%                       INPUT FROM USER:
% #########################################################################

% Assignment of currents per power supply:
% =========================================================================

switch 3
    case 1
        numData=6;
        % ii = 1;% Only for Josh Proto-MPEX Expt. # 28950
        % coilCurrents{ii}.TR1 = 542;
        % coilCurrents{ii}.TR2 = 6597;
        % coilCurrents{ii}.PS1 = 3469;
        % coilCurrents{ii}.PS2 = 6017;
        % coilCurrents{ii}.PS3 = 138;
        
        
        %B Scan
        ii = 1;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3720;
        coilCurrents{ii}.PS1 = 2000;
        coilCurrents{ii}.PS2 = 4000;
        coilCurrents{ii}.PS3 = 200;
        
        
        
        
        ii = 2;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3720;
        coilCurrents{ii}.PS1 = 3000;
        coilCurrents{ii}.PS2 = 4000;
        coilCurrents{ii}.PS3 = 200;
        
        ii = 3;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3720;
        coilCurrents{ii}.PS1 = 4000;
        coilCurrents{ii}.PS2 = 4000;
        coilCurrents{ii}.PS3 = 200;
        
        ii = 4;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3720;
        coilCurrents{ii}.PS1 = 5000;
        coilCurrents{ii}.PS2 = 4000;
        coilCurrents{ii}.PS3 = 200;
        
        ii = 5;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3720;
        coilCurrents{ii}.PS1 = 6000;
        coilCurrents{ii}.PS2 = 4000;
        coilCurrents{ii}.PS3 = 200;
        
        ii = 6;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3720;
        coilCurrents{ii}.PS1 = 7000;
        coilCurrents{ii}.PS2 = 4000;
        coilCurrents{ii}.PS3 = 200;


    case 2
% ICH Proto-MPEX Expt. # 24775
        numData=10;

        % MPEX-like limiter: case 6 ( EC+Impurity case 05/06/2020, Shot #: 30000)
        % coilCurrents{1}.TR1 = 530;
        % coilCurrents{1}.TR2 = 2100;
        % coilCurrents{1}.PS1 = 6800;
        % coilCurrents{1}.PS2 = 3500;
        % coilCurrents{1}.PS3 = 400;

        ii = 1;
        
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 2000;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 2;
        
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 2500;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 3;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 2800;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 4;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 3500;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 5;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 4000;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 6;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 4500;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 7;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 5000;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 8;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 5500;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 9;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 6000;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        ii = 10;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 7000;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        case 3
% ICH Proto-MPEX Expt. test Bfield with no gradient
        numData=1;

   
        
        ii = 1;
        coilCurrents{ii}.TR1 = 520;
        coilCurrents{ii}.TR2 = 3500;
        coilCurrents{ii}.PS1 = 2800;
        coilCurrents{ii}.PS2 = 4500;
        coilCurrents{ii}.PS3 = 220;
        
        
          case 4
% ICH Proto-MPEX Expt. Josh case 28950
        numData=1;

   
        
        ii = 1;
        coilCurrents{ii}.TR1 = 542;
        coilCurrents{ii}.TR2 = 6597;
        coilCurrents{ii}.PS1 = 3469;
        coilCurrents{ii}.PS2 = 6017;
        coilCurrents{ii}.PS3 = 138;
end




% Magnetic configuration of interest:
confType = 'conf_E';
        
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

% Correction:
coilSetup.datum = 5*ones(size(coilSetup.datum));

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
r1D = linspace(1e-3,1e-3,400);
z1D = linspace(z_Dump,z_Target,500);

% Calculate the magnetic field and magnetic vector potential:
% =========================================================================
dum1 = tic;
disp(['Calculating magnetic field for ',num2str(numel(coilCurrents)),' cases ...'])
for ii = 1:numData
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

% SECTION 3: plot magnetic field
%%=========================================================================

figure('color','w')

yyaxis left
imagesc(z1D, r1D, Bz2D')
set(gca,'YDir','normal')
yyaxis right
hold on
for ii = 1:1
    plot(z1D,B2D{ii})
end

set(gcf,'position',[482   144   560   363])
xlabel('x [m]')
ylabel('B [T]')


% SECTION 4: Interpolate data
% =========================================================================

% Select subset:
xL = 0.5;
xR = 4.2;

% Select size:
N = 200;

% Create query vector:
xq = linspace(xL,xR,N)';

% Assign data:
for ii = 1:numData
    B(:,ii) = interp1(z1D',B2D{ii}(:,1),xq);
    plot(xq,B(:,ii),'lineWidth',2)
end

% B(find(xq>4.04))=0.839531;




%% SECTION 5: Normalize data:
% =========================================================================

% Search B reference:
xmin = 1.5;
xmax = 1.9;

% Select spatial search range:
rng = find(xq >= xmin & xq <= xmax);

% Extract B refernce:
[B0,in] = min(B(rng,1));
x0 = xq(rng(in));

% Normalize data:
B_norm = B/B0;

figure('color','w')
hold on
for ii = 1:numData
    plot(xq,B(:,ii))
end
plot(x0,1,'ko')
plot(x0,1,'k.')
box on

set(gcf,'position',[482   144   560   363])
xlabel('x [m]')
ylabel('B/B_0')
titleText = ['B_0 = ',num2str(B0),' , X_0 = ',num2str(x0),' , X_L = ',num2str(xL),...
    ' , X_R = ',num2str(xR),' , N = ',num2str(N)];
title(titleText) 

% % Save figure:
% % =========================================================================
% if saveFig
%     figureName = 'ProtoMPEX_B_dataset_1';
% 
%     % PDF figure:
%     exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 
% 
%     % TIFF figure:
%     exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
% end

% Save data to text file:
% =========================================================================
if saveData
    fileName = 'ProtoMPEX_B_dataset_1';

    for ii = 1:numel(B2D)
        f1 = [B_norm(:,ii)];
        save([fileName,'_PS1_',num2str(coilCurrents{ii}.PS1),'.txt'],'f1','-ascii');
    end

    fileName = ['x_',fileName];
    f1 = [xq];
    save([fileName,'.txt'],'f1','-ascii');
end

% return;
%% Calculate B-field and write it in .csv and .nc format
disp('Writing the B-field profiles');

nR=length(r1D);
nZ=length(z1D);

r=r1D; z=z1D;
Br=Br2D'; Bz=Bz2D'; Bt=zeros(size(Bz2D))';

ncid = netcdf.create(('./bfield_protoMPEX.nc'),'NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);

dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);

gridZnc = netcdf.defVar(ncid,'z','float',dimZ);

brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,r);
netcdf.putVar(ncid,gridZnc,z);


netcdf.putVar(ncid,brnc,Br);
netcdf.putVar(ncid,btnc,Bt);
netcdf.putVar(ncid,bznc,Bz);

netcdf.close(ncid);

hold off

% Quiver plot for Br, Bz
figure
fileB = 'bfield_protoMPEX.nc';
% file1='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/protoMPEX/parametricScan/no_diffision_flag2/densityScan/te8to8ne1e18to1e19/input/profiles_protoMPEX_SOLPS.nc';
xB = ncread(fileB,'x');
zB = ncread(fileB,'z');

br0 = ncread(fileB,'br');
bt0 = ncread(fileB,'bt');
bz0 = ncread(fileB,'bz');

B_0 = sqrt(br0.^2+bz0.^2);
figure; plot(z,br0(400,:)./B_0)

figure;
imagesc(z,r,Bz)
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24)
xlabel('$z$ [m]','interpreter','latex','fontSize',24)
xlim([0.5 4.2])

disp('Calculated the B-field profile');





disp('End of script')