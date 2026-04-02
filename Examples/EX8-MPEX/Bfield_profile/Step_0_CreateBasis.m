% Step_0:
% For a given MPEX coil positioning and geometry, the total magnetic field
% is expressed as a linear combination of so-called "magnetic field basis".
% Given a set of currents for each coil, the total magnetic field can be
% calculated using the "magnetic field basis" as unit vectors and the
% currents are coefficients
% The magnetic field basis is produced by computing the magnetic field vector
% produced by each coil for a unit current over the entire computational
% domain of interest.

 clear all
 close all
 clc
  
 compute = 0;

if compute
    t2 = tic;
    disp('Calculating "magnetic field basis...')
    
    % Define computational grid:
    % =========================================================================
    LZ_max = +10;
    LZ_min = -2;
    NZ     = 200;
    LR_max = 1.5;
    LR_min = 1e-3;
    NR     = 100;

    % Add path to repo with the magnetic field functions:
    % =========================================================================

    % Generate coil geometric setup from spreadsheet:
    % =========================================================================
    t1 = tic;
    disp('Reading "coilSetup" spreadsheet...')
    coilSetup = readtable('CoilSetup_MPEX.xlsx','Sheet','conf_1');
    disp('Reading complete!')

    % Convert to metric:
    % =========================================================================
    coilSetup.z       = coilSetup.z/1000;
    coilSetup.dz      = coilSetup.dz/1000;
    coilSetup.r_inner = coilSetup.r_inner/1000;
    coilSetup.r_outer = coilSetup.r_outer/1000;

    % Assignment of unit currents per power supply:
    % =========================================================================
    NC = numel(coilSetup.ps);
    for cc = 1:NC
        fieldName = ['C',num2str(cc)];
        coilCurrents.(fieldName) = 1;
    end

    % Create "coil" structure: 
    % =========================================================================
    disp('Creating "coil" structure...')
    [AllCoils] = CreateCoilStructure(coilSetup,coilCurrents);
    disp(['Complete! Elapsed time: ',num2str(toc(t1)),' s'])
    clearvars dum*

    if 1
        % =========================================================================
        % Plot Magnetic coils and current filaments:
        figure('color','w','Tag','coilSetup')
        hold on
        for ii = 1:numel(AllCoils)
            plot(AllCoils{ii}.zfil,+AllCoils{ii}.rfil,'r.');
            plot(AllCoils{ii}.zfil,-AllCoils{ii}.rfil,'r.');
        end
        % Formatting:
        set(gca,'FontName','times','FontSize',11)
        xlabel('z [m]','Interpreter','latex','FontSize',13)
        ylabel('r [m]','Interpreter','latex','FontSize',13)
        box on
        grid on
        axis image
        xlim([LZ_min,LZ_max])
        ylim([-1,1]*LR_max)
    end

    % Define the area to evaluate the fields at:
    % =========================================================================
    tic
    r1D = linspace(LR_min,LR_max,NR);
    z1D = linspace(LZ_min,LZ_max,NZ);

    % Calculate magnetic field basis::
    % =========================================================================
    dum1 = tic;
    for cc = 1:NC
        disp(['Coil ',num2str(cc), ', calculating magnetic field...'])
        coil{1} = AllCoils{cc};
        [br2D,bz2D,atheta2D,phi2D,z2D,r2D] = CalculateMagField(coil,z1D,r1D,'grid');
         b_hat{cc}.r = br2D;
         b_hat{cc}.z = bz2D;
         at_hat{cc}  = atheta2D;
         phi_hat{cc} = phi2D;
    end
    disp(['Complete! Elapsed time: ',num2str(toc(dum1)),' s'])
    clearvars dum*

    disp('Completed calculating magnetic field basis!')
    disp(['Calculation time: ',num2str(toc(t2)),' [s]'])
    
    % Save data:
    % =====================================================================
    fileName = ['MagneticField_Basis_',num2str(date)];
    save(fileName,'b_hat','at_hat','phi_hat','z2D','r2D','z1D','r1D','AllCoils');
    % Full calculation takes about 500 seconds: 200x100 grid
    
else
    fileName = 'MagneticField_Basis_05-Apr-2021';
    load(fileName);
end


%% Testing output data:
figure
hold on
surf(z1D,r1D,b_hat{1}.z','lineStyle','none')

figure
hold on
surf(z1D,r1D,phi_hat{1}','lineStyle','none')


bz  = zeros(size(b_hat{1}.z));
phi = zeros(size(b_hat{1}.z));
for cc = 1:numel(AllCoils)
    bz  = bz  + b_hat{cc}.z;
    phi = phi + phi_hat{cc};    
end

figure
surf(z1D,r1D,bz','lineStyle','none')

figure
contourf(z1D,r1D,phi',30,'lineStyle','none')