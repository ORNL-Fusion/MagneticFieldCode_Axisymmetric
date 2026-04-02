% Step_1a_CreateFieldsAndValidate
% Using the magnetic field basis created in step_0, we compute the magnetic
% field scenarios.
% In addition, we save the magnetic field profiles.

clear all
close all
clc

saveData = 1;
saveFig  = 1;

% Load magnetic field basis:
% =========================================================================
fileName = 'Step_0_MagneticField_Basis_05-Apr-2021';
load(fileName);
    
% Select file containing the coil currents:
% =========================================================================
targetFolder = '../MagneticFieldProfiles_2021_04_15';
fileName = '/MPEX Coil Parameters - 2021-04-15';      
date = '2021_04_15';

% Load coil currents for all scenarios:
% =========================================================================
currents = readmatrix([targetFolder,fileName],'Sheet','MPEX Preliminary Design - Field','Range','K5:AH27');
NP = size(currents,2);

% Compute scenarios:
% =========================================================================
NC = numel(AllCoils);
for pp = 1:NP
    Bz{pp}  = zeros(size(b_hat{1}.z));
    phi{pp} = zeros(size(b_hat{1}.z));
    I{pp} = currents(:,pp);
    for cc = 1:NC
        a = 1;
        if (sum(cc == [7:17]) == 1)
            a = 63.5/65;
        end
        Bz{pp}  = Bz{pp}  + a*I{pp}(cc)*b_hat{cc}.z;
        phi{pp} = phi{pp} + a*I{pp}(cc)*phi_hat{cc};
    end
    Bz_axis{pp} = Bz{pp}(:,1);
end

% Comparing magnetic field calculations:
% =========================================================================
figure('color','w')
for pp = 1:NP
    subplot(5,5,pp)
    hold on
    hB(1) = plot(z1D,Bz_axis{pp},'k','lineWidth',2);

    xlim([z1D(1),z1D(end)])
    box on
    title(['scenario: ',num2str(pp)])
    xlabel('z [m]','interpreter','latex')
end

set(gcf,'position',[0.0010    0.0410    1.2800    0.6073]*1000)

% Save figure:
% =========================================================================
if saveFig
    fileName = ['Step_1a_FieldScenarios_',date,'.pdf'];
    exportgraphics(gcf,fileName,'Resolution',300) 
end

% Select scenario and select subset:
% =========================================================================
zDump   = -2;
zTarget = 8; 
rng = find(z1D >= zDump & z1D <= zTarget);
pp = 14;
figure('color','w')
plot(z1D(rng),Bz_axis{pp}(rng),'k','lineWidth',2)

% Create data for PICOS:
% =========================================================================
NX = 200;
picos.x = linspace(zDump,zTarget,NX)';
picos.Bx = interp1(z1D,Bz_axis{pp},picos.x); 

hold on
plot(picos.x,picos.Bx,'r.','lineWidth',1)

% Extract reference magnetic field
% =========================================================================
x0 = 0;
i_x0 = find(picos.x >= 0,1);
Bx_0 = picos.Bx(i_x0);

% Normalized magnetic field
% =========================================================================
picos.Bx_norm = picos.Bx./Bx_0;
picos.n_norm = ones(size(picos.Bx_norm));
picos.Tpar_norm = ones(size(picos.Bx_norm));
picos.Tper_norm = ones(size(picos.Bx_norm));

% Plot normalized magnetic field:
% =========================================================================
figure('color','w')
plot(picos.x,picos.Bx_norm,'k-','lineWidth',2)
xlabel('x [m]','interpreter','latex','fontSize',15)
ylabel('normalized B','interpreter','latex','fontSize',15)
title(['$x_0 = $',num2str(x0),' , $B_0 = $',num2str(Bx_0),' [T]'],'interpreter','latex','fontSize',15)

% Save normalized magnetic field data:
% =========================================================================
if saveFig
    fileName = ['MPEX_B_norm_PICOS_scenario_',num2str(pp),'.pdf'];
    exportgraphics(gcf,fileName,'Resolution',300) 
end

% Save data to text file:
% =========================================================================
if saveData
    fileName = 'MPEX_B_norm_PICOS_scenario_14';
    f1 = [picos.Bx_norm];
    save([fileName,'.txt'],'f1','-ascii');
    
    fileName = 'MPEX_x_B_PICOS_scenario_14';
    f1 = [picos.x];
    save([fileName,'.txt'],'f1','-ascii');
     
%     fileName2 = figureName2  
%     f2 = [picos.n_norm];
%     save([fileName2,'.txt'],'f2','-ascii');
%     
%     fileName3 = figureName3  
%     f3 = [picos.Tpar_norm];
%     save([fileName3,'.txt'],'f3','-ascii');
%     
%     fileName4 = figureName4  
%     f4 = [picos.Tper_norm];
%     save([fileName4,'.txt'],'f4','-ascii');    
end