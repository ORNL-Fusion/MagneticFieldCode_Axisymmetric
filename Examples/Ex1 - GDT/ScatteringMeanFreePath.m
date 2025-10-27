% Plotting the mean free path for 90 degree coulomb scattering in a plasma
% 2019_12_08

clear all
close all

% =========================================================================
% Define parameter space:
nb = [0.5,1,2,4]*1e20; 
logA = 13;
E_a = logspace(2,4,100);

% =========================================================================
% Define mirror ratio:
R = 20;

% Calculate scattering length:
% =========================================================================
% Based on: A. A. Ivanov and V. V. Prikhodko, “Gas-dynamic trap: an
% overview of the concept and experimental results,” Plasma Phys. Control.
% Fusion, vol. 55, no. 6, p. 063001, 2013.
for s = 1:numel(nb)
    L_perp_ab{s} = ((E_a.^2)/(nb(s)*e_c*e_c))*((8*pi*e_0^2)/logA)*(log(R)/R);
end

% =========================================================================
% Plot results:
figure('color','w','position',[513.0000  309.0000  407.3333  308.6667])
figName = 'GDT length scaling';
C = {'k','r','bl','g','m','c','k:','r:'};
hold on
for s = 1:numel(nb)
    h(s) = plot(E_a*1e-3,L_perp_ab{s},C{s},'LineWidth',2);
    legendText{s} = ['$n_b$ = ',num2str(nb(s)*1e-20),' $\times10^{20}$ $[m^{-3}]$'];
end
% -------------------------------------------------------------------------
% Figure formatting:
box on
grid on
ylim([1e-1,1e3])
xlim([0.1,10])
set(gca,'XScale','log','Yscale','log','FontName','Times')
set(gca,'XTickLabel',[0.1,1,10])
set(gca,'YTickLabel',[0.1,1,10,100,1000])
title('$\lambda_{\perp}^{ii}$ $\frac{\ln R}{R}$ [m]','Interpreter','latex','FontSize',14)
ylabel('$\lambda_{\perp}^{ii}$ $\frac{\ln R}{R}$ [m]','Interpreter','latex','FontSize',14)
xlabel('$E_{i}$ [keV]','Interpreter','latex','FontSize',13,'FontName','times')
% -------------------------------------------------------------------------
% Legend:
hL = legend(h,legendText);
hL.Interpreter = 'latex';
hL.FontSize = 10;
hL.Location = 'NorthWest';
hL.Box = 'on';
% -------------------------------------------------------------------------
% Mirror ratio text:
textText = ['R = ',num2str(R)];
text(2,2e-1,textText,'Interpreter','latex','FontSize',13,'BackgroundColor','w','EdgeColor','k')

% -------------------------------------------------------------------------
% Rectangle to show region of interest, from 1 m to 10 m
if 0
    rectangle('Position',[min(E_a*1e-3),1,1e4,10],'EdgeColor','g','FaceColor','g'); % [x y w h]
end

saveFig = 0;
if saveFig
    saveas(gcf,figName,'tiffn')
end

