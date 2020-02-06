% In this script we provide an alternate method to determine the coil setup
% for a GDT plasma

% =========================================================================
DefineCoilTypes

% =========================================================================
% Define lengths of coil assemblies:
lf = 4; % Length of the fast ion confinement zone 
lt = 4; % Length of the transition zone to mirror magnet
lm = 0; % Length of the mirror zone magnet
le = 0; % Length of the expander zone magnet

% =========================================================================
% Define gaps between coil assemblies:
gap_ft = 0.8; % gap between fast ion confinement and transition zone
gap_tm = 0.01; % gap between transition and mirror coil
gap_me1 = 2.0; % gap between mirror coil and expander1 coil

% Derived quantities:
length_m2m   = (lf + 2*gap_ft) + 2*(lt + gap_tm);
length_total = length_m2m + 2*gap_me1;

%%
% =========================================================================
% Define coil spacing:
% -------------------------------------------------------------------------
% Central cell:
% -------------------------------------------------------------------------
n = 1;
fctr(n) = sqrt(2);
coilGap_desired = coilType{n}.r1*fctr(n)*0.6167;
num_coils(n) = 1 + round(lf/coilGap_desired);
zf_1 = -lf/2;
zf_2 =  lf/2;
coilPosition{n} = linspace(zf_1,zf_2,num_coils(n));

% -------------------------------------------------------------------------
% transition coils:
% -------------------------------------------------------------------------
n = 2;
fctr(n) = sqrt(2);
coilGap_desired = coilType{n}.r1*fctr(n)*0.6167;
num_coils(n) = 1 + round(lf/coilGap_desired);
zt_1 = lf/2 + gap_ft;
zt_2 = zt_1 + lt - gap_tm;
coilPosition{n} = linspace(zt_1,zt_2,num_coils(n));
coilPosition{n} = [-fliplr(coilPosition{n}),coilPosition{n}];

% -------------------------------------------------------------------------
% mirror coils:
% -------------------------------------------------------------------------
n = 3;
fctr(n) = sqrt(2);
coilGap_desired = coilType{n}.r1*fctr(n)*0.6167;
num_coils(n) = 1 + round(lm/coilGap_desired);
zm_1 = zt_2 + gap_tm;
zm_2 = zm_1 + lm;
coilPosition{n} = linspace(zm_1,zm_2,num_coils(n));
coilPosition{n} = [-fliplr(coilPosition{n}),coilPosition{n}];

% -------------------------------------------------------------------------
% Expander coils:
% -------------------------------------------------------------------------
n = 4;
fctr(n) = sqrt(2);
coilGap_desired = coilType{n}.r1*fctr(n)*0.6167;
num_coils(n) = 1 + round(lm/coilGap_desired);
ze_1 = zm_2 + gap_me1;
ze_2 = ze_1;
coilPosition{n} = linspace(ze_1,ze_2,num_coils(n));
coilPosition{n} = [-fliplr(coilPosition{n}),coilPosition{n}];

% =========================================================================
% Produce coil setup structure:
k = 1;
for n = 1:numel(coilPosition) % For all coil types
    for ii = 1:length(coilPosition{n})
        coilSetup.ps{k} = coilType{n}.ps;
        coilSetup.datum(k) = coilType{n}.datum;
        coilSetup.z(k) = coilPosition{n}(ii);
        coilSetup.dz(k) = coilType{n}.dz;
        coilSetup.r_inner(k) = coilType{n}.r1;
        coilSetup.r_outer(k) = coilType{n}.r2;
        coilSetup.layers_z(k) = coilType{n}.layers_z;
        coilSetup.layers_r(k) = coilType{n}.layers_r;
        k = k + 1;
    end
end

% =========================================================================
% Sort coils:
[~,b]   = sort(coilSetup.z);
coilSetup.ps = coilSetup.ps(b);
coilSetup.datum = coilSetup.datum(b);
coilSetup.z = coilSetup.z(b);
coilSetup.dz = coilSetup.dz(b);
coilSetup.r_inner = coilSetup.r_inner(b);
coilSetup.r_outer = coilSetup.r_outer(b);
coilSetup.layers_z = coilSetup.layers_z(b);
coilSetup.layers_r = coilSetup.layers_r(b);

%% Test
% =========================================================================
% Here we plot the location of the finite size coils and we use the
% filamentary current loop model to calculate the magnetic field profile to
% provide a quick check for the correct position of the coils

% Test prompt:
InputStructure.prompt = {['Test "coilSetup"?, Yes [1], No [0]']};
testCoils = GetUserInput(InputStructure);

% =========================================================================
if testCoils
    disp('Testing "coilSetup"...')
    % ---------------------------------------------------------------------
    % Assignment of currents per power supply:
    coilCurrents.expander_1  = +0e3; % Expander coils
    coilCurrents.expander_2  = -0.5e3; % Expander coils
    coilCurrents.mirror      = +15e3; % Mirror coils
    coilCurrents.centralCell = +5e3; % Central cell coils
    coilCurrents.transition  = +8e3; 

    % ---------------------------------------------------------------------
    % Produce coil strucutre:
    [coil] = CreateCoilStructure(coilSetup,coilCurrents);

    % ---------------------------------------------------------------------
    % Equation for the on-axis B field of a current loop:
    f = @(X) 1./( (1 + (X.^2)).^(3/2) );

    % ---------------------------------------------------------------------
    % Calculate the magnetic field:
    z = linspace(-16,16,1e3);
    fm = 0;
    mu_0 = 4*pi*1e-7;
    for ii = 1:numel(coil)
        % Offsetted argument:
        arg{ii} = (z - coil{ii}.z)/coil{ii}.r;
        % Evaluate magnetic field of single current loop:
        fn{ii} = 0.5*f(arg{ii})*coil{ii}.nfil*coil{ii}.current*mu_0/coil{ii}.r;
        % Add the contributios from all coils:
        fm = fm + fn{ii};
    end

    % ---------------------------------------------------------------------
    % Plot data:
    figureName = 'Testing Coil setup';
    figure('color','w'); 
    hold on
    for ii = 1:numel(coil)
        plot(z,fn{ii})
    end
    for ii = 1:numel(coil)
        plot(coil{ii}.zfil,coil{ii}.rfil,'r.');
    end
    plot(z,fm,'linewidth',2)
    ylim([0,20])
    xlim([-10,10])
    box on
    grid on
    xlabel('z [m]','Interpreter','latex','FontSize',14)
    
    disp('Testing complete')
end
