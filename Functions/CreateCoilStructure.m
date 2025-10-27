function [coil] = CreateCoilStructure(coilGeometry,coilCurrents)
% #########################################################################
% Created 2019_12_09, JF Caneses
% =========================================================================
% CreateCoilStructure creates a structure called "coil" which contains all
% the geometric information, currents and location of coil filaments based
% on two inputs "coilGeometry" and "coilCurrents"
% =========================================================================
%                              INPUT:
% =========================================================================
% coilGeometry:
% -------------------------------------------------------------------------
% Table that contains all the information to build a
% physical coil. 
% An example is shown below:
%     coil_index         ps          datum      z      dz     r_inner    r_outer    layers_z    layers_r    setup    dateEffective    comment
%     __________    _____________    _____    _____    ___    _______    _______    ________    ________    _____    _____________    _______ 
%      1            'expander_1'     5           -6    0.1     0.8        0.9       10          10            1      Add date         None    
%      2            'mirror'         5        -5.25    0.1    0.15       0.25       20          20          NaN      NaN              NaN    
%      3            'centralCell'    5         -4.5    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%      4            'centralCell'    5        -3.75    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%      5            'centralCell'    5           -3    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%      6            'centralCell'    5        -2.25    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%      7            'centralCell'    5         -1.5    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%      8            'centralCell'    5        -0.75    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%      9            'centralCell'    5            0    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%     10            'centralCell'    5         0.75    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%     11            'centralCell'    5          1.5    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%     12            'centralCell'    5         2.25    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%     13            'centralCell'    5            3    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%     14            'centralCell'    5         3.75    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%     15            'centralCell'    5          4.5    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
%     16            'mirror'         5         5.25    0.1    0.15       0.25       20          20          NaN      NaN              NaN    
%     17            'expander_2'     5            6    0.1     0.8        0.9       10          10          NaN      NaN              NaN    
% The columns "setup", "dateEffective" and "comment" are optional
% All other columns are required
% -------------------------------------------------------------------------
% coilCurrents:
% -------------------------------------------------------------------------
% Structure that holds the current associated with each power supply "ps"
% described in the table above.
% An example showing the fields and values of each field associated with
% the table shown above:
% coilCurrents =
%      expander_1: 0
%      expander_2: -500
%          mirror: 8000
%     centralCell: 2000
% =========================================================================
%                               OUTPUT:
% =========================================================================
% coil:
% Structure that contains the following fields:
% coil{1}
%      current: 0
%        datum: 5
%           dr: 0.1000
%        drfil: 0.0100
%         dum1: 0.0100
%           dz: 0.1000
%        dzfil: 0.0100
%     layers_r: 10
%     layers_z: 10
%         nfil: 100
%           ps: {'expander_1'}
%            r: 0.8500
%           r1: 0.8000
%           r2: 0.9000
%            z: -6
%         zfil: [10x10 double]
%         rfil: [10x10 double]
% =========================================================================
% #########################################################################


% START OF FUNCTION:
% =========================================================================
% Define coil setup and geometry:
dum1 = who;
z         = coilGeometry.z;
dz        = coilGeometry.dz;
r1        = coilGeometry.r_inner;
r2        = coilGeometry.r_outer;
layers_z  = coilGeometry.layers_z;
layers_r  = coilGeometry.layers_r;
datum     = coilGeometry.datum; % Defines how "z" is interpreted
ps        = coilGeometry.ps;
current  = []  ; % Variable to hold the value of the coil current

% =========================================================================
% Derived quantities
dr    = r2-r1;
r     = (r1+r2)/2;
nfil  = layers_z.*layers_r; % Total number filaments per coil
dzfil = dz./layers_z;
drfil = dr./layers_r;
dum2 = who;

% =========================================================================
% Record the name of all the variables that will be part of the "coil"
% structure:
% For each physical coil, a structure "coil" will be created with fields
% given by "coil_fieldnames"
coil_fieldnames = setdiff(dum2,dum1);

% =========================================================================
% Definition of coil datum:

% 7=====8=====9
% |     |     |
% 4=====5=====6
% |     |     |
% 1=====2=====3

% =========================================================================
% Assign coil currents to an array
psType = unique(ps);
for ii = 1:numel(psType)
    dum1 = find(strcmp(ps,psType{ii}));
    try
        current(dum1) = coilCurrents.(psType{ii});
    catch
        current(dum1) = 0;
    end
end
clearvars dum* ii

% =========================================================================
% Structure that holds coil geometric setup and currents:
for ii = 1:numel(z)
    for jj = 1:numel(coil_fieldnames)
        dum1 = eval(coil_fieldnames{jj});
        coil{ii}.(coil_fieldnames{jj}) = dum1(ii);
    end
end
clearvars ii jj dum*

% =========================================================================
% Create the current filaments for each coil:
% Each filament is in fact a single current loop
for ii = 1:numel(coil)
    z_ii     = coil{ii}.z;
    dz_ii    = coil{ii}.dz;
    dzfil_ii = coil{ii}.dzfil;
    r_ii     = coil{ii}.r;
    dr_ii    = coil{ii}.dr;
    drfil_ii = coil{ii}.drfil;

    dum1 = (dzfil_ii/2):dzfil_ii:dz_ii;
    dum2 = (drfil_ii/2):drfil_ii:dr_ii;
    [zfil_rel,rfil_rel] = meshgrid(dum1,dum2);
    clearvars dum*

    switch coil{ii}.datum
        case 1
        case 2
        case 3
        case 4
        case 5
            coil{ii}.zfil = zfil_rel + z_ii - dz_ii/2;
            coil{ii}.rfil = rfil_rel + r_ii - dr_ii/2;
        case 6
    end
end
clearvars ii jj
end

