% Central cell coils:
% =========================================================================
n = 1;
coilType{n}.dz = 0.1;
coilType{n}.r1 = 0.8;
coilType{n}.dr = 0.1;
coilType{n}.r2 = coilType{n}.r1 + coilType{n}.dr;
coilType{n}.layers_z = 10;
coilType{n}.layers_r = 10;
coilType{n}.datum = 5;
coilType{n}.ps = 'centralCell';
coilType{n}.current = 0.5;

% Transition coils:
% =========================================================================
n = 2;
coilType{n}.dz = 0.1;
coilType{n}.r1 = 0.7;
coilType{n}.dr = 0.1;
coilType{n}.r2 = coilType{n}.r1 + coilType{n}.dr;
coilType{n}.layers_z = 10;
coilType{n}.layers_r = 10;
coilType{n}.datum = 5;
coilType{n}.ps = 'transition';
coilType{n}.current = 0.8;

% Mirror coils:
% =========================================================================
n = 3;
coilType{n}.dz = 0.1;
coilType{n}.r1 = 0.15;
coilType{n}.dr = 0.1;
coilType{n}.r2 = coilType{n}.r1 + coilType{n}.dr;
coilType{n}.layers_z = 20;
coilType{n}.layers_r = 20;
coilType{n}.datum = 5;
coilType{n}.ps = 'mirror';
coilType{n}.current = 4;

% Expander1 coils:
% =========================================================================
n = 4;
coilType{n}.dz = 0.1;
coilType{n}.r1 = 1.0;
coilType{n}.dr = 0.1;
coilType{n}.r2 = coilType{n}.r1 + coilType{n}.dr;
coilType{n}.layers_z = 10;
coilType{n}.layers_r = 10;
coilType{n}.datum = 5;
coilType{n}.ps = 'expander_1';
coilType{n}.current = -0.1;
