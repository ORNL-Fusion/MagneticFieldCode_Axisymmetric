function [Br2D,Bz2D,Atheta2D,phi2D,z2D,r2D] = CalculateMagField(coil,z1D,r1D,evalType)
% #########################################################################
% Created 2019_12_09, JF Caneses
% =========================================================================
% CALCULATEMAGFIELD:
% Calculates the magnetic field in 2D space using filametary current loops
% defined by the structure coil
% The structure "coil" needs to be created by the function
% "CreateCoilStructure"
% =========================================================================
%                               INPUT:
% =========================================================================
% coil:
% Run function "CreateCoilStructure" to create "coil"
% The "coil" structure contains the following fields:
% coil{ii}
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
% where ii refers to the "ith" coil in the geometry
% -------------------------------------------------------------------------
% z1D and r1D:
% one-dimensional arrays that define the 2D space at which the magnetic
% field will be computee
% #########################################################################

% START OF FUNCTION:
% =========================================================================
% Define the area to evaluate the fields at:
if strcmpi(evalType,'grid')
    % Generate a 2D grid for evaluating field
    [r2D,z2D] = meshgrid(r1D,z1D);
elseif strcmpi(evalType,'contour')
    % Generate a contour over which field is to be evaluated such as along
    % a boundary
    r2D = r1D;
    z2D = z1D;
else
    error('Please define "evalType"')
end

% =========================================================================
% Calculate the magnetic field and magnetic vector potential:
Br2D = 0;
Bz2D = 0;
Atheta2D = 0;
for ii = 1:numel(coil)
    current = coil{ii}.current;
    Br_n{ii} = 0;
    Bz_n{ii} = 0;
    Atheta_n{ii} = 0;
    % Calculate field produced by each current filament loop:
    for jj = 1:coil{ii}.nfil
        [Br0,Bz0,Atheta0] = bfield_circular_coil_analytic(coil{ii}.rfil(jj),coil{ii}.zfil(jj),r2D,z2D);
        Br_n{ii}     = Br_n{ii}     + Br0*current;
        Bz_n{ii}     = Bz_n{ii}     + Bz0*current;
        Atheta_n{ii} = Atheta_n{ii} + Atheta0*current;
    end
    % Add the contribution from all coils:
    Br2D = Br2D + Br_n{ii};
    Bz2D = Bz2D + Bz_n{ii};
    Atheta2D = Atheta2D + Atheta_n{ii};
end

% =========================================================================
% Calculate the magnetic flux:
phi2D = 2*pi*Atheta2D.*r2D;

end

