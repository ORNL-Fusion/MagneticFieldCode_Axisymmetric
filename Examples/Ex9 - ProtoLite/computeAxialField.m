function B_total = computeAxialField(data1, Z, mu0, I)
% computeAxialField computes the total axial magnetic field along the Z-axis.

% Inputs:
% - data1: Table containing coil parameters
% - Z: Axial grid vector
% - mu0: Magnetic permeability of free space
% - I: Vector of currents for each coil

% Output:
% - B_total: Axial magnetic field on axial grid

B_total = zeros(size(Z)); % Initialize total field
lim = height(data1);      % Number of coils

for i = 1:length(Z)
    for j = 1:lim
        Xj = data1.z(j) - data1.dz(j)/2; % Axial start position for this coil, Z: center loc
        dR = (data1.r_outer(j) - data1.r_inner(j)) / data1.layers_r(j);
        dZ = data1.dz(j) / data1.layers_z(j);
        for l = 1:data1.layers_z(j)            % Iterate axial layers
            R = data1.r_inner(j);              % Start at inner radius
            for k = 1:data1.layers_r(j)        % Iterate radial layers
                B = mu0 * R^2 * I(j) / (2 * (( (Z(i) - Xj)^2 + R^2)^(3/2)));
                B_total(i) = B_total(i) + B;   % Add contribution to total field
                R = R + dR;                    % Increment radius for next layer
            end
            Xj = Xj + dZ;                      % Increment axial position
        end
    end
end
end
