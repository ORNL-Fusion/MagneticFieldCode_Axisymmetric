function A_phi_total = computeMagneticFlux(data1, Y, Z_matrix, I)
% computeMagneticFlux computes the azimuthal magnetic potential A_phi_total.

% Inputs:
% - data1: Table with coil parameters (e.g., r_inner, r_outer, dz, z, ...)
% - y: Radial grid vector
% - Z: Axial grid vector
% - I: Current vector for coils

% Outputs:
% - A_phi_total: Total magnetic potential
% - Y, Z_matrix: Meshgrid arrays for plotting


A_phi_total = zeros(size(Y));   % Initialize total potential
lim = height(data1); % Number of coils (data rows)

for j = 1:lim
    Xj = data1.z(j) - data1.dz(j)/2; % Axial start position for this coil, Z-is center of the coil
    dR = (data1.r_outer(j) - data1.r_inner(j)) / data1.layers_r(j);
    dZ = data1.dz(j) / data1.layers_z(j);      % Incremental size per layer
    for l = 1:data1.layers_z(j)                % Iterate over axial layers
        R = data1.r_inner(j);                  % Start at inner radius
        for m = 1:data1.layers_r(j)            % Iterate over radial layers
            Dr = (R + Y).^2 + (Z_matrix - Xj).^2; % Compute distance matrix
            k2 = 4 * R * Y ./ Dr;
            k2(Dr == 0) = 0;                   % Prevent division by zero
            k = sqrt(k2);
            [K, E] = ellipke(k2);              % Compute elliptic integrals
            term1 = sqrt(R ./ Y);
            term1(Y <= 0) = 0;                 % Avoid invalid sqrt operations
            term2 = (2./k - k).*K - (2./k).*E;
           %
            A_phi = 2e-7 * I(j) .* term1 .* term2; % Compute A_phi contribution
            A_phi(isnan(A_phi)) = 0;           % Replace NaNs with zeros
            A_phi_total = A_phi_total + A_phi; % Update total potential
            R = R + dR;                        % Increment radius for next layer
        end
        Xj = Xj + dZ;                          % Increment axial position
    end
end
end
