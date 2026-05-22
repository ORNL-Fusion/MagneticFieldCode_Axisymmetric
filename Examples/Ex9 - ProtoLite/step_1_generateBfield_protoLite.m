% main_program.m - Main script to compute and plot magnetic flux and axial field

clc;
clear all;

% ----- Read in coil data and assign coil parameters ----------
data1 = readtable("Coil_setup_I_locs_8inbackward.xlsx"); % Import coil setup data    
%data1 = readtable("Coil_setup_I_locs.xlsx"); % Import coil setup data
% Coil parameters:
    nR = data1.layers_r;       % Number of radial layers
    nZ = data1.layers_z;       % Number of axial layers (turns)
    R1 = data1.r_inner;        % Inner radius
    R2 = data1.r_outer;        % Outer radius
    wZ = data1.dz;             % Coil thickness (axially)
    X = data1.z;               % Axial location of coil start
    ps = data1.ps;             % String names for power supplies
    z_start = -0.25; % begining of the dump chamber
    z_end = 3.5; % end of the device
% ----- Map the string current names to numeric current values ----------
    mapping = containers.Map({'PS1', 'TR2', 'TR1'}, [1575, 220, 1575]); % Conversion map for currents
    I = zeros(size(ps));       % Preallocate current vector
for idx = 1:length(ps)
    I(idx) = mapping(ps{idx}); % Map string names to numerical currents
end

% ----- Define grid for flux calculation (radial & axial) -----------------
    y = linspace(0.001, 0.13, 150); % Radial grid (meters)
    Z = linspace(z_start, z_end, 500);  %linspace(0, 2.8616, 300);     % Axial grid (meters)
    [Y, Z_matrix] = meshgrid(y, Z); % Create meshgrid for radial & axial grid
    dist = sqrt((Z_matrix).^2 + (Y).^2);

% ----- Compute the magnetic flux function (A_phi_total) -----------------
    A_phi_total  = computeMagneticFlux(data1, Y, Z_matrix, I);


    % axial limit for B-field:
    ZB = linspace(z_start, z_end, 500);
    %YB = linspace(1e-3, 0.05, 5); % this radial coordinate really doesnt matter
% ----- Compute the axial magnetic field profile --------------------------
    B_total = computeAxialField(data1, ZB, 4 * pi * 1e-7, I);
b_data = [ZB',B_total'];
% ----- Plot axial field profile with vertical coil-location lines --------
    figure;
    subplot(2, 1, 1);
    plot(ZB, B_total, 'LineWidth', 2.5); % Plot axial field
    % datawrite = [ZB', B_total'];
    % writematrix(datawrite, 'B_field_1575_150_500.xlsx')
    hold on;
for j = 1:length(X)
    xline(X(j), 'r--', ['Coil ' num2str(j)]); % Mark coil positions
end
    grid on;
    hold off;
    xlabel('Z (m)');
    ylabel('Axial magnetic field');
    title('Axial Magnetic Field Profile');
    %ylim([0,1]);
    % xlim([0,4]);

% ----- Plot the device geometry and flux contours ------------------------
    subplot(2,1,2);
    Find_LCFS(A_phi_total, y, Z, X, data1,z_start) %, Y, Z_matrix);
%plotDeviceBoundary(A_phi_total, y, Z, X, Y, Z_matrix);
    hold on
% Plotting coil locations as per their positions
 for j=1:length(X)
    dR = R2(j)-R1(j);
    rectangle('Position',[(X(j)-wZ(j)/2), R1(j), wZ(j), dR],'FaceColor',[1 0 0],'EdgeColor', 'r', 'LineWidth', 1);
    rectangle('Position',[(X(j)-wZ(j)/2), -R2(j), wZ(j), dR],'FaceColor',[1 0 0],'EdgeColor', 'r', 'LineWidth', 1);
 end
 hold off
    % figure;
    % Find_LCFS(A_phi_total, y, Z, X, Y, Z_matrix);