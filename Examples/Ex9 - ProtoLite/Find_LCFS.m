function Find_LCFS(A_phi_total, y, Z, X, data1,z_start)
    % plotDeviceBoundary plots the physical device geometry and the magnetic flux lines.

    %% Step 1: Physical Device Geometry
        in2m = 0.0254; % Inch to meter conversion factor

    % Define dump plate
    rDump = 11*in2m/2 + in2m; % Im adding 2" to radius of dump plate %0.178002;
    zDump1 = z_start; %-0.25; %0.0;
    
    rDumpplate = 11*in2m/2;
    rDumpplate2 = 0.0;
    zDumplplate = X(1) - 8*in2m; % John email 1/29/2026, the current dump plate position
    % 2.3.2026 setup   %15*in2m; %
    zDumplplate2 = zDumplplate - 0.0005;
    Dumpplate.z = [zDumplplate,zDumplplate,zDumplplate2,zDumplplate2];
    Dumpplate.r = [rDumpplate2,rDumpplate,rDumpplate,rDumpplate2];
    
    
    zDump2 = X(1) - data1.dz(1); %0.67;

    % Define vacuum vessel
    rVac = 0.0745;
    zVac1 = zDump2;
    zVac2 = X(8) + 11.5*in2m; % 11.17.25, John's new drawing of Proto_lite 10.15.25  %3.2731;

    % Define end chamber
    rEnd = 11.5*in2m/2;  % 0.1255;
    zEnd1 = zVac2;
    zEnd2 = zEnd1+ 19.5*in2m;  %3.8351;

    % Define helicon window dimensions
    L = 11.8 * in2m; % Window length
    r1 = 4.75*in2m/2; %(5.426 - (2 * 0.25)) * in2m / 2;
    r2 = rVac;
    center_w = X(4) - 4.938*0.0254;
    z2 = center_w + L/2;
    z1 = center_w - L/2;
    % z1 = 0.5 * (X(3) + X(4)) - L / 2;
    % z2 = 0.5 * (X(3) + X(4)) + L / 2;
    heliconWindow.r = [r2, r1, r1, r2];
    heliconWindow.z = [z1, z1, z2, z2];

    % Define limiter dimensions
    limiterLength = 0.1778;  % 30e-2;
    %limiterWidth = 2.55e-3;
    z2 = X(4) + 10.129*0.0254;    % z1 + limiterLength;
    z1 = z2 - limiterLength;   %heliconWindow.z(end) + 1e-2;
    r1 = 0.120904/2; %2.5 * in2m - limiterWidth;
    r2 = rVac;
    limiter.r = [r2, r1, r1, r2];
    limiter.z = [z1, z1, z2, z2];

    % Skimmer 1: GD: Needs to be updated, removing for now
    r1 = (7 / 100) / 2;
    r2 = rVac; % GD: radius of the vacuum vessel
    z1 = X(4) + 18.434*0.0254 - 0.5e-2;   %X(5) - 0.075 - 0.02 - 0.5e-2;
    z2 =  z1 + 0.5e-2;                    %X(5) - 0.075 - 0.02 + 0.5e-2;
    skimmer1.r = [r2, r1, r1, r2];
    skimmer1.z = [z1, z1, z2, z2];
    % define Target location
        z1 =  X(8) - 5.438*in2m;  %(X(8)+ X(7))/2;   %2.8616; % targetlocation
        z2 = z1+ 0.01;
        r1 = 0.02; %radius of the target cover %0.01; % approximately 1 cm radius
        r2 = 0;
        Target.z = [z1,z1,z2,z2];
        Target.r = [r2,r1,r1,r2];

    % Combine parts into a unified boundary
        vessel_U.r = [rDump, rDump, rVac, rVac, rEnd, rEnd, 0];
        vessel_U.z = [zDump1, zDump2, zVac1, zVac2, zEnd1, zEnd2, zEnd2];

    % Plot geometry
        plot(vessel_U.z, vessel_U.r, 'k-','LineWidth',1);
        hold on;
        plot(heliconWindow.z, heliconWindow.r, 'b-','LineWidth',1);
        plot(limiter.z, limiter.r, 'g-','LineWidth',1);
        plot(vessel_U.z, -vessel_U.r, 'k-','LineWidth',1);
        plot(heliconWindow.z, -heliconWindow.r, 'b-','LineWidth',1);
        plot(limiter.z, -limiter.r, 'g-','LineWidth',1);
        plot(skimmer1.z, skimmer1.r, 'k-','LineWidth',1);
        plot(skimmer1.z, -skimmer1.r, 'k-','LineWidth',1);
        plot(Target.z,Target.r, 'k-', 'LineWidth',2);
        plot(Target.z,-Target.r, 'k-', 'LineWidth',2);
        plot(Dumpplate.z,Dumpplate.r, 'k-', 'LineWidth',2);
        plot(Dumpplate.z,-Dumpplate.r, 'k-', 'LineWidth',2);
        %
        fill(heliconWindow.z, heliconWindow.r, 'blue');
        fill(heliconWindow.z, -heliconWindow.r, 'blue');
        fill(limiter.z, limiter.r, 'green');
        fill(limiter.z, -limiter.r, 'green');
        %
        yline(0, '--k','LineWidth',1);
        %xlim([0, 3.5]);
        %ylim([-0.2, 0.2]);

    %% Step 2: Magnetic Flux - Adjusted Filtering for Limiter
        dataLog = log10(abs(2 * pi * y' .* A_phi_total')); % Logarithmic flux data

    % Generate contour matrix
        C = contourc(Z, y, dataLog, 150); % Extract contour segments
        idx = 1;

    % Initialize contour level limit
        contourLevel = -5.2; % Default minimum contour level
        stepSize = 0.1; % Increment step for contour level adjustments
        maxRadialLocation = -Inf; % To track the maximum radial location of filtered segments

    % Iteratively filter contours and adjust contour level
    while idx < size(C, 2)
        nPoints = C(2, idx); % Number of points in the segment
        xPoints = C(1, idx+1 : idx+nPoints); % Z-coordinates of the points
        yPoints = C(2, idx+1 : idx+nPoints); % R-coordinates of the points

        % Segment-based boundary checks (your specific boundary values)
    if any(xPoints >= zDump1 & xPoints <= zDump2 & yPoints >= rDump)
    elseif any(xPoints >= zVac1 & xPoints <= zVac2 & yPoints >= rVac)
    elseif any(xPoints >= zEnd1 & xPoints <= zEnd2 & yPoints >= rEnd)
    elseif any(xPoints >= heliconWindow.z(1) & xPoints <= heliconWindow.z(3) & ...
                   yPoints >= heliconWindow.r(2))
    elseif any(xPoints >= limiter.z(1) & xPoints <= limiter.z(3) & ...
               yPoints >= limiter.r(2))
%        Update maximum radial extent from valid segments
        maxRadialLocation = max(maxRadialLocation, max(abs(yPoints)));
    else
        contourLevel = C(1, idx); % Update the maximum allowed contour level
    end

        idx = idx + nPoints + 1; % Move to the next segment
    end

    % Adjust contour level to reach limiter radius if necessary
    % heliconlim=heliconWindow.r(2);    
    limiterRadius = limiter.r(2); % Minimum limiter radius
    while maxRadialLocation < limiterRadius
        %while maxRadialLocation < limiterRadius
        contourLevel = contourLevel + stepSize; % Increase contour level
        %fluxLevels = linspace(-12, contourLevel, 60); % Update flux levels
        maxRadialLocation = limiterRadius; % Assume it now touches the limiter
    end

    % Plot filtered contours
    fluxLevels = linspace(-12, contourLevel, 30); % Final flux levels    
    %fluxLevels = linspace(-12, -1, 100); % Final flux levels
        contour(Z, y, dataLog, fluxLevels, 'LineWidth', 1); % Upper half
        contour(Z, -y, dataLog, fluxLevels, 'LineWidth', 1); % Lower half
        colormap("jet");
        colorbar;
        title('Magnetic Flux and Device Geometry');
        xlabel('Z (m)');
        ylabel('Radial Distance (m)');
        grid on;
        hold off;
end

