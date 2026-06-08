function [C_min,C_max, max_p_hel, max_p_lim] = Find_BDXY(A_phi_total, y, Z, X, data1,z_start)

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

    mpex_bd = [vessel_U.z;vessel_U.r];
    tmp = [heliconWindow.z;heliconWindow.r];
    tmp2 = [limiter.z;limiter.r];
    mpex_bd(1,1) = 0;
    mpex_bd = [[0;0] mpex_bd(:,1:3) tmp tmp2 mpex_bd(:,4:end)];

    %% Step 2: Magnetic Flux - Adjusted Filtering for Limiter
    dataLog = log10(abs(2 * pi * y' .* A_phi_total')); % Logarithmic flux data

    n_C = 10;
    C_min = min(dataLog(:));
    C_max = max(dataLog(:));

    while (C_min-C_max)/C_max>1e-8
        c_line = linspace(C_min,C_max,n_C);

        C = contourc(Z, y, dataLog,c_line); % Extract contour segments
        c_ind = 1;
        cl = [];
        while c_ind<size(C,2)
            cl_c = [C(1,c_ind);0];
            cp = C(:,c_ind+1:c_ind+C(2,c_ind));
            c_ind = c_ind+C(2,c_ind)+1;
            if min(cp(1,:))<zDump2
                if max(cp(1,:))>zEnd1
                    if min(cp(2,:))<rVac
                        t_ind = find(cp(1,:)>=rVac & cp(1,:)<=zEnd1);
                        [max_v,max_ind] = max(cp(2,t_ind));
                        if max_v<rVac
                            t_ind = find(cp(1,:)>=heliconWindow.z(2) & cp(1,:)<=heliconWindow.z(3));
                            [max_v,max_ind] = max(cp(2,t_ind));
                            if max_v<heliconWindow.r(2)
                                max_p_hel = cp(:,t_ind(max_ind));
                                t_ind = find(cp(1,:)>=limiter.z(2) & cp(1,:)<=limiter.z(3));
                                [max_v,max_ind] = max(cp(2,t_ind));
                                if max_v<limiter.r(2)
                                    max_p_lim = cp(:,t_ind(max_ind));
                                    t_ind = find(cp(1,:)>=skimmer1.z(2) & cp(1,:)<=skimmer1.z(3));
                                    [max_v,max_ind] = max(cp(2,t_ind));
                                    if max_v<skimmer1.r(2)
                                        cl_c(2) = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            cl = [cl cl_c];
        end
        us_ind = find(cl(2,:)==1);
        if isempty(us_ind)
            max_p_hel=[0;0];
            max_p_lim = [0;0];
            break
        end
        C_min = cl(1,us_ind(end));
        if us_ind(end)<n_C
            C_max = cl(1,us_ind(end)+1);
        end
        if isempty(cl)
            max_p_hel=[0;0];
            max_p_lim = [0;0];
            break
        end
    end
    




