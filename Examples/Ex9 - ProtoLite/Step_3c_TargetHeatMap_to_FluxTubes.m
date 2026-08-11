%% Step_2_MapTargetHeatFlux_ByPsiFluxTubes_3D.m
% Constant-psi target heat-flux mapping with Proto-Lite geometry.
% black = vessel/skimmers, blue = helicon window, green = limiter

clc
clear
close all
format long

%% ========================= USER INPUTS ================================

targetMapFile = "IR_TargetHeatMap.mat";
bfieldFile    = "bfield_protoLite.nc";
geometryFile  = "Coil_setup_I_locs_8inbackward.xlsx";

z_start  = -0.25;
z_end    =  3.50;
z_target =  2.30;

dz_map = -2e-3;

q_threshold_MWm2 = 0.005;
max_tubes_to_plot = 2000;
tube_selection_mode = "weighted";

outputFile = "TargetHeatFlux_3D_PsiFluxTubes.mat";

%% ===================== DEVICE GEOMETRY ================================

vessel_top_z = [-0.25  0.09  0.09  2.73  2.73  3.23  3.23  3.50];
vessel_top_R = [ 0.165 0.165 0.075 0.075 0.145 0.145 0.145 0.145];

% Blue helicon window centered near z = 1 m
helicon_window_z = [0.88 1.18];
helicon_window_R = 0.063;
helicon_window_thick = 0.010;

% Green limiter / aperture with finite axial width
limiter_z = [1.28 1.45];
limiter_R = 0.063;
limiter_thick = 0.010;

% Black skimmers / apertures
left_aperture_z  = 0.00;
left_aperture_R  = 0.140;

mid_aperture_z   = 1.60;
mid_aperture_R   = 0.035;

right_aperture_z = 2.32;
right_aperture_R = 0.020;

%% ===================== READ COIL GEOMETRY =============================

data1 = readtable(geometryFile);

R1 = data1.r_inner;
R2 = data1.r_outer;
wZ = data1.dz;
X  = data1.z;
ps = data1.ps;

fprintf("Read geometry from %s\n", geometryFile);
fprintf("Number of coils = %d\n", numel(X));

%% ===================== LOAD TARGET HEAT MAP ===========================

S = load(targetMapFile);

q_avg = double(S.q_avg);
Xtarget_m = double(S.Xtarget_m);
Ytarget_m = double(S.Ytarget_m);
Rtarget_m = double(S.Rtarget_m);
ThetaTarget_deg = double(S.ThetaTarget_deg);
roi_mask = logical(S.roi_mask);

valid_map = roi_mask & isfinite(q_avg) & q_avg > q_threshold_MWm2;
idx_valid = find(valid_map);

if isempty(idx_valid)
    error("No valid heated target pixels found. Lower q_threshold_MWm2.");
end

idx_trace = select_target_pixels( ...
    idx_valid, q_avg, max_tubes_to_plot, tube_selection_mode);

nTrace = numel(idx_trace);

fprintf("Valid heated pixels = %d\n", numel(idx_valid));
fprintf("Selected tubes      = %d\n", nTrace);
fprintf("Target z = %.3f m, dump z = %.3f m\n", z_target, z_start);

%% ===================== LOAD B FIELD / PSI =============================

r = double(ncread(bfieldFile, "r"));
z = double(ncread(bfieldFile, "z"));
psi = double(ncread(bfieldFile, "psi"));

r = r(:);
z = z(:);

if ~isequal(size(psi), [numel(r), numel(z)])
    psi = permute(psi, [2 1]);
end

Fpsi = griddedInterpolant({r,z}, psi, "linear", "none");

Br_field = double(ncread(bfieldFile, "br"));
Bt_field = double(ncread(bfieldFile, "bt"));
Bz_field = double(ncread(bfieldFile, "bz"));

if ~isequal(size(Bz_field), [numel(r), numel(z)])
    Br_field = permute(Br_field, [2 1]);
    Bt_field = permute(Bt_field, [2 1]);
    Bz_field = permute(Bz_field, [2 1]);
end

Bmag = sqrt(Br_field.^2 + Bt_field.^2 + Bz_field.^2);
FBmag = griddedInterpolant({r, z}, Bmag, "linear", "none");

%% ===================== BUILD Z LINE ==================================

zline = z_target:dz_map:z_start;
zline = zline(:);

%% ===================== MAP TARGET PIXELS ==============================

x_path_all = cell(nTrace,1);
y_path_all = cell(nTrace,1);
z_path_all = cell(nTrace,1);
r_path_all = cell(nTrace,1);

q_trace = q_avg(idx_trace);
B_target_all = NaN(nTrace, 1);

for i = 1:nTrace

    r0 = Rtarget_m(idx_trace(i));
    th = ThetaTarget_deg(idx_trace(i));

    psi0 = Fpsi(r0, z_target);
    B_target_all(i) = FBmag(r0, z_target);

    if ~isfinite(psi0)
        continue
    end

    rline = map_psi_to_rline(psi0, r, zline, Fpsi);

    good = isfinite(rline);

    rline = rline(good);
    zuse  = zline(good);

    x_path_all{i} = rline .* cosd(th);
    y_path_all{i} = rline .* sind(th);
    z_path_all{i} = zuse;
    r_path_all{i} = rline;
end

%% ===================== COLOR SETTINGS =================================

cmap = turbo(256);
qmin = min(q_trace, [], "omitnan");
qmax = max(q_trace, [], "omitnan");

%% ===================== 3D FLUX-TUBE PLOT ==============================

figure("Color","w","Position",[60 60 1500 900])
hold on
grid on
box on

plot_vessel_3D(vessel_top_z, vessel_top_R)

plot_aperture_3D(left_aperture_z, left_aperture_R)
plot_aperture_3D(mid_aperture_z, mid_aperture_R)
plot_aperture_3D(right_aperture_z, right_aperture_R)

plot_helicon_window_3D(helicon_window_z, helicon_window_R)
plot_limiter_3D(limiter_z, limiter_R)

text(0.075, 0, mean(helicon_window_z), ...
    "Helicon Window", ...
    "Color", "b", ...
    "FontWeight", "bold", ...
    "FontSize", 12)

text(0.075, 0, mean(limiter_z), ...
    "Limiter", ...
    "Color", "g", ...
    "FontWeight", "bold", ...
    "FontSize", 12)

plot_coils_3D(X, R1, R2, wZ, ps)

plot_target_plane_3D(z_target, 0.145, [0 0 1])
plot_dump_plane_3D(z_start, 0.165, [1 0 0])

plot3([0 0], [0 0], [z_start z_end], ...
    "k--", "LineWidth", 1.0)

for i = 1:nTrace

    if isempty(x_path_all{i})
        continue
    end

    ci = color_index(q_trace(i), qmin, qmax, 256);

    plot3(x_path_all{i}, y_path_all{i}, z_path_all{i}, ...
        "Color", cmap(ci,:), ...
        "LineWidth", 0.8)
end

scatter3(Xtarget_m(idx_trace), Ytarget_m(idx_trace), ...
    z_target * ones(size(idx_trace)), ...
    18, q_trace, "filled")

xlabel("x [m]")
ylabel("y [m]")
zlabel("z [m]")
title("3D Constant-\psi Flux Tubes Carrying Target Heat Flux")
axis equal
view(-35,22)
camlight headlight
lighting gouraud
xlim([-0.18 0.18])
ylim([-0.18 0.18])
zlim([z_start z_end])
colormap(turbo)

cb = colorbar;
ylabel(cb, "q_{target} [MW/m^2]")

%% ===================== R-Z PROJECTION =================================

figure("Color","w","Position",[80 80 1550 850])
hold on
grid on
box on

colormap(turbo)
clim([qmin qmax])

plot_vessel_RZ(vessel_top_z, vessel_top_R)

plot_aperture_RZ(left_aperture_z, left_aperture_R)
plot_aperture_RZ(mid_aperture_z, mid_aperture_R)
plot_aperture_RZ(right_aperture_z, right_aperture_R)

plot_helicon_window_RZ(helicon_window_z, helicon_window_R, helicon_window_thick)
plot_limiter_RZ(limiter_z, limiter_R, limiter_thick)

plot_coils_RZ(X, R1, R2, wZ, ps)

xline(z_target, "b--", "target", ...
    "LineWidth", 1.5, ...
    "HandleVisibility","off")

xline(z_start, "r--", "dump", ...
    "LineWidth", 1.5, ...
    "HandleVisibility","off")

yline(0, "k:", "HandleVisibility","off")

for i = 1:nTrace

    if isempty(r_path_all{i})
        continue
    end

    plot(z_path_all{i}, r_path_all{i}, ...
        "Color", [0.55 0.55 0.55], ...
        "LineWidth", 0.45, ...
        "HandleVisibility","off")

    plot(z_path_all{i}, -r_path_all{i}, ...
        "Color", [0.55 0.55 0.55], ...
        "LineWidth", 0.45, ...
        "HandleVisibility","off")
end

z_sample_values = [z_target, mean(helicon_window_z), z_start];

sample_labels = ["target", ...
                 "helicon window", ...
                 "dump"];

marker_styles = ["o", "s", "d"];

for is = 1:numel(z_sample_values)

    zslice = z_sample_values(is);

    [~, ~, rs, ~, qs] = sample_flux_tubes_at_z( ...
        x_path_all, y_path_all, z_path_all, q_trace, zslice, FBmag, B_target_all);

    if isempty(rs)
        continue
    end

    scatter(zslice * ones(size(rs)), rs, ...
        30, qs, marker_styles(is), ...
        "filled", ...
        "MarkerEdgeColor", "k", ...
        "DisplayName", sample_labels(is))

    scatter(zslice * ones(size(rs)), -rs, ...
        30, qs, marker_styles(is), ...
        "filled", ...
        "MarkerEdgeColor", "k", ...
        "HandleVisibility","off")
end

xlabel("z [m]")
ylabel("Radial Distance [m]")
title("Flux Tubes and Transported Target Heat Flux")
xlim([z_start z_end])
ylim([-0.18 0.18])

cb = colorbar;
ylabel(cb, "q_{target} [MW/m^2]")

legend("Location","bestoutside")

%% ===================== TARGET HEAT MAP WITH ROI ========================

figure("Color","w","Position",[100 100 750 650])

imagesc(S.x_cm, S.y_cm, S.q_avg_plot)
set(gca,"YDir","reverse")
axis image
xlabel("target x [cm]")
ylabel("target y [cm]")
title(sprintf("Target Heat-Flux Map with ROI at z = %.2f m", z_target))

cb = colorbar;
ylabel(cb, "q_{target} [MW/m^2]")

hold on

if isfield(S, "roi_cx_cm") && isfield(S, "roi_cy_cm") && isfield(S, "roi_radius_cm")
    plot_circle_cm(S.roi_cx_cm, S.roi_cy_cm, S.roi_radius_cm, "k", 2.5);
end

clim([0 0.8])
colormap(turbo)

%% ===================== CROSS-SECTION SLICES ===========================

z_slice_values = [mean(helicon_window_z), mean(limiter_z), z_start];

z_slice_names  = ["Helicon window", ...
                  "Limiter", ...
                  "Dump"];

figure("Color","w","Position",[120 120 1500 500])

for is = 1:numel(z_slice_values)

    zslice = z_slice_values(is);

    [xs, ys, ~, ~, qs] = sample_flux_tubes_at_z( ...
        x_path_all, y_path_all, z_path_all, q_trace, zslice, FBmag, B_target_all);

    subplot(1,3,is)

    scatter(xs, ys, 28, qs, "filled")
    hold on

    if is == 1
        plot_circle(helicon_window_R, "b", 2.0)
    elseif is == 2
        plot_circle(limiter_R, "g", 2.0)
    else
        plot_circle(0.165, "r", 2.0)
    end

    axis equal
    grid on
    box on
    xlabel("x [m]")
    ylabel("y [m]")
    title(z_slice_names(is))
    xlim([-0.16 0.16])
    ylim([-0.16 0.16])
    colormap(turbo)
    clim([qmin qmax])
    colorbar
end

sgtitle("Target Heat Flux Carried by Flux Tubes at Selected Axial Planes")

%% ===================== SAVE ===========================================

save(outputFile, ...
    "x_path_all", "y_path_all", "z_path_all", "r_path_all", ...
    "q_trace", "B_target_all", "idx_trace", ...
    "z_target", "z_start", "z_end", "zline", ...
    "X", "R1", "R2", "wZ", "ps", ...
    "vessel_top_z", "vessel_top_R", ...
    "helicon_window_z", "helicon_window_R", ...
    "limiter_z", "limiter_R", ...
    "left_aperture_z", "left_aperture_R", ...
    "mid_aperture_z", "mid_aperture_R", ...
    "right_aperture_z", "right_aperture_R", ...
    "-v7.3");

fprintf("Saved %s\n", outputFile);

%% ===================== LOCAL FUNCTIONS ================================

function idx_trace = select_target_pixels(idx_valid, q_avg, maxN, mode)

    idx_valid = idx_valid(:);

    if numel(idx_valid) <= maxN
        idx_trace = idx_valid;
        return
    end

    mode = string(mode);

    switch mode
        case "top"
            [~, order] = sort(q_avg(idx_valid), "descend");
            idx_trace = idx_valid(order(1:maxN));

        case "random"
            rng(1)
            order = randperm(numel(idx_valid), maxN);
            idx_trace = idx_valid(order);

        case "weighted"
            nTop = round(0.45 * maxN);
            nRand = maxN - nTop;

            [~, orderTop] = sort(q_avg(idx_valid), "descend");
            idx_top = idx_valid(orderTop(1:nTop));

            remaining = setdiff(idx_valid, idx_top);
            rng(1)

            if numel(remaining) > nRand
                idx_rand = remaining(randperm(numel(remaining), nRand));
            else
                idx_rand = remaining;
            end

            idx_trace = [idx_top(:); idx_rand(:)];

        otherwise
            error("Unknown tube_selection_mode.");
    end
end

function rline = map_psi_to_rline(psi0, rgrid, zline, Fpsi)

    rline = NaN(size(zline));

    for k = 1:numel(zline)

        zk = zline(k);

        psi_profile = Fpsi(rgrid, zk * ones(size(rgrid)));
        good = isfinite(psi_profile);

        if nnz(good) < 2
            continue
        end

        rg = rgrid(good);
        pg = psi_profile(good);

        diffpsi = pg - psi0;

        cross_idx = find(diffpsi(1:end-1).*diffpsi(2:end) <= 0, ...
            1, "first");

        if isempty(cross_idx)
            [~, imin] = min(abs(diffpsi));
            rline(k) = rg(imin);
        else
            r1 = rg(cross_idx);
            r2 = rg(cross_idx+1);
            p1 = pg(cross_idx);
            p2 = pg(cross_idx+1);

            if abs(p2-p1) < 1e-14
                rline(k) = r1;
            else
                rline(k) = r1 + (psi0-p1) * (r2-r1) / (p2-p1);
            end
        end
    end
end

function [xs, ys, rs, ths, qs] = sample_flux_tubes_at_z( ...
    x_path_all, y_path_all, z_path_all, q_trace, zslice, FBmag, B_target_all)

    xs = [];
    ys = [];
    rs = [];
    ths = [];
    qs = [];

    for i = 1:numel(z_path_all)

        zpath = z_path_all{i};

        if isempty(zpath)
            continue
        end

        zmin = min(zpath);
        zmax = max(zpath);

        if zslice < zmin || zslice > zmax
            continue
        end

        xpath = x_path_all{i};
        ypath = y_path_all{i};

        xq = interp1(zpath, xpath, zslice, "linear");
        yq = interp1(zpath, ypath, zslice, "linear");

        if ~isfinite(xq) || ~isfinite(yq)
            continue
        end

        rq = hypot(xq, yq);
        thq = mod(atan2d(yq, xq), 360);

        xs(end+1,1) = xq;
        ys(end+1,1) = yq;
        rs(end+1,1) = rq;
        ths(end+1,1) = thq;
        B_local = FBmag(rq, zslice);
        B_ref   = B_target_all(i);
        if isfinite(B_local) && isfinite(B_ref) && B_ref > 0
            qs(end+1,1) = q_trace(i) * B_local / B_ref;
        else
            qs(end+1,1) = q_trace(i);
        end
    end
end

function ci = color_index(value, vmin, vmax, nColor)

    if ~isfinite(value) || vmax <= vmin
        ci = 1;
    else
        ci = round(1 + (nColor-1)*(value-vmin)/(vmax-vmin));
        ci = max(1, min(nColor, ci));
    end
end

function c = coil_color(ps_value)

    ps_string = string(ps_value);

    if ps_string == "PS1"
        c = [1 0 0];
    elseif ps_string == "TR1"
        c = [0 0 1];
    elseif ps_string == "TR2"
        c = [0 0.8 0];
    else
        c = [0.7 0.7 0.7];
    end
end

function label = coil_label(ps_value)

    ps_string = string(ps_value);

    if ps_string == "PS1"
        label = "5000 A";
    elseif ps_string == "TR1"
        label = "3000 A";
    elseif ps_string == "TR2"
        label = "600 A";
    else
        label = ps_string;
    end
end

function plot_vessel_RZ(zTop, rTop)

    plot(zTop, rTop, ...
        "k-", ...
        "LineWidth", 2.6, ...
        "HandleVisibility","off")

    plot(zTop, -rTop, ...
        "k-", ...
        "LineWidth", 2.6, ...
        "HandleVisibility","off")

    plot([zTop(1) zTop(1)], [-rTop(1) rTop(1)], ...
        "k-", "LineWidth", 2.6, "HandleVisibility","off")

    plot([zTop(end) zTop(end)], [-rTop(end) rTop(end)], ...
        "k-", "LineWidth", 2.6, "HandleVisibility","off")
end

function plot_vessel_3D(zTop, rTop)

    theta = linspace(0, 2*pi, 180);

    for k = 1:numel(zTop)-1

        z1 = zTop(k);
        z2 = zTop(k+1);
        R1 = rTop(k);
        R2 = rTop(k+1);

        if abs(z2-z1) > 1e-12 && abs(R2-R1) < 1e-12

            [TH, ZZ] = meshgrid(theta, linspace(z1, z2, 25));
            X = R1*cos(TH);
            Y = R1*sin(TH);

            surf(X, Y, ZZ, ...
                "FaceColor", [0.75 0.75 0.75], ...
                "FaceAlpha", 0.05, ...
                "EdgeColor", "none")

            plot3(R1*cos(theta), R1*sin(theta), z1*ones(size(theta)), ...
                "k-", "LineWidth", 1.5)

            plot3(R1*cos(theta), R1*sin(theta), z2*ones(size(theta)), ...
                "k-", "LineWidth", 1.5)

        elseif abs(z2-z1) < 1e-12

            for th = linspace(0, 2*pi, 24)
                plot3([R1 R2]*cos(th), ...
                      [R1 R2]*sin(th), ...
                      [z1 z1], ...
                      "k-", "LineWidth", 1.2)
            end
        end
    end
end

function plot_aperture_RZ(z0, Ropen)

    plot([z0 z0], [-Ropen Ropen], ...
        "k-", "LineWidth", 3.0, "HandleVisibility","off")
end

function plot_aperture_3D(z0, Ropen)

    theta = linspace(0, 2*pi, 180);

    plot3(Ropen*cos(theta), Ropen*sin(theta), ...
        z0*ones(size(theta)), ...
        "k-", "LineWidth", 3.0)

    for th = linspace(0, 2*pi, 18)
        plot3([0 Ropen*cos(th)], [0 Ropen*sin(th)], [z0 z0], ...
            "k-", "LineWidth", 0.7)
    end
end

function plot_helicon_window_RZ(zwin, Rwin, thick)

    rectangle("Position", [zwin(1), Rwin, diff(zwin), thick], ...
        "FaceColor", [0 0 1], "EdgeColor", "b", "LineWidth", 1.4, ...
        "HandleVisibility","off")

    rectangle("Position", [zwin(1), -Rwin-thick, diff(zwin), thick], ...
        "FaceColor", [0 0 1], "EdgeColor", "b", "LineWidth", 1.4, ...
        "HandleVisibility","off")
end

function plot_helicon_window_3D(zwin, Rwin)

    theta = linspace(0,2*pi,120);
    [TH,ZZ] = meshgrid(theta,zwin);

    X = Rwin*cos(TH);
    Y = Rwin*sin(TH);

    surf(X,Y,ZZ, ...
        "FaceColor",[0 0 1], ...
        "FaceAlpha",0.35, ...
        "EdgeColor","none")

    plot3(Rwin*cos(theta), ...
          Rwin*sin(theta), ...
          zwin(1)*ones(size(theta)), ...
          "b-","LineWidth",3)

    plot3(Rwin*cos(theta), ...
          Rwin*sin(theta), ...
          zwin(2)*ones(size(theta)), ...
          "b-","LineWidth",3)
end

function plot_limiter_RZ(zlim, Rlim, thick)

    rectangle("Position", [zlim(1), Rlim, diff(zlim), thick], ...
        "FaceColor", [0 0.8 0], ...
        "EdgeColor", "g", ...
        "LineWidth", 1.5, ...
        "HandleVisibility","off")

    rectangle("Position", [zlim(1), -Rlim-thick, diff(zlim), thick], ...
        "FaceColor", [0 0.8 0], ...
        "EdgeColor", "g", ...
        "LineWidth", 1.5, ...
        "HandleVisibility","off")
end

function plot_limiter_3D(zlim, Rlim)

    theta = linspace(0, 2*pi, 160);

    z1 = zlim(1);
    z2 = zlim(2);

    [TH, ZZ] = meshgrid(theta, linspace(z1, z2, 2));

    X = Rlim*cos(TH);
    Y = Rlim*sin(TH);

    surf(X, Y, ZZ, ...
        "FaceColor", [0 0.8 0], ...
        "FaceAlpha", 0.25, ...
        "EdgeColor", "none")

    plot3(Rlim*cos(theta), Rlim*sin(theta), ...
        z1*ones(size(theta)), ...
        "g-", "LineWidth", 4)

    plot3(Rlim*cos(theta), Rlim*sin(theta), ...
        z2*ones(size(theta)), ...
        "g-", "LineWidth", 4)

    for th = linspace(0, 2*pi, 24)

        plot3([Rlim 0.075]*cos(th), ...
              [Rlim 0.075]*sin(th), ...
              [z1 z1], ...
              "g-", "LineWidth", 1.0)

        plot3([Rlim 0.075]*cos(th), ...
              [Rlim 0.075]*sin(th), ...
              [z2 z2], ...
              "g-", "LineWidth", 1.0)
    end
end

function plot_target_plane_3D(z0, Rmax, colorSpec)

    th = linspace(0, 2*pi, 150);

    for rr = linspace(0, Rmax, 6)
        plot3(rr*cos(th), rr*sin(th), z0*ones(size(th)), ...
            "Color", colorSpec, "LineWidth", 1.0)
    end

    text(0, 0, z0, " target", "Color", colorSpec, "FontWeight", "bold")
end

function plot_dump_plane_3D(z0, Rmax, colorSpec)

    th = linspace(0, 2*pi, 150);

    for rr = linspace(0, Rmax, 6)
        plot3(rr*cos(th), rr*sin(th), z0*ones(size(th)), ...
            "Color", colorSpec, "LineWidth", 1.0)
    end

    text(0, 0, z0, " dump", "Color", colorSpec, "FontWeight", "bold")
end

function plot_circle(R, colorSpec, lineWidth)

    th = linspace(0, 2*pi, 300);

    plot(R*cos(th), R*sin(th), ...
        "Color", colorSpec, ...
        "LineWidth", lineWidth)
end

function plot_circle_cm(cx, cy, radius, colorSpec, lineWidth)

    th = linspace(0, 2*pi, 300);

    plot(cx + radius*cos(th), ...
         cy + radius*sin(th), ...
         "Color", colorSpec, ...
         "LineWidth", lineWidth, ...
         "LineStyle", "--");
end

function plot_coils_RZ(X, R1, R2, wZ, ps)

    for j = 1:numel(X)

        z1 = X(j) - wZ(j)/2;
        dZ = wZ(j);
        dR = R2(j) - R1(j);

        c = coil_color(ps{j});

        rectangle("Position", [z1, R1(j), dZ, dR], ...
            "FaceColor", c, "EdgeColor", c, "LineWidth", 1.0, ...
            "HandleVisibility","off")

        rectangle("Position", [z1, -R2(j), dZ, dR], ...
            "FaceColor", c, "EdgeColor", c, "LineWidth", 1.0, ...
            "HandleVisibility","off")

        text(X(j), 0.5*(R1(j)+R2(j)), coil_label(ps{j}), ...
            "Color", "w", "FontWeight", "bold", "FontSize", 8, ...
            "HorizontalAlignment", "center", "Rotation", 90)

        text(X(j), -0.5*(R1(j)+R2(j)), coil_label(ps{j}), ...
            "Color", "w", "FontWeight", "bold", "FontSize", 8, ...
            "HorizontalAlignment", "center", "Rotation", 90)
    end
end

function plot_coils_3D(X, R1, R2, wZ, ps)

    theta = linspace(0, 2*pi, 90);

    for j = 1:numel(X)

        z1 = X(j) - wZ(j)/2;
        z2 = X(j) + wZ(j)/2;

        c = coil_color(ps{j});

        plot_annular_cylinder(R1(j), R2(j), z1, z2, theta, c, 0.35)
        plot_coil_edges(R1(j), R2(j), z1, z2, theta, c)
    end
end

function plot_annular_cylinder(Rin, Rout, z1, z2, theta, faceColor, alphaVal)

    [TH, ZZ] = meshgrid(theta, [z1 z2]);

    Xo = Rout*cos(TH);
    Yo = Rout*sin(TH);
    Xi = Rin*cos(TH);
    Yi = Rin*sin(TH);

    surf(Xo, Yo, ZZ, ...
        "FaceColor", faceColor, "EdgeColor", "none", "FaceAlpha", alphaVal)

    surf(Xi, Yi, ZZ, ...
        "FaceColor", faceColor, "EdgeColor", "none", "FaceAlpha", alphaVal)

    [TH2, RR] = meshgrid(theta, [Rin Rout]);

    Xcap = RR.*cos(TH2);
    Ycap = RR.*sin(TH2);

    surf(Xcap, Ycap, z1*ones(size(Xcap)), ...
        "FaceColor", faceColor, "EdgeColor", "none", "FaceAlpha", alphaVal)

    surf(Xcap, Ycap, z2*ones(size(Xcap)), ...
        "FaceColor", faceColor, "EdgeColor", "none", "FaceAlpha", alphaVal)
end

function plot_coil_edges(Rin, Rout, z1, z2, theta, colorSpec)

    for R = [Rin Rout]
        plot3(R*cos(theta), R*sin(theta), z1*ones(size(theta)), ...
            "Color", colorSpec, "LineWidth", 0.8)

        plot3(R*cos(theta), R*sin(theta), z2*ones(size(theta)), ...
            "Color", colorSpec, "LineWidth", 0.8)
    end

    for th = linspace(0, 2*pi, 12)
        plot3([Rin Rout]*cos(th), [Rin Rout]*sin(th), [z1 z1], ...
            "Color", colorSpec, "LineWidth", 0.6)

        plot3([Rin Rout]*cos(th), [Rin Rout]*sin(th), [z2 z2], ...
            "Color", colorSpec, "LineWidth", 0.6)
    end
end