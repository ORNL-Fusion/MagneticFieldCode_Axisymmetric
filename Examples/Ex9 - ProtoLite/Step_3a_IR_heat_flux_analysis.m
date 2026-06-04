clc
clear
close all
format long

%% ========================= USER INPUTS ================================

datadir   = "2026-03-11";
shotnum   = 1002187;

dt_frame  = 20e-3;       % Frame spacing [s]
dT_ref    = 60.9;        % Reference temperature rise [K]

% Average only over this time window
t_start = 1.1;          % s
t_end   = 1.3;          % s

frame_chk = 120;

% Initial box only defines ROI center
roi_y = 151:200;
roi_x = 351:400;

roi_radius_cm = 3.0;     % Circular ROI radius [cm]

xlim_img = [];
ylim_img = [];

clim_temp = [28 40];
clim_qavg = [];          % Example: [0 0.6]

show_roi_rectangle = false;
show_roi_circle    = true;

%% ========================= LOAD DATA =================================

data_folder = fullfile(pwd, "data", datadir, string(shotnum));
files = dir(fullfile(data_folder, "*.csv"));

if isempty(files)
    error("No CSV files found in: %s", data_folder);
end

frame_nums = nan(numel(files),1);

for i = 1:numel(files)
    token = regexp(files(i).name, "(\d+)\.csv$", "tokens", "once");
    if ~isempty(token)
        frame_nums(i) = str2double(token{1});
    end
end

valid = ~isnan(frame_nums);
files = files(valid);
frame_nums = frame_nums(valid);

if isempty(files)
    error("No files with valid frame numbers found in: %s", data_folder);
end

[frame_nums, sort_idx] = sort(frame_nums);
files = files(sort_idx);

first_img = readmatrix(fullfile(data_folder, files(1).name));

[ny, nx] = size(first_img);
nt = numel(files);

temperature = zeros(ny, nx, nt);
temperature(:,:,1) = first_img;

for k = 2:nt
    img = readmatrix(fullfile(data_folder, files(k).name));

    if ~isequal(size(img), [ny nx])
        error("Frame size mismatch in file: %s", files(k).name);
    end

    temperature(:,:,k) = img;
end

if nt < 2
    error("At least two frames are required.");
end

%% ====================== TIME DERIVATIVE ==============================

dTdt = zeros(size(temperature));

dTdt(:,:,1) = (temperature(:,:,2) - temperature(:,:,1)) / dt_frame;

for k = 2:nt-1
    dTdt(:,:,k) = (temperature(:,:,k+1) - temperature(:,:,k-1)) / ...
                  (2 * dt_frame);
end

dTdt(:,:,nt) = (temperature(:,:,nt) - temperature(:,:,nt-1)) / dt_frame;

%% ====================== MATERIAL PROPERTIES ==========================

cp = 0.125 * 4.184 / 1e-3;    % J/(kg K)
rho = 7894;                   % kg/m^3
thickness = 0.06 * 0.0254;    % m

heatflux_factor = rho * thickness * cp / 1e6;

qss = heatflux_factor * dTdt;             % MW/m^2
q0  = heatflux_factor * (dT_ref / 0.5);   % MW/m^2

fprintf("\nReference heat flux q0 = %.6f MW/m^2\n", q0);

%% ======================== TIME VECTOR ================================

t = (0:nt-1) * dt_frame;

integration_idx = find(t >= t_start & t <= t_end);

if isempty(integration_idx)
    error("No frames found in requested time window %.3f to %.3f s.", ...
        t_start, t_end);
end

dt_window = t(integration_idx(end)) - t(integration_idx(1));

idx_chk = frame_chk + 1;

%% ====================== SPATIAL COORDINATES ==========================

px = 63 / 0.0254;     % pixels/m
py = 63 / 0.0254;     % pixels/m

x_cm = ((0:nx-1) / px) * 100;
y_cm = ((0:ny-1) / py) * 100;

[Xcm, Ycm] = meshgrid(x_cm, y_cm);

fprintf("Image size          = %d x %d pixels\n", nx, ny);
fprintf("Full x extent       = %.3f to %.3f cm\n", x_cm(1), x_cm(end));
fprintf("Full y extent       = %.3f to %.3f cm\n", y_cm(1), y_cm(end));

%% ====================== CIRCULAR ROI =================================

if max(roi_y) > ny || max(roi_x) > nx
    error("ROI exceeds image dimensions.");
end

roi_x0_cm = x_cm(roi_x(1));
roi_x1_cm = x_cm(roi_x(end));
roi_y0_cm = y_cm(roi_y(1));
roi_y1_cm = y_cm(roi_y(end));

roi_center_x_cm = 0.5 * (roi_x0_cm + roi_x1_cm);
roi_center_y_cm = 0.5 * (roi_y0_cm + roi_y1_cm);

roi_mask = ((Xcm - roi_center_x_cm).^2 + ...
            (Ycm - roi_center_y_cm).^2) <= roi_radius_cm^2;

roi_area_cm2 = sum(roi_mask(:)) * mean(diff(x_cm)) * mean(diff(y_cm));

fprintf("Circular ROI center = (%.3f, %.3f) cm\n", ...
    roi_center_x_cm, roi_center_y_cm);
fprintf("Circular ROI radius = %.3f cm\n", roi_radius_cm);
fprintf("Circular ROI area   = %.3f cm^2\n", roi_area_cm2);

roi_info.x0 = roi_x0_cm;
roi_info.y0 = roi_y0_cm;
roi_info.width = roi_x1_cm - roi_x0_cm;
roi_info.height = roi_y1_cm - roi_y0_cm;
roi_info.cx = roi_center_x_cm;
roi_info.cy = roi_center_y_cm;
roi_info.radius = roi_radius_cm;
roi_info.show_rectangle = show_roi_rectangle;
roi_info.show_circle = show_roi_circle;

%% ======================== ROI ANALYSIS ===============================

roi_q = zeros(nt,1);
roi_T = zeros(nt,1);

for k = 1:nt
    q_frame = qss(:,:,k);
    T_frame = temperature(:,:,k);

    roi_q(k) = mean(q_frame(roi_mask), "omitnan");
    roi_T(k) = mean(T_frame(roi_mask), "omitnan");
end

trap_q = trapz(t(integration_idx), roi_q(integration_idx));
avg_q  = trap_q / dt_window;

avg_T_window = mean(roi_T(integration_idx), "omitnan");

fprintf("Time window          = %.3f to %.3f s\n", ...
    t(integration_idx(1)), t(integration_idx(end)));
fprintf("Integrated heat flux = %.6f MW/m^2*s\n", trap_q);
fprintf("Average heat flux    = %.6f MW/m^2\n", avg_q);
fprintf("Average temperature  = %.6f C\n", avg_T_window);
fprintf("Integration window   = %.6f s\n", dt_window);

%% ====================== AVERAGE 2D MAPS ==============================

temperature_avg_window = mean(temperature(:,:,integration_idx), 3, "omitnan");
qss_avg_window = mean(qss(:,:,integration_idx), 3, "omitnan");

qss_avg_window_plot = qss_avg_window;
qss_avg_window_plot(qss_avg_window_plot < 0) = 0;

%% ============================ PLOTS ==================================

figure("Color","w", ...
       "Position",[100 100 1300 900])

subplot(2,2,1)
plot_temperature_frame( ...
    x_cm, y_cm, temperature_avg_window, ...
    xlim_img, ylim_img, clim_temp, ...
    sprintf("Average Temperature: %.2f-%.2f s", t_start, t_end), ...
    roi_info);

subplot(2,2,2)
plot_heatflux_frame( ...
    x_cm, y_cm, qss_avg_window_plot, ...
    xlim_img, ylim_img, clim_qavg, ...
    sprintf("Average Heat Flux: %.2f-%.2f s", t_start, t_end), ...
    roi_info);

subplot(2,2,3)
plot(t, roi_T, "LineWidth", 1.7)
hold on
xline(t(integration_idx(1)), "k--", "Alpha", 0.5)
xline(t(integration_idx(end)), "k--", "Alpha", 0.5)
grid on
box on
xlim([1 1.5])
ylim padded
xlabel("t [s]")
ylabel("T [C]")
title(sprintf("ROI-Averaged Temperature, Radius = %.1f cm", roi_radius_cm))

subplot(2,2,4)
plot(t, roi_q, "LineWidth", 1.7)
hold on
xline(t(integration_idx(1)), "k--", "Alpha", 0.5)
xline(t(integration_idx(end)), "k--", "Alpha", 0.5)
yline(avg_q, "r--", "LineWidth", 1.4)
grid on
box on
xlim([1 1.5])
ylim padded
xlabel("t [s]")
ylabel("q_{ss} [MW/m^2]")
title(sprintf("ROI-Averaged Heat Flux, Average = %.4f MW/m^2", avg_q))

sgtitle(sprintf("Average IR Temperature and Heat Flux Analysis, ROI Radius = %.1f cm", ...
    roi_radius_cm))

%% ======================= OPTIONAL CHECK FRAME ========================

if idx_chk <= nt

    figure("Color","w", ...
           "Position",[150 150 1200 500])

    subplot(1,2,1)
    plot_temperature_frame( ...
        x_cm, y_cm, temperature(:,:,idx_chk), ...
        xlim_img, ylim_img, clim_temp, ...
        sprintf("Temperature Frame %d", frame_chk), roi_info);

    subplot(1,2,2)
    plot_heatflux_frame( ...
        x_cm, y_cm, qss(:,:,idx_chk), ...
        xlim_img, ylim_img, clim_qavg, ...
        sprintf("Heat Flux Frame %d", frame_chk), roi_info);

end

%% ======================== LOCAL FUNCTIONS ============================

function plot_temperature_frame( ...
    x_cm, y_cm, img, ...
    xlimits, ylimits, climits, ttl, roi_info)

    imagesc(x_cm, y_cm, img)

    if isempty(xlimits)
        xlim([x_cm(1), x_cm(end)])
    else
        xlim(xlimits)
    end

    if isempty(ylimits)
        ylim([y_cm(1), y_cm(end)])
    else
        ylim(ylimits)
    end

    set(gca, "YDir", "reverse")

    if ~isempty(climits)
        clim(climits)
    end

    xlabel("x [cm]")
    ylabel("y [cm]")
    title(ttl)

    cb = colorbar;
    ylabel(cb, "Temperature [C]")

    hold on
    draw_roi(roi_info)

    axis tight
end

function plot_heatflux_frame( ...
    x_cm, y_cm, img, ...
    xlimits, ylimits, climits, ttl, roi_info)

    imagesc(x_cm, y_cm, img)

    if isempty(xlimits)
        xlim([x_cm(1), x_cm(end)])
    else
        xlim(xlimits)
    end

    if isempty(ylimits)
        ylim([y_cm(1), y_cm(end)])
    else
        ylim(ylimits)
    end

    set(gca, "YDir", "reverse")

    if ~isempty(climits)
        clim(climits)
    end

    xlabel("x [cm]")
    ylabel("y [cm]")
    title(ttl)

    cb = colorbar;
    ylabel(cb, "q_{ss} [MW/m^2]")

    hold on
    draw_roi(roi_info)

    axis tight
end

function draw_roi(roi_info)

    if roi_info.show_rectangle
        rectangle( ...
            "Position", [roi_info.x0, roi_info.y0, ...
                         roi_info.width, roi_info.height], ...
            "EdgeColor", "w", ...
            "LineWidth", 2.2, ...
            "LineStyle", "-");
    end

    if roi_info.show_circle
        rectangle( ...
            "Position", [roi_info.cx - roi_info.radius, ...
                         roi_info.cy - roi_info.radius, ...
                         2 * roi_info.radius, ...
                         2 * roi_info.radius], ...
            "Curvature", [1 1], ...
            "EdgeColor", "k", ...
            "LineWidth", 2.0, ...
            "LineStyle", "--");
    end

    text( ...
        roi_info.cx, roi_info.cy, "ROI", ...
        "Color", "w", ...
        "FontSize", 11, ...
        "FontWeight", "bold", ...
        "HorizontalAlignment", "center", ...
        "VerticalAlignment", "middle", ...
        "BackgroundColor", "k", ...
        "Margin", 2);

end