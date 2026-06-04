clc
clear
close all
format long

%% ========================= USER INPUTS ================================

datadir   = "2026-03-11";
shotnum   = 1002187;

dt_frame  = 20e-3;       % Frame spacing [s]
dT_ref    = 60.9;        % Reference temperature rise [K]

frame1    = 54;
frame2    = 70;
frame_chk = 120;

roi_y = 151:200;
roi_x = 351:400;

% Use empty arrays [] for full spatial extent.
% Example crop:
% xlim_img = [8 14];
% ylim_img = [9 15];
xlim_img = [];
ylim_img = [];

clim_temp = [28 40];

show_roi_rectangle = true;
show_roi_circle    = true;

%% ========================= LOAD DATA =================================

data_folder = fullfile(pwd, "data", datadir, string(shotnum));

files = dir(fullfile(data_folder, "*.csv"));

if isempty(files)
    error("No CSV files found in: %s", data_folder);
end

frame_nums = nan(numel(files),1);

for i = 1:numel(files)

    token = regexp(files(i).name, "(\d+)\.csv$", ...
        "tokens", "once");

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

% temperature is stored as y x x x time.
% Therefore dT/dt must be computed along dimension 3.
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

idx1 = frame1 + 1;
idx2 = frame2 + 1;
idx_chk = frame_chk + 1;

if idx2 > nt
    error("Requested frame exceeds available frames.");
end

dt_window = t(idx2) - t(idx1);

%% ======================== ROI ANALYSIS ===============================

if max(roi_y) > ny || max(roi_x) > nx
    error("ROI exceeds image dimensions.");
end

roi_q = squeeze(mean(qss(roi_y, roi_x, :), [1 2]));

% Integration uses frames frame1 through frame2-1.
% This is equivalent to a half-open interval [frame1, frame2).
integration_idx = (frame1 + 1):frame2;

trap_q = trapz(t(integration_idx), roi_q(integration_idx));
avg_q  = trap_q / dt_window;

fprintf("Integrated heat flux = %.6f MW/m^2*s\n", trap_q);
fprintf("Average heat flux    = %.6f MW/m^2\n", avg_q);
fprintf("Integration window   = %.6f s\n", dt_window);

%% ====================== SPATIAL COORDINATES ==========================

px = 63 / 0.0254;     % pixels/m
py = 63 / 0.0254;     % pixels/m

x_cm = ((0:nx-1) / px) * 100;
y_cm = ((0:ny-1) / py) * 100;

fprintf("Image size          = %d x %d pixels\n", nx, ny);
fprintf("Full x extent       = %.3f to %.3f cm\n", x_cm(1), x_cm(end));
fprintf("Full y extent       = %.3f to %.3f cm\n", y_cm(1), y_cm(end));

%% ====================== ROI COORDINATES ==============================

roi_x0_cm = x_cm(roi_x(1));
roi_x1_cm = x_cm(roi_x(end));
roi_y0_cm = y_cm(roi_y(1));
roi_y1_cm = y_cm(roi_y(end));

roi_width_cm  = roi_x1_cm - roi_x0_cm;
roi_height_cm = roi_y1_cm - roi_y0_cm;

roi_center_x_cm = 0.5 * (roi_x0_cm + roi_x1_cm);
roi_center_y_cm = 0.5 * (roi_y0_cm + roi_y1_cm);

roi_radius_cm = 0.5 * max(abs(roi_width_cm), abs(roi_height_cm));

fprintf("ROI x extent        = %.3f to %.3f cm\n", roi_x0_cm, roi_x1_cm);
fprintf("ROI y extent        = %.3f to %.3f cm\n", roi_y0_cm, roi_y1_cm);
fprintf("ROI center          = (%.3f, %.3f) cm\n", ...
    roi_center_x_cm, roi_center_y_cm);

roi_info.x0 = roi_x0_cm;
roi_info.y0 = roi_y0_cm;
roi_info.width = roi_width_cm;
roi_info.height = roi_height_cm;
roi_info.cx = roi_center_x_cm;
roi_info.cy = roi_center_y_cm;
roi_info.radius = roi_radius_cm;
roi_info.show_rectangle = show_roi_rectangle;
roi_info.show_circle = show_roi_circle;

%% ============================ PLOTS ==================================

figure("Color","w", ...
       "Position",[100 100 1200 1000])

subplot(3,1,1)

plot_temperature_frame( ...
    x_cm, y_cm, temperature(:,:,idx1), ...
    xlim_img, ylim_img, clim_temp, ...
    sprintf("Temperature Frame %d", frame1), roi_info);

subplot(3,1,2)

plot_temperature_frame( ...
    x_cm, y_cm, temperature(:,:,idx2), ...
    xlim_img, ylim_img, clim_temp, ...
    sprintf("Temperature Frame %d", frame2), roi_info);

subplot(3,1,3)

plot(t, roi_q, "LineWidth", 1.7)

hold on

xline(t(idx1), "k--", "Alpha", 0.5)
xline(t(idx2), "k--", "Alpha", 0.5)

grid on
box on

xlim([1 1.5])
ylim padded

xlabel("t [s]")
ylabel("q_{ss} [MW/m^2]")

title("ROI-Averaged Heat Flux")

sgtitle("IR Heat Flux Analysis - Full Spatial Extent with ROI")

%% ======================= OPTIONAL FRAME ==============================

if idx_chk <= nt

    figure("Color","w")

    plot_temperature_frame( ...
        x_cm, y_cm, temperature(:,:,idx_chk), ...
        xlim_img, ylim_img, clim_temp, ...
        sprintf("Temperature Frame %d", frame_chk), roi_info);

end

%% ======================== LOCAL FUNCTION =============================

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

    clim(climits)

    xlabel("x [cm]")
    ylabel("y [cm]")
    title(ttl)

    cb = colorbar;
    ylabel(cb, "Temperature [C]")

    hold on

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

    axis tight

end