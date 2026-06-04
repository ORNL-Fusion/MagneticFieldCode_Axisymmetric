%% Generate_IR_TargetHeatMap.m
% Generates: IR_TargetHeatMap.mat

clc
clear
close all
format long

%% ========================= USER INPUTS ================================

datadir = "2026-03-11";
shotnum = 1002187;

dt_frame = 20e-3;        % [s]

t_start = 1.1;          % [s]
t_end   = 1.27;          % [s]

roi_y = 151:200;
roi_x = 351:400;
roi_radius_cm = 3.0;     % [cm]

outputFile = "IR_TargetHeatMap.mat";

%% ========================= LOAD IR CSV DATA ===========================

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

%% ====================== TIME DERIVATIVE ===============================

dTdt = zeros(size(temperature));

dTdt(:,:,1) = (temperature(:,:,2) - temperature(:,:,1)) / dt_frame;

for k = 2:nt-1
    dTdt(:,:,k) = (temperature(:,:,k+1) - temperature(:,:,k-1)) / ...
                  (2 * dt_frame);
end

dTdt(:,:,nt) = (temperature(:,:,nt) - temperature(:,:,nt-1)) / dt_frame;

%% ====================== HEAT FLUX ====================================

cp = 0.125 * 4.184 / 1e-3;    % [J/kg/K]
rho = 7894;                   % [kg/m^3]
thickness = 0.06 * 0.0254;    % [m]

heatflux_factor = rho * thickness * cp / 1e6;

qss = heatflux_factor * dTdt; % [MW/m^2]

%% ====================== TIME WINDOW ==================================

t = (0:nt-1) * dt_frame;

integration_idx = find(t >= t_start & t <= t_end);

if isempty(integration_idx)
    error("No frames found between %.3f and %.3f s.", t_start, t_end);
end

%% ====================== SPATIAL COORDINATES ===========================

px = 63 / 0.0254;     % [pixels/m]
py = 63 / 0.0254;     % [pixels/m]

x_cm = ((0:nx-1) / px) * 100;
y_cm = ((0:ny-1) / py) * 100;

[Xcm, Ycm] = meshgrid(x_cm, y_cm);

%% ====================== ROI ==========================================

roi_cx_cm = 0.5 * (x_cm(roi_x(1)) + x_cm(roi_x(end)));
roi_cy_cm = 0.5 * (y_cm(roi_y(1)) + y_cm(roi_y(end)));

roi_mask = ((Xcm - roi_cx_cm).^2 + ...
            (Ycm - roi_cy_cm).^2) <= roi_radius_cm^2;

%% ====================== AVERAGE MAPS =================================

T_avg = mean(temperature(:,:,integration_idx), 3, "omitnan");
q_avg = mean(qss(:,:,integration_idx), 3, "omitnan");

q_avg_plot = q_avg;
q_avg_plot(q_avg_plot < 0) = 0;

%% ====================== TARGET COORDINATES ============================

Xtarget_m = (Xcm - roi_cx_cm) * 1e-2;
Ytarget_m = (Ycm - roi_cy_cm) * 1e-2;

Rtarget_m = hypot(Xtarget_m, Ytarget_m);
ThetaTarget_deg = mod(atan2d(Ytarget_m, Xtarget_m), 360);

%% ====================== ROI AVERAGES =================================

roi_T = zeros(nt,1);
roi_q = zeros(nt,1);

for k = 1:nt
    T_frame = temperature(:,:,k);
    q_frame = qss(:,:,k);

    roi_T(k) = mean(T_frame(roi_mask), "omitnan");
    roi_q(k) = mean(q_frame(roi_mask), "omitnan");
end

avg_T = mean(roi_T(integration_idx), "omitnan");
avg_q = mean(roi_q(integration_idx), "omitnan");

fprintf("\nSaved target heat-map variables:\n");
fprintf("Average temperature = %.6f C\n", avg_T);
fprintf("Average heat flux   = %.6f MW/m^2\n", avg_q);
fprintf("ROI radius          = %.3f cm\n", roi_radius_cm);
fprintf("Time window         = %.3f to %.3f s\n", t_start, t_end);

%% ====================== SAVE OUTPUT ==================================

save(outputFile, ...
    "T_avg", "q_avg", "q_avg_plot", ...
    "temperature", "qss", ...
    "x_cm", "y_cm", ...
    "Xtarget_m", "Ytarget_m", ...
    "Rtarget_m", "ThetaTarget_deg", ...
    "roi_mask", "roi_radius_cm", ...
    "roi_cx_cm", "roi_cy_cm", ...
    "roi_T", "roi_q", ...
    "avg_T", "avg_q", ...
    "t", "t_start", "t_end", "integration_idx", ...
    "heatflux_factor", ...
    "-v7.3");

fprintf("Output saved: %s\n", outputFile);

%% ====================== QUICK CHECK PLOT ==============================

figure("Color","w","Position",[100 100 1200 500])

subplot(1,2,1)
imagesc(x_cm, y_cm, T_avg)
set(gca,"YDir","reverse")
axis image
xlabel("x [cm]")
ylabel("y [cm]")
title("Average Temperature")
colorbar
hold on
plot_circle(roi_cx_cm, roi_cy_cm, roi_radius_cm)

subplot(1,2,2)
imagesc(x_cm, y_cm, q_avg_plot)
set(gca,"YDir","reverse")
axis image
xlabel("x [cm]")
ylabel("y [cm]")
title("Average Heat Flux")
cb = colorbar;
ylabel(cb, "q [MW/m^2]")
hold on
plot_circle(roi_cx_cm, roi_cy_cm, roi_radius_cm)

%% ====================== LOCAL FUNCTION ================================

function plot_circle(cx, cy, radius)

    th = linspace(0, 2*pi, 300);
    plot(cx + radius*cos(th), ...
         cy + radius*sin(th), ...
         "k--", "LineWidth", 2.0);

end