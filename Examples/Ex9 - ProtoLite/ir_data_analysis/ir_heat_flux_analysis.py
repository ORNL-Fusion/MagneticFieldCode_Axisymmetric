import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ========================= USER INPUTS =========================

datadir = "2026-03-11"
shotnum = 1002187

dt_frame = 20e-3      # Frame spacing [s]
dT_ref = 60.9         # Reference temperature rise [K]

frame1 = 54
frame2 = 70
frame_chk = 120

roi_y = slice(150, 200)
roi_x = slice(350, 400)

# Use full spatial extent. Set to tuples such as (8, 14) or (9, 15) to crop.
xlim_img = None
ylim_img = None

clim_temp = (28, 40)


# ========================= LOAD DATA ===========================

data_folder = Path.cwd() / "data" / datadir / str(shotnum)

files = sorted(data_folder.glob("*.csv"))

if not files:
    raise FileNotFoundError(f"No CSV files found in: {data_folder}")

frame_nums = []
valid_files = []

for file in files:
    match = re.search(r"(\d+)\.csv$", file.name)

    if match:
        frame_nums.append(int(match.group(1)))
        valid_files.append(file)

if not valid_files:
    raise FileNotFoundError(f"No files with valid frame numbers found in: {data_folder}")

sort_idx = np.argsort(frame_nums)

files = [valid_files[i] for i in sort_idx]
frame_nums = np.array(frame_nums)[sort_idx]

first_img = np.loadtxt(files[0], delimiter=",")

ny, nx = first_img.shape
nt = len(files)

temperature = np.zeros((nt, ny, nx))
temperature[0, :, :] = first_img

for k, file in enumerate(files[1:], start=1):
    img = np.loadtxt(file, delimiter=",")

    if img.shape != (ny, nx):
        raise ValueError(f"Frame size mismatch in file: {file.name}")

    temperature[k, :, :] = img

if nt < 2:
    raise ValueError("At least two frames are required.")


# ====================== TIME DERIVATIVE ========================

dTdt = np.gradient(temperature, dt_frame, axis=0)


# ====================== MATERIAL PROPERTIES ====================

cp = 0.125 * 4.184 / 1e-3     # J/(kg K)
rho = 7894                    # kg/m^3
thickness = 0.06 * 0.0254     # m

heatflux_factor = rho * thickness * cp / 1e6

qss = heatflux_factor * dTdt
q0 = heatflux_factor * (dT_ref / 0.5)

print(f"\nReference heat flux q0 = {q0:.6f} MW/m^2")


# ======================== TIME VECTOR ==========================

t = np.arange(nt) * dt_frame

idx1 = frame1
idx2 = frame2
idx_chk = frame_chk

if idx2 >= nt:
    raise IndexError(f"Requested frame2={frame2} exceeds available frames. nt={nt}")

dt_window = t[idx2] - t[idx1]


# ======================== ROI ANALYSIS =========================

if roi_y.stop > ny or roi_x.stop > nx:
    raise IndexError(
        f"ROI exceeds image dimensions. Image size = {ny} x {nx}, "
        f"ROI y stop = {roi_y.stop}, ROI x stop = {roi_x.stop}"
    )

roi_q = np.mean(qss[:, roi_y, roi_x], axis=(1, 2))

integration_idx = slice(frame1, frame2)

trap_q = np.trapz(roi_q[integration_idx], t[integration_idx])
avg_q = trap_q / dt_window

print(f"Integrated heat flux = {trap_q:.6f} MW/m^2*s")
print(f"Average heat flux    = {avg_q:.6f} MW/m^2")
print(f"Integration window   = {dt_window:.6f} s")


# ====================== SPATIAL COORDINATES ====================

px = 63 / 0.0254
py = 63 / 0.0254

x_cm = (np.arange(nx) / px) * 100
y_cm = (np.arange(ny) / py) * 100

extent = [x_cm[0], x_cm[-1], y_cm[-1], y_cm[0]]

print(f"Image size          = {nx} x {ny} pixels")
print(f"Full x extent       = {x_cm[0]:.3f} to {x_cm[-1]:.3f} cm")
print(f"Full y extent       = {y_cm[0]:.3f} to {y_cm[-1]:.3f} cm")


# ======================== PLOT FUNCTION ========================

def plot_temperature_frame(ax, img, title):
    im = ax.imshow(
        img,
        extent=extent,
        aspect="equal",
        vmin=clim_temp[0],
        vmax=clim_temp[1],
    )

    if xlim_img is not None:
        ax.set_xlim(xlim_img)
    else:
        ax.set_xlim(x_cm[0], x_cm[-1])

    if ylim_img is not None:
        ax.set_ylim(ylim_img[::-1])
    else:
        ax.set_ylim(y_cm[-1], y_cm[0])

    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_title(title)

    cb = plt.colorbar(im, ax=ax)
    cb.set_label("Temperature [C]")

    return im


# ============================ PLOTS ============================

fig, axes = plt.subplots(3, 1, figsize=(12, 12))

plot_temperature_frame(
    axes[0],
    temperature[idx1, :, :],
    f"Temperature Frame {frame1}",
)

plot_temperature_frame(
    axes[1],
    temperature[idx2, :, :],
    f"Temperature Frame {frame2}",
)

axes[2].plot(t, roi_q, linewidth=1.7)

axes[2].axvline(t[idx1], color="k", linestyle="--", alpha=0.5)
axes[2].axvline(t[idx2], color="k", linestyle="--", alpha=0.5)

axes[2].grid(True)
axes[2].set_xlim(1, 1.5)
axes[2].set_xlabel("t [s]")
axes[2].set_ylabel(r"$q_{ss}$ [MW/m$^2$]")
axes[2].set_title("ROI-Averaged Heat Flux")

fig.suptitle("IR Heat Flux Analysis - Full Spatial Extent")
fig.tight_layout()

plt.show()


# ======================= OPTIONAL FRAME ========================

if idx_chk < nt:
    fig2, ax2 = plt.subplots(figsize=(8, 6))

    plot_temperature_frame(
        ax2,
        temperature[idx_chk, :, :],
        f"Temperature Frame {frame_chk}",
    )

    fig2.tight_layout()
    plt.show()
