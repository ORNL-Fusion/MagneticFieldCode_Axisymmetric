clear; clc; close all;

%% ---------------- Constants ----------------
e    = 1.60217662e-19;   % [C]
me   = 9.10938356e-31;   % [kg]
eps0 = 8.854187e-12;     % [F/m]
kB   = 1.380649e-23;     % [J/K]

lnLambda = 15;

%% ---------------- Plasma inputs ----------------
ne = 1e20;       % electron density [m^-3]
Te_eV = 10;      % electron temperature [eV]

%% ---------------- Neutral inputs ----------------
nn = 1e19;       % neutral density [m^-3]

%% ---------------- Convert temperature ----------------
Te_J = Te_eV * e;

%% ---------------- Electron thermal velocity ----------------
vth_e = sqrt(2 * Te_J / me);

%% ---------------- Electron-ion collisions (FULL SI) ----------------
nu_ei = (ne * e^4 * lnLambda) / ...
        (3 * (2*pi)^(3/2) * eps0^2 * sqrt(me) * (Te_J)^(3/2));

lambda_ei = vth_e / nu_ei;

%% ---------------- Electron-neutral (from <σv>) ----------------
% Fit based on your figure (D2)
rate_en = 1.5e-13 * (1 - exp(-Te_eV/1.5));   % [m^3/s]

nu_en = nn * rate_en;   % <-- NO velocity multiplication
lambda_en = vth_e / nu_en;

%% ---------------- Ionization (optional, from green curve) ----------------
rate_ion = 1e-15 * Te_eV;   % crude linear fit for 5–15 eV range
nu_ion = nn * rate_ion;

%% ---------------- Total collision frequency ----------------
nu_total = nu_ei + nu_en + nu_ion;

lambda_total = vth_e / nu_total;

%% ---------------- Output ----------------
fprintf('--- Electron-Ion ---\n');
fprintf('nu_ei     = %.3e s^-1\n', nu_ei);
fprintf('lambda_ei = %.3e m\n', lambda_ei);

fprintf('\n--- Electron-Neutral (D2, from <σv>) ---\n');
fprintf('<sigma v> = %.3e m^3/s\n', rate_en);
fprintf('nu_en     = %.3e s^-1\n', nu_en);
fprintf('lambda_en = %.3e m\n', lambda_en);

fprintf('\n--- Ionization ---\n');
fprintf('nu_ion    = %.3e s^-1\n', nu_ion);

fprintf('\n--- Total ---\n');
fprintf('nu_total  = %.3e s^-1\n', nu_total);
fprintf('lambda_total = %.3e m\n', lambda_total);

%% ---------------- Regime comparison ----------------
fprintf('\n--- Dominant mechanism ---\n');

if nu_en > nu_ei
    fprintf('Electron-neutral collisions dominate (edge / gas puff region)\n');
else
    fprintf('Electron-ion collisions dominate (core plasma)\n');
end