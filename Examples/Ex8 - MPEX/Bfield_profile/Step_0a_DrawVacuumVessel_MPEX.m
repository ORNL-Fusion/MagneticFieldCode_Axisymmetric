function [vessel, window, target, vessel_raw] = Step_0a_DrawVacuumVessel_MPEX(boundaryFile, doPlot)
% Step_0a_DrawVacuumVessel_MPEX
% Read MPEX vacuum-vessel coordinates carefully from MPEX_innerbound.xlsx.
%
% Outputs
%   vessel.z, vessel.r          : upper inner-wall branch only [m]
%   vessel_raw.full_z/full_r    : full closed boundary from spreadsheet [m]
%   vessel_raw.upper_z/upper_r  : upper branch [m]
%   vessel_raw.lower_z/lower_r  : lower branch [m]

if nargin < 1 || isempty(boundaryFile)
    boundaryFile = 'MPEX_innerbound.xlsx';
end
if nargin < 2
    doPlot = true;
end

A = readcell(boundaryFile, 'Sheet', 1);

% Main boundary: columns B:E, starting below unit rows
mainBlock = A(4:end, 2:5);

z_cm = nan(size(mainBlock,1),1);
r_cm = nan(size(mainBlock,1),1);
for ii = 1:size(mainBlock,1)
    z_cm(ii) = localToDouble(mainBlock{ii,2});
    r_cm(ii) = localToDouble(mainBlock{ii,4});
end

mask = isfinite(z_cm) & isfinite(r_cm);
z_full = z_cm(mask) / 100;
r_full = r_cm(mask) / 100;

% Upper branch = first contiguous nonnegative-r branch
idxFirstNeg = find(r_full < 0, 1, 'first');
if isempty(idxFirstNeg)
    iUpperEnd = numel(r_full);
else
    iUpperEnd = idxFirstNeg - 1;
end

z_upper = z_full(1:iUpperEnd);
r_upper = r_full(1:iUpperEnd);
if numel(z_upper) >= 2
    if abs(z_upper(end) - z_upper(1)) < 1e-10 && abs(r_upper(end) - r_upper(1)) < 1e-10
        z_upper(end) = [];
        r_upper(end) = [];
    end
end

if isempty(idxFirstNeg)
    z_lower = [];
    r_lower = [];
else
    z_lower = z_full(idxFirstNeg:end);
    r_lower = r_full(idxFirstNeg:end);
    if numel(z_lower) >= 1 && abs(z_lower(end) - z_upper(1)) < 1e-10 && abs(r_lower(end) - r_upper(1)) < 1e-10
        z_lower(end) = [];
        r_lower(end) = [];
    end
end

vessel.z = z_upper(:)';
vessel.r = r_upper(:)';

vessel_raw.full_z  = z_full(:)';
vessel_raw.full_r  = r_full(:)';
vessel_raw.upper_z = z_upper(:)';
vessel_raw.upper_r = r_upper(:)';
vessel_raw.lower_z = z_lower(:)';
vessel_raw.lower_r = r_lower(:)';

% Window/source reference: force helicon window center at z = 0
window.z = 0.0;

if window.z < min(vessel.z) || window.z > max(vessel.z)
    error('window.z = 0 lies outside the vessel z-range.');
end

% Stepped vessel geometry can contain repeated z values, so avoid interp1.
tol = 1e-8;
idxExact = find(abs(vessel.z - window.z) < tol);

if ~isempty(idxExact)
    % If z=0 exists multiple times, use the smallest positive radius there
    window.r = min(vessel.r(idxExact));
else
    % Find nearest points on each side of z=0
    idxL = find(vessel.z < window.z, 1, 'last');
    idxR = find(vessel.z > window.z, 1, 'first');

    if isempty(idxL) || isempty(idxR)
        [~, idxNear] = min(abs(vessel.z - window.z));
        window.r = vessel.r(idxNear);
    else
        % Conservative choice for stepped wall: use smaller of the two radii
        window.r = min(vessel.r([idxL idxR]));
    end
end

window.label = 'helicon window center';

% Target plane from dedicated target columns, if present.
targetBlock = A(4:end, 7:8);
targetZ = [];
targetR = [];
for ii = 1:size(targetBlock,1)
    zt = localToDouble(targetBlock{ii,1});
    rt = localToDouble(targetBlock{ii,2});
    if isfinite(zt) && isfinite(rt)
        targetZ(end+1,1) = zt / 100;
        targetR(end+1,1) = rt / 100;
    end
end

if ~isempty(targetZ)
    target.z = median(targetZ);
    target.r = max(abs(targetR));
else
    target.z = vessel.z(end);
    target.r = 0.025;
end
target.label = 'downstream target plane';

if doPlot
    figure('Color', 'w'); hold on
    plot(vessel_raw.full_z, vessel_raw.full_r, 'r.-', 'DisplayName', 'raw full boundary')
    plot(vessel.z,  vessel.r, 'k-', 'LineWidth', 1.5, 'DisplayName', 'upper branch used')
    plot(vessel.z, -vessel.r, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off')
    plot(window.z,  window.r, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'window ref')
    plot(window.z, -window.r, 'bo', 'MarkerFaceColor', 'b', 'HandleVisibility', 'off')
    plot([target.z target.z], [-target.r target.r], 'm-', 'LineWidth', 2, 'DisplayName', 'target')
    axis equal
    grid on
    box on
    xlabel('z [m]')
    ylabel('r [m]')
    title('MPEX vacuum-vessel inner boundary (careful read)')
    legend('Location', 'best')
end
end

function x = localToDouble(v)
if isnumeric(v)
    x = double(v);
elseif isstring(v) || ischar(v)
    x = str2double(v);
else
    x = NaN;
end
end