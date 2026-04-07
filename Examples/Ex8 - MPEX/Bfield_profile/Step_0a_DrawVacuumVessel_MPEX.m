function [vessel, window, target, vessel_raw] = Step_0a_DrawVacuumVessel_MPEX(boundaryFile, doPlot)
% DrawVacuumVessel_MPEX_careful
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

% Read without assuming clean headers.
A = readcell(boundaryFile, 'Sheet', 1);

% Main boundary lives in columns B:E, starting below the unit rows.
%   col 2 = Z [in], col 3 = Z [cm], col 4 = r [in], col 5 = r [cm]
mainBlock = A(4:end, 2:5);

z_cm = nan(size(mainBlock,1),1);
r_cm = nan(size(mainBlock,1),1);
for ii = 1:size(mainBlock,1)
    z_cm(ii) = localToDouble(mainBlock{ii,2});
    r_cm(ii) = localToDouble(mainBlock{ii,4});
end

mask = isfinite(z_cm) & isfinite(r_cm);
z_full = z_cm(mask) / 100;  % cm -> m
r_full = r_cm(mask) / 100;


% axisymmetric field calculation.
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

% Lower branch from the first negative-r point onward, excluding any final
% repeated closure point if present.
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

% Window reference: use the narrow source opening near z ~ 0.02-0.06 m on
% the upper branch. Fall back to the smallest radius in the near-source zone.
sourceMask = vessel.z > 0.015 & vessel.z < 0.060;
if any(sourceMask)
    z_source = vessel.z(sourceMask);
    r_source = vessel.r(sourceMask);
    [window.r, idxMin] = min(r_source);
    window.z = z_source(idxMin);
else
    nearSource = vessel.z > -0.10 & vessel.z < 0.30;
    if any(nearSource)
        z_source = vessel.z(nearSource);
        r_source = vessel.r(nearSource);
        [window.r, idxMin] = min(r_source);
        window.z = z_source(idxMin);
    else
        window.z = 0.0224;
        window.r = 0.0622;
    end
end
window.label = 'source / window reference';

% Target plane: use the dedicated target columns if present.
% In the file, the values appear in cm even though the note says [mm].
% Example: 782.7 corresponds to the 7.827 m downstream end, so divide by 100.
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
    plot(vessel.z,  vessel.r, 'k-', 'LineWidth', 1.5)
    plot(vessel.z, -vessel.r, 'k-', 'LineWidth', 1.5)
    plot(window.z,  window.r, 'bo', 'MarkerFaceColor', 'b')
    plot(window.z, -window.r, 'bo', 'MarkerFaceColor', 'b')
    line([target.z target.z], [-target.r target.r], 'Color', 'r', 'LineWidth', 2)
    axis equal
    grid on
    box on
    xlabel('z [m]')
    ylabel('r [m]')
    title('MPEX vacuum-vessel inner boundary (careful read)')
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
