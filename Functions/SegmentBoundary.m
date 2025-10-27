function [vessel_segmented] = SegmentBoundary(vessel,ds)
% SegmentBoundary:
% -------------------------------------------------------------------------
% This function is used to divide a vacuum boundary defined by "vessel" into
% smaller segments between vertices. This is useful for evaluating the
% magnetic flux at the vacuum vessel boundary in to order to determine the
% last limiting surface and this the radial extent of the plasma.
% -------------------------------------------------------------------------
% INPUT:
% vessel: structure the defines the vertices that make up the vacuum
% boundary of the device. It has at least two fields: "z" and "r"
% ds: Each line that makes up the vacuum boundary as defined by "vessel" is
% to be divided into n segments. "ds" corresponds to the desired segment
% length.
% -------------------------------------------------------------------------
% OUTPUT:
%  vessel_segmented: structure that defines the vacuum boundary of the
%  vessel. Each point is separated by its next point by a distance of
%  approximately "ds".
% -------------------------------------------------------------------------

% Start of function:
% =========================================================================
% Calculate the length of each segment:
arcLength = sqrt(diff(vessel.r).^2 + diff(vessel.z).^2);

% =========================================================================
% Initialize the segmented boundary:
vessel_segmented.z = vessel.z(1);
vessel_segmented.r = vessel.r(1);

% =========================================================================
% Perform segmentation:
for s = 1:numel(arcLength)
    nSegments = round(arcLength(s)/ds);
    if nSegments < 3
        nSegments = 2;
    end
    % Segment the line into "nSegments" elements:
    dum.r = linspace(vessel.r(s),vessel.r(s+1),nSegments);
    dum.z = linspace(vessel.z(s),vessel.z(s+1),nSegments);
    % ---------------------------------------------------------------------
    % Append segmented line:
    vessel_segmented.z = [vessel_segmented.z,dum.z(2:end)];
    vessel_segmented.r = [vessel_segmented.r,dum.r(2:end)];
end


end

