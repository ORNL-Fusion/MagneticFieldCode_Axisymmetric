function [vessel_new] = AddComponent(vessel,component)
% AddComponent:
% -------------------------------------------------------------------------
% This function is used to insert "component.z" and "component.r" into the
% structure "vessel". This is useful when one defines the general shape of
% the vacuum vessel and one needs to add components such as limiters, coil
% indentation, ports, etc.
% -------------------------------------------------------------------------
% INPUT:
% -------------------------------------------------------------------------
% vessel: structure that contains at least the following two fields:
%     .z: 
%     .r:
% Together, these two fields define the outline of the vacuum vessel in a
% 2D axisymmetric geometry such as in Proto-MPEX
% -------------------------------------------------------------------------
% OUTPUT:
% -------------------------------------------------------------------------
% vessel_new: structure that defines the updated vessel profile. It was two
% fields: "z" and "r"
% -------------------------------------------------------------------------

% Start of function:
% =========================================================================
% Append the "component" structure
vessel_new.z = [vessel.z,component.z];
vessel_new.r = [vessel.r,component.r];

% =========================================================================
% Arrange the elements in order:
[~,b] = sort(vessel_new.z);

% =========================================================================
% Compute the final version of "vessel_new"
vessel_new.z = vessel_new.z(b);
vessel_new.r = vessel_new.r(b);

end

