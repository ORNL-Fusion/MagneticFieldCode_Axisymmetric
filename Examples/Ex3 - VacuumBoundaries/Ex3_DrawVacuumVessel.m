% This script demonstrates how to draw a vacuum vessel, add components,
% perform operations on the boundary, perform reflections and addition
% using the functions:
% "AddComponent"
% "SegmentBoundary"

clearvars 
close all
clc

% =========================================================================
% Define the search path:
homeFolder = cd;
cd ..\..
currentFolder = cd;
functionFolder = [currentFolder,'\Functions'];
p = addpath(genpath(functionFolder));
cd(homeFolder)

% =========================================================================
% Start by defining som arbitrary vacuum vessel:
z0 = 0.5;
zE = 2.8;
d0 = 1.5;
vessel.z = [0 ,z0,2.3,zE];
vessel.r = [2 ,d0,d0 ,0 ];

figure; 
hold on
plot(vessel.z,vessel.r,'k-')
xlim([0,3])
ylim([0,2.5])
box on 
grid on

% =========================================================================
% Perform an operation on the vessel boundary:
d1 = 1.3;
cut.z = [1 ,1 ,2 ,2 ];
cut.r = [d0,d1,d1,d0];

% =========================================================================
% Add element/ perform operation:
vessel = AddComponent(vessel,cut);

figure; 
hold on
plot(vessel.z,vessel.r,'k-')
xlim([0,3])
ylim([0,2.5])
box on 
grid on

% =========================================================================
% Define skimmer geometry
skimmer.z1 = 1.3;
skimmer.w  = 0.1;
skimmer.r1 = 0.3;
skimmer.r2 = d1;
% -------------------------------------------------------------------------
% Derived quantities:
skimmer.z2 = skimmer.z1 + skimmer.w;
% -------------------------------------------------------------------------
% Skimmer profile:
skimmer.z = [skimmer.z1, skimmer.z1, skimmer.z2, skimmer.z2];
skimmer.r = [skimmer.r2, skimmer.r1, skimmer.r1, skimmer.r2];

% =========================================================================
% Add skimmer to vessel:
vessel = AddComponent(vessel,skimmer);

% =========================================================================
% Define second skimmer:
skimmer.z1 = 1.7;
skimmer.w  = 0.15;
skimmer.r1 = 0.9;
skimmer.r2 = d1;
% -------------------------------------------------------------------------
% Derived quantities:
skimmer.z2 = skimmer.z1 + skimmer.w;
% -------------------------------------------------------------------------
% Skimmer profile:
skimmer.z = [skimmer.z1, skimmer.z1, skimmer.z2, skimmer.z2];
skimmer.r = [skimmer.r2, skimmer.r1, skimmer.r1, skimmer.r2];

% =========================================================================
% Add skimmer to vessel:
vessel = AddComponent(vessel,skimmer);

figure; 
hold on
plot(vessel.z,vessel.r,'k-')
xlim([0,3])
ylim([0,2.5])
box on 
grid on

% =========================================================================
% Mirror image of the vessel:
mirror_vessel.z = -fliplr(vessel.z);
mirror_vessel.r = +fliplr(vessel.r);

% =========================================================================
% Combine vessel and its mirror image:
vessel = AddComponent(vessel,mirror_vessel);

% =========================================================================
% Segment the boundary:
% This may be necesary in case one needs to evaluate the magnetic flux over
% the boundary to determine the limiting surface.

% Approximate segment length:
ds = 0.1;
vessel_segmented = SegmentBoundary(vessel,ds);

%% Plot the final result:

% =========================================================================
% Plot final result:
figure; 
hold on
% Plot vessel boundary:
plot(vessel.z ,vessel.r,'ko-')
xlim([-3,3])
ylim([0,2.5])
box on
grid on
% Plot segmented vessel boundary:
for ii = 1:length(vessel_segmented.z)
    plot(vessel_segmented.z(ii),vessel_segmented.r(ii),'r.')
    drawnow
    pause(0.01)
end
