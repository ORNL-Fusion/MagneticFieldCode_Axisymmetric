% Inch to meter conversion factor
in2m = 0.0254;

% Dump plate geometry:
rDump = 15.75*in2m/2;
zDump = 0.5;

% Central chamber:
zCC1 = coil{6}.z + 0.5*coil{6}.dz;
zCC2 = coil{7}.z - 0.5*coil{7}.dz;
rCC  = 24*in2m/2;

% Vaccum vessel:
rVac1 = 5.834*in2m/2; % Vacuum vessel rad that is on either side of ALN window
rVac2a = 19.25*in2m/2; % Wide sections dowmstream of central chamber to end
rVac2b = 4.272*2*in2m/2;  % Narrow sections are just inside of coil
zVac1  = zCC1;
zVac2a = zCC2;
zVac2b = coil{end}.z + 0.5*coil{end}.dz;

% Helicon window:
L =  11.8*in2m; % L is 11.8" from Meitner
r1 = 4.95*in2m/2; % diameter is 5.1 inches --> update 4/21/16.  Nominal is ~4.95;
r2 = rVac1;
z1 = 0.5*(coil{3}.z + coil{4}.z) - L/2;
z2 = 0.5*(coil{3}.z + coil{4}.z) + L/2;
heliconWindow.r = [r2,r1,r1,r2];
heliconWindow.z = [z1,z1,z2,z2];

% MPEX-like limiter:
limiterLength = 30e-2;
limiterWidth  = 3e-3;
z1 = heliconWindow.z(end) + 1e-2;
z2 = z1 + limiterLength;
r1 = 2.5*in2m - limiterWidth;
r2 = rVac1;

limiter.r = [r2,r1,r1,r2];
limiter.z = [z1,z1,z2,z2];

% Skimmer:
r1 = (7/100)/2;
r2 = rVac1;
z1 = coil{5}.z + 0.0702 - 0.5e-2;
z2 = coil{5}.z + 0.0702 + 0.5e-2;
skimmer.r = [r2,r1,r1,r2];
skimmer.z = [z1,z1,z2,z2];

% ECH heating region:
r1 = rVac2b;
r2 = rVac2a;
z1 = coil{8}.z + 0.5*coil{8}.dz;
z2 = coil{9}.z - 0.5*coil{9}.dz;
echSection.r = [r1,r2,r2,r1];
echSection.z = [z1,z1,z2,z2];

% ICH Sleeve: -------------------------------------------------------------
r1 = rVac2b;
r2 = 0.08/2; % ID = 80 mm
z1 = coil{09}.z - 0.5*coil{09}.dz;
z2 = coil{10}.z + 0.5*coil{10}.dz;
ichSleeve.r = [r1,r2,r2,r1];
ichSleeve.z = [z1,z1,z2,z2];

% Space between coil 10-11
r1 = rVac2b;
r2 = rVac2a;
z1 = coil{10}.z + 0.5*coil{10}.dz;
z2 = coil{11}.z - 0.5*coil{11}.dz;
space1.r = [r1,r2,r2,r1];
space1.z = [z1,z1,z2,z2];

% Space between coil 11-12
r1 = rVac2b;
r2 = rVac2a;
z1 = coil{11}.z + 0.5*coil{11}.dz;
z2 = coil{12}.z - 0.5*coil{12}.dz;
space2.r = [r1,r2,r2,r1];
space2.z = [z1,z1,z2,z2];

% Space between coil 12-13
r1 = rVac2b;
r2 = rVac2a;
z1 = coil{end-1}.z + 0.5*coil{end-1}.dz;
z2 = coil{end}.z - 0.5*coil{end}.dz;
space3.r = [r1,r2,r2,r1];
space3.z = [z1,z1,z2,z2];

vessel.r = [0    ,rDump,rDump    ,rVac1     ,rVac1,rCC  ,rCC   ,rVac2b,rVac2b];
vessel.z = [zDump,zDump,zDump+0.2,zDump+0.2 ,zVac1,zVac1,zVac2a,zVac2a,zVac2b];

vessel = AddComponent(vessel,echSection);
vessel = AddComponent(vessel,ichSleeve);
vessel = AddComponent(vessel,space1);
vessel = AddComponent(vessel,space2);
vessel = AddComponent(vessel,space3);
vessel = AddComponent(vessel,heliconWindow);
vessel = AddComponent(vessel,skimmer);
vessel_wLim = AddComponent(vessel,limiter);

vessel = SegmentBoundary(vessel,0.02);
vessel_wLim = SegmentBoundary(vessel_wLim,0.02);

figure;
plot(vessel.z,vessel.r,'k.-')
xlim([0,5])
ylim([0,1])