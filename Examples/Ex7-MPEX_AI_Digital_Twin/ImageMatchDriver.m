close all;
clear all;

ntheta = 100;
nz = 50;

zHeliconMin = 1.6;
zHeliconMax = 1.9;
rHelicon = 0.06;

zTarget = 4.2;
rTargetMin = 0;
rTargetMax = 0.05;

currents = [530, 2100, 6800, 3500, 430];

[r,z,Qtot_r, Qtot_z, Qtot_r_Helicon, Qtot_z_Helicon,rHeliconIndex,zGrid,thetaGrid] = ImageMatch(currents,ntheta,nz,zHeliconMin,zHeliconMax,rHelicon);

figure;
hold on;
plot(z,Qtot_r(rHeliconIndex,:));
plot(zGrid,Qtot_r_Helicon);

figure
imagesc(zGrid, thetaGrid, Qtot_r_Helicon'*ones(1,ntheta))
