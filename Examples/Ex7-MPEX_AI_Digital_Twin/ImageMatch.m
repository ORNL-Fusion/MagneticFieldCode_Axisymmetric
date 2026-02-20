function[r,z,Qtot_r, Qtot_z, Qtot_r_Helicon, Qtot_z_Helicon,rHeliconIndex,zGrid,thetaGrid] = ImageMatch(currents,ntheta,nz,zHeliconMin,zHeliconMax,rHelicon)

[r,z,Qtot_r, Qtot_z] = generate_plasmaProfile_opt(currents);

[rDistance, rHeliconIndex] = min(abs(r - rHelicon));
rHeliconNearest = r(rHeliconIndex);

thetaGrid = linspace(0,360,ntheta+1);
thetaGrid = (thetaGrid(1:end-1) + thetaGrid(2:end))/2;
zGrid = linspace(zHeliconMin,zHeliconMax,nz+1);
zGrid = (zGrid(1:end-1) + zGrid(2:end))/2;

[Z, Theta] = meshgrid(zGrid,thetaGrid);

Qtot_r_Helicon = interp1(z,Qtot_r(rHeliconIndex,:),zGrid);
Qtot_z_Helicon = interp1(z,Qtot_z(rHeliconIndex,:),zGrid);


