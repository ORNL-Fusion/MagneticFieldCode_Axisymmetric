function [objectives, weightedObjective] ... 
    = ObjectiveOracle(currents,zHeliconMin,zHeliconMax,rHelicon, zTarget, rTargetMin,rTargetMax)

objectives = [0 0 0]';
weights = [1 1 1]';
weights = weights/sum(weights);

% get the r, z coordinates and the r,z components of the heat flux

[r,z,Qtot_r, Qtot_z] = generate_plasmaProfile_opt(currents);

% recover the mesh size for r and z
dr = r(2)-r(1);
dz = z(2)-z(1);

% find the value of r and the associated index that is closest to the helicon radius 

[rDistance, rHeliconIndex] = min(abs(r - rHelicon));
rHeliconNearest = r(rHeliconIndex);

% find the value of z and the associated index that is closest to the target

[zDistance, zTargetIndex] = min(abs(z - zTarget));
zTargetNearest = z(zTargetIndex);


% compute three objectives

% objective 1 is the L^2 value of the heat flux along the helicon window

integrand = Qtot_r.^2;
integrand = integrand(z > zHeliconMin - dz/2 & z < zHeliconMax + dz/2);
objectives(1) = 2*pi*dz*sum(integrand);

% objective 1 is the L^2 value of the z derivate of the heat flux along the helicon window

integrand = (gradient(Qtot_r(rHeliconIndex,:),dz)).^2;
integrand = integrand(z > zHeliconMin - dz/2 & z < zHeliconMax + dz/2);
objectives(2) = 2*pi*dz*sum(integrand);



% objective 3 is the L^2 value of the heat flux along the target

integrand = (Qtot_z(:,zTargetIndex)).^2;
integrand = integrand(r > rTargetMin - dr/2 & r < rTargetMax + dr/2);
objectives(3) = -pi*dr*sum(integrand.*r(r > rTargetMin - dr/2 & r < rTargetMax + dr/2));

weightedObjective = weights'*objectives;





