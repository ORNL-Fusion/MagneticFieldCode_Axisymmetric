function [r,z,Qperp_r,Qperp_z] = ObjectiveOracle(params)


zHeliconMin = 1.6;
zHeliconMax = 1.9;
rHelicon = 0.06;

zTarget = 4.2
rTargetMin = 0;
rTargetMax = 0.05;




cd 'MagneticFieldCode_Axisymmetric-master/Examples/Ex4-ProtoMPEX_FluxMapping';
bFieldFileName = '../../../analytical_model/bfield_protoMPEX.nc';
bFieldFileNameLoc = 'bfield_protoMPEX.nc';
Ex4_FluxMappingVarious_Opt2(params,bFieldFileName);
cd '../../../analytical_model'
[r,z,Qperp_r,Qperp_z] = profiles_protoMPEX_wPerpFlux_Opt2(bFieldFileNameLoc);

dr = r(2)-r(1);
dz = z(2)-z(1);

% Find first objective
% 1) find rindex that is closest to rHelicon

[rDistance, rHeliconIndex] = min(abs(r - rHelicon));
rHeliconNearest = r(rHeliconIndex);

[zDistance, zTargetIndex] = min(abs(z - zTarget));
zTargetNearest = z(zTargetIndex);

integrand_1 = (gradient(Qperp_r(rHeliconIndex,:),dr)).^2;
integrand_1 = integrand_1(z > zHeliconMin & z < zHeliconMax);

integrand_2 = (Qperp_z(:,zTargetIndex)).^2;
integrand_2 = integrand_2(r > rTargetMin & r < rTargetMax);

objective_1 = dz*sum(integrand_1);
objective_2 = dr*sum(integrand_2.*r(r > rTargetMin & r < rTargetMax));

objective = objective_1 + objective_2
;


% dQperp_dz = gradient(Qperp,dz);


