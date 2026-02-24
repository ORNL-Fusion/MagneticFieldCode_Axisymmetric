zHeliconMin = 1.6;
zHeliconMax = 1.9;
rHelicon = 0.06;

zTarget = 4.2;
rTargetMin = 0;
rTargetMax = 0.05;

currents = [530, 2100, 6800, 3500, 430];

[objectives, weightedObjective] = ObjectiveOracle(currents,zHeliconMin,zHeliconMax,rHelicon, zTarget, rTargetMin,rTargetMax)
