%% This code does not work.


function [currentList, ObjectiveList] = ObjectiveSampler(currents,zHeliconMin,zHeliconMax,rHelicon, zTarget, rTargetMin,rTargetMax)
zHeliconMin = 1.6;
zHeliconMax = 1.9;
rHelicon = 0.06;

zTarget = 4.2;
rTargetMin = 0;
rTargetMax = 0.05;

sampleRate = 2;
width = 0.2;

coilCurrentsNominal = [530 2100 6800 3500 430]';
numCoils = length(coilCurrentsNominal);
coilCurrentsRangees = zeros(numCoils, sampleRate);
for i=1:numCoils
    nominal = coilCurrentsNominal(i);
    temp = (linspace((1-width/2)*nominal, (1+width/2)*nominal, sampleRate))';
    coilCurrentsRangees(i,:) = temp; 
end

newMin = Inf;
currents = zeros(numCoils,1);
counter = 0;
 for i=1:sampleRate
    for ii=1:sampleRate
        for iii=1:sampleRate
            for iv=1:sampleRate
                for v = 1:sampleRate
                    index = [i ii iii iv v]';
                    for j=1:numCoils
                        currents(j) = coilCurrentsRangees(j,index(j));
                    end
                    [~, weightedObjective] ...
                        = ObjectiveOracle(currents,zHeliconMin,zHeliconMax,rHelicon, zTarget, rTargetMin,rTargetMax);
                    counter = counter +1
                    current = currents'
                    currentMin = newMin
                    nextObjective = weightedObjective
                    newMin = min(currentMin,weightedObjective)
                end
            end
        end
    end
 end
    
    

