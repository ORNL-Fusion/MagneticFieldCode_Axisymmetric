function [currentList, objectiveList] = ObjectiveSampler(currentList,zHeliconMin,zHeliconMax,rHelicon, zTarget, rTargetMin,rTargetMax)


[numSamples,numCoils]  = size(currentList);

objectiveList = zeros(numSamples,1);

for i=1:numSamples
    currents = currentList(i,:);
    [~, weightedObjective] ... 
        = ObjectiveOracle(currents,zHeliconMin,zHeliconMax,rHelicon, zTarget, rTargetMin,rTargetMax);
    
    objectiveList(i) = weightedObjective;

end



    
    

