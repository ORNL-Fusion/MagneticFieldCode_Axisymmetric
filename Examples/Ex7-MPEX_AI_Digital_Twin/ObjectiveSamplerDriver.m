
zHeliconMin = 1.6;
zHeliconMax = 1.9;
rHelicon = 0.06;

zTarget = 4.2;
rTargetMin = 0;
rTargetMax = 0.05;

numSamples = 2;
width = 0.2;

coilCurrentsNominal = [530 2100 6800 3500 430];
numCoils = length(coilCurrentsNominal);
coilCurrentsRangees = zeros(numSamples,numCoils);
for i=1:numCoils
    nominal = coilCurrentsNominal(i);
    temp = linspace((1-width/2)*nominal, (1+width/2)*nominal, numSamples);
    coilCurrentsRangees(:,i) = temp; 
end

counter = 0;

currentList = zeros(numSamples,numCoils);

for i=1:numSamples
    for ii=1:numSamples
        for iii=1:numSamples
            for iv=1:numSamples
                for v = 1:numSamples
                    index = [i ii iii iv v]';
                    counter = counter +1;
                    for j=1:numCoils
                        currentList(counter,j) = coilCurrentsRangees(index(j),j);
                    end    
                end
            end
        end
    end
 end

[currentList, objectiveList] = ObjectiveSampler(currentList,zHeliconMin,zHeliconMax,rHelicon, zTarget, rTargetMin,rTargetMax);

newMin = Inf;
for i = 1:counter    
    oldMin = newMin
    nextObjective = objectiveList(i)
    newMin = min(oldMin,nextObjective)
end
    
    

