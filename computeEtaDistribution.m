function [subs] = computeEtaDistribution(etaMean,SD,nSubjects)
%% computeEtaDistribution 
% draws values from normal distribution with given mean and SD, restricting
% them to be within 3 standard deviations

    subs = [99]; %just a value making function go into while loop
    while ~all(subs>etaMean-(3*SD) & subs<etaMean+(3*SD)) %values should not be
        %more than 3 SDs away from mean
        subs = SD.*randn(nSubjects,1) + etaMean; %draw 20 values from normal 
        %distribution around the eta mean and with a given SD
    end
end