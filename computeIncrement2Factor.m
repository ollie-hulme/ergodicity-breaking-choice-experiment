function [ growthFactors ] = computeIncrement2Factor(growthIncrements,currentWealth)
%% computeIncrement2Factor
% computeIncrement2Factor computes growth factors from additive growth 
% increments and current wealth

for lp=1:length(currentWealth)
growthFactors(:,:,lp)=(growthIncrements+currentWealth(lp))./currentWealth(lp);
end

