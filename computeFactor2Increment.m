function [ growthIncrements ] = computeFactor2Increment(growthfactors,currentWealth)
%%  computeFactor2Increment
% computeFactor2Increment computes linear changes in wealth as a function 
% of current wealth and growth rate

for lp=1:length(currentWealth)
growthIncrements(:,:,lp) =(growthfactors*currentWealth(lp))-currentWealth(lp);
end


