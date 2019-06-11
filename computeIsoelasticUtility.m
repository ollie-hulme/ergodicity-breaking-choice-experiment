function [isoU, v1, v2] = computeIsoelasticUtility (eta,wealth,outcome)
%% computeIsoelasticUtility
% This takes two outcome spaces outcome{1} and outcome{2} and the initial 
% wealth and computes isoelastic utility according given an eta value

if eta == 1.0
    uInit=log(wealth);
    u1 = log(outcome{1}+wealth); u2 = log(outcome{2}+wealth);
else
    uInit=((wealth.^(1-eta))-1)./(1-eta);%utility of initial wealth
    u1= (((outcome{1}+wealth).^(1-eta))-1)./(1-eta); u2=(((outcome{2}+wealth).^(1-eta))-1)./(1-eta);%utilities of wealth following outcomes
end

%% changes in utilities and their expectation value
v1=u1-uInit;v2=u2-uInit;%changes in utility for each outcome
isoU = (v1*0.5) + (v2*0.5);%expected change in utility

end