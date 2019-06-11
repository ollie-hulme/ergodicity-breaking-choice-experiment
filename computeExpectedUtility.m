function [ eu ] = computeExpectedUtility(outcome,exponentGain,exponentLoss,lossLambda,probVals)
%% computeExpectedUtility
% This takes two outcome spaces outcome{1} and outcome{2} and computes
% expected utility according to prospect theory givne the input exponents,
% loss aversion parameters (lambda), and probability values (probVals), via
% a weighted sum, as per cumulative prospect theory.

p1=probVals(1);% probability of outcome 1
try
    p2=probVals(2);% probability of outcome 2
    if p1+p2 > 1 % ensure probabilites do not exceed 1
        disp('Error: probabilities sum > 1')
    end
end

%% Values for outcome 1
v1=nan(size(outcome{1}));
v1(outcome{1}>=0)=outcome{1}(outcome{1}>=0).^exponentGain;%gains
v1(outcome{1}<0)=-lossLambda.*(((-outcome{1}(outcome{1}<0)).^exponentLoss));%losses

eu = (v1*p1);

try
    
    %% Values for outcome 2
    v2=nan(size(outcome{2}));
    v2(outcome{2}>=0)=outcome{2}(outcome{2}>=0).^exponentGain;%gains
    v2(outcome{2}<0)=-lossLambda.*(((-outcome{2}(outcome{2}<0)).^exponentLoss));%losses
    
    %% Weighted sum
    eu = (v1*p1) + (v2*p2);
    
end


end

