function [choice] = computeChoiceRealisation(choiceProb)
%% computeChoiceRealisation
% computes probabilistically realised choices (0's or 1's) based on choice
% probabilities
    if rand < choiceProb
        choice = 1; %e.g. if rand returns 0.6 and CP for  
        %this trial is 0.8, then choice for this trial is left (=1)
    else
        choice = 0;
    end
end