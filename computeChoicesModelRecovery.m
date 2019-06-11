function [CP_add, CP_multi, choice_add, choice_multi] = computeChoicesModelRecovery(utFunc,params,subjList,gambs)
%% computeChoicesModelRecovery
% function loops over additive and multiplicative session and 
% 1. computes the utility for the gambles (and their difference) given the 
% utility function and the relevant parameter values (calling 
% computeIsoelasticUtility.m or computeExpectedUtility.m if not time model
% in which case the difference in utility is equal to the difference in
% change in wealth precomputed in variables LinU_GamX_X_add/..._multi)
% 2. computes choice probabilities (for choosing gamble on the left side)
% given a beta value
% 3. computes choices (0 = right, 1 = left) probabilistically realised on
% the basis of the choice probability (calling computeChoiceRealisation.m)

    CP_add = {}; %cell array to be filled with choice probabilities (for 
    %left gamble) for every subject, every trial, additive session
    CP_multi = {}; %as above, for multiplicative session
    choice_add = {}; %probabilistically realised choices (left=1, right=0)
    %based on choices probabilities, for every subject, every trial,
    %additive session
    choice_multi = {}; %as above, for multiplicative session
    p =  params(end); %probability of gamble outcomes (always 0.5)
    for i = subjList 
        for j = 1:2 %1 = additive session, 2 = multiplicative session
            %Prospect theory
            if strcmpi(utFunc,'PT')
                alpha = params(1);
                lambda = params(2);
                if j == 1 %additive session
                    beta = params(3);
                    CP_add{i} = nan(size(gambs{1}{i}));
                    for t = 1:numel(gambs{5}{i})
                        Ut_PT_left = computeExpectedUtility({gambs{1}{i}(t),gambs{2}{i}(t)},alpha,alpha,lambda,[p,p]);
                        Ut_PT_right = computeExpectedUtility({gambs{3}{i}(t),gambs{4}{i}(t)},alpha,alpha,lambda,[p,p]);
                        CP_add{i}(t) = 1/(1+exp(-((Ut_PT_left-Ut_PT_right)*beta)));
                        choice_add{i}(t) = computeChoiceRealisation(CP_add{i}(t));
                    end
                else %multiplicative session
                    beta = params(4);
                    CP_multi{i} = nan(size(gambs{5}{i}));
                    for t = 1:numel(gambs{1}{i})
                        Ut_PT_left = computeExpectedUtility({gambs{5}{i}(t),gambs{6}{i}(t)},alpha,alpha,lambda,[p,p]);
                        Ut_PT_right = computeExpectedUtility({gambs{7}{i}(t),gambs{8}{i}(t)},alpha,alpha,lambda,[p,p]);
                        CP_multi{i}(t) = 1/(1+exp(-((Ut_PT_left-Ut_PT_right)*beta)));
                        choice_multi{i}(t) = computeChoiceRealisation(CP_multi{i}(t));
                    end
                end
            %Isoelastic utility
            elseif strcmpi(utFunc,'isoUt')
                eta = params(1);
                if j == 1 %additive session
                    beta = params(2);
                    CP_add{i} = nan(size(gambs{1}{i}));
                    for t = 1:numel(gambs{1}{i})
                        Ut_iso_left = computeIsoelasticUtility(eta,gambs{9}(i),{gambs{1}{i}(t),gambs{2}{i}(t)});
                        Ut_iso_right = computeIsoelasticUtility(eta,gambs{9}(i),{gambs{3}{i}(t),gambs{4}{i}(t)});
                        CP_add{i}(t) = 1/(1+exp(-((Ut_iso_left-Ut_iso_right)*beta)));
                        choice_add{i}(t) = computeChoiceRealisation(CP_add{i}(t));
                    end
                else %multiplicative session
                    beta = params(3);
                    CP_multi{i} = nan(size(gambs{5}{i}));
                    for t = 1:numel(gambs{5}{i})
                        Ut_iso_left = computeIsoelasticUtility(eta,gambs{10}(i),{gambs{5}{i}(t),gambs{6}{i}(t)});
                        Ut_iso_right = computeIsoelasticUtility(eta,gambs{10}(i),{gambs{7}{i}(t),gambs{8}{i}(t)});
                        CP_multi{i}(t) = 1/(1+exp(-((Ut_iso_left-Ut_iso_right)*beta)));
                        choice_multi{i}(t) = computeChoiceRealisation(CP_multi{i}(t));
                    end
                end
            %Time
            else 
                if j == 1 %additive session
                    beta = params(1);
                    CP_add{i} = nan(size(gambs{1}{i}));
                    for t = 1:numel(gambs{1}{i})
                        CP_add{i}(t) = 1/(1+exp(-((gambs{1}{i}(t))*beta)));
                        choice_add{i}(t) = computeChoiceRealisation(CP_add{i}(t));
                    end
                else %multiplicative session
                    beta = params(2);
                    CP_multi{i} = nan(size(gambs{2}{i}));
                    for t = 1:numel(gambs{2}{i})
                        CP_multi{i}(t) = 1/(1+exp(-((gambs{2}{i}(t))*beta)));
                        choice_multi{i}(t) = computeChoiceRealisation(CP_multi{i}(t));
                    end
                end
            end
        end
    end
end