function [CP_add, CP_multi, choice_add, choice_multi] = computeChoicesParameterRecovery(params,gambs,plotUt,counter)
%% computeChoicesParameterRecovery
% function loops over additive and multiplicative session and computes 
% 1. computes the utility for the gambles (and their difference) given an 
% eta value (calling computeIsoelasticUtility.m)
% 2. computes choice probabilities (for choosing gamble on the left side)
% given a beta value
% 3. computes choices (0 = right, 1 = left) probabilistically realised on
% the basis of the choice probability (calling computeChoiceRealisation.m)

    CP_add = {}; %cell array to be filled with choice probabilities (for 
    %left gamble) for every subject, every trial, additive session
    CP_multi = {}; %as above, for multiplicative session
    choice_add = {}; %probabilistically realised choices (left=1, right=0)
    %based on choices probabilities, for every subject, every trial
    choice_multi = {};
    Ut_left_add = {};
    Ut_right_add = {};
    Ut_left_multi = {};
    Ut_right_multi = {};
    p =  params(end); %probability of gamble outcomes (always 0.5)
    for j = 1:2 %1 = additive session, 2 = multiplictive
        if j == 1 %additive session
            CP_add{1} = nan(size(gambs{1})); %pre-allocate empty array to be filed with choice probs
            choice_add{1} = nan(size(gambs{1})); %pre-allocate empty array to be filed with choices
            Ut_left_add{1} = nan(3,numel(gambs{1}));
            Ut_right_add{1} = nan(3,numel(gambs{1}));
            eta = params(1); %params contains eta add, eta multi, beta add, beta multi, p (latter not used in this script)
            beta = params(3);
            for t = 1:numel(gambs{1}) %for every trial, retrieve iso 
                %utils for left and right gamble, then apply softmax
                %to calculate choice probabilities, then call
                %"computeChoiceRealisation" to realize choices probabilistically
                [Ut_left, Ut_left_top, Ut_left_bot] = computeIsoelasticUtility(eta,...
                    gambs{9},{gambs{1}(t),gambs{2}(t)});
                [Ut_right, Ut_right_top, Ut_right_bot] = computeIsoelasticUtility(eta,...
                    gambs{9},{gambs{3}(t),gambs{4}(t)});
                Ut_left_add{1}(1,t) = Ut_left;
                Ut_left_add{1}(2,t) = Ut_left_top;
                Ut_left_add{1}(3,t) = Ut_left_bot;
                Ut_right_add{1}(1,t) = Ut_right;
                Ut_right_add{1}(2,t) = Ut_right_top;
                Ut_right_add{1}(3,t) = Ut_right_bot;
                CP_add{1}(t) = 1/(1+exp(-((Ut_left-Ut_right)*beta))); %softmax
                choice_add{1}(t) = computeChoiceRealisation(CP_add{1}(t)); %probabilistically
                %realized choice, see function below
            end
        else %multiplicative session
            CP_multi{1} = nan(size(gambs{5})); %pre-allocate empty array to be filed with choice probs
            choice_multi{1} = nan(size(gambs{5})); %pre-allocate empty array to be filed with choices
            Ut_left_multi{1} = nan(3,numel(gambs{5}));
            Ut_right_multi{1} = nan(3,numel(gambs{5}));
            eta = params(2); %params contains eta add, eta multi, beta add, beta multi, p (latter not used in this script)
            beta = params(4);
            for t = 1:numel(gambs{5}) %for every trial, retrieve iso 
                %utils for left and right gamble, then apply softmax
                %to calculate choice probabilities, then call
                %"computeChoiceRealisation" to realize choices probabilistically
                [Ut_left, Ut_left_top, Ut_left_bot] = computeIsoelasticUtility(eta,...
                    gambs{10},{gambs{5}(t),gambs{6}(t)});
                [Ut_right, Ut_right_top, Ut_right_bot] = computeIsoelasticUtility(eta,...
                    gambs{10},{gambs{7}(t),gambs{8}(t)});
                Ut_left_multi{1}(1,t) = Ut_left;
                Ut_left_multi{1}(2,t) = Ut_left_top;
                Ut_left_multi{1}(3,t) = Ut_left_bot;
                Ut_right_multi{1}(1,t) = Ut_right;
                Ut_right_multi{1}(2,t) = Ut_right_top;
                Ut_right_multi{1}(3,t) = Ut_right_bot;
                CP_multi{1}(t) = 1/(1+exp(-((Ut_left-Ut_right)*beta))); %softmax
                choice_multi{1}(t) = computeChoiceRealisation(CP_multi{1}(t)); %probabilistically
                %realized choice, see function below
            end
        end
    end
    if plotUt
        figure(6);
        subplot(5,5,counter);
        yyaxis left
        p1 = plot(gambs{1},Ut_left_add{1}(2,:),'.',gambs{2},Ut_left_add{1}(3,:),'.',...
            gambs{3},Ut_right_add{1}(2,:),'.',gambs{4},Ut_right_add{1}(3,:),'.');
        ylabel('utility add');
        hold on;
        yyaxis right
        p2 = plot(gambs{5},Ut_left_multi{1}(2,:),'.',gambs{6},Ut_left_multi{1}(3,:),'.',...
            gambs{7},Ut_right_multi{1}(2,:),'.',gambs{8},Ut_right_multi{1}(3,:),'.');
        ylabel('utility multi');
        xlabel('change in wealth')
        if counter == 9
            legend([p1(1),p2(1)],'additive','multiplicative')
        end
        title(sprintf('Eta add: %.3f, eta multi: %.3f',params(1),params(2)));
        set(figure(6),'Position',[0 50 1920 1000]);
        set(figure(6),'color','w');
    end
end