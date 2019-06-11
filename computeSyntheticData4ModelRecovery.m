%% Model Recovery script
%This scripts "presents" synthetic agents with the gambles presented to the
%subjects. The agents operate with utility functions with different
%parameter values. The script generates choice probabilities based on these
%functions, and then probabilistically realises left-right choices. Can
%then be used to feed into HBM script to see whether parameters can be
%recovered

%The script computeSyntheticData4ModelRecovery.m is used to create 
%synthetic agents with different utility functions.
%The script computeSyntheticData4ParameterRecovery.m is used to create 
%synthetic agents with only one type of utility function (isoelastic), but
%where we loop through different combinations of eta parameter values for 
%the two different sessions.

%Script now generates 27 synthetic agents, making decisions with different
%utility models and different parameter values, based on the gambles
%presented to subject 1. 
%The file allData_synth_modelRecov.mat will contain choices in following
%order (we use beta = exp(-1) for additive session and beta = exp(5) for 
%multiplicative session):
% 1) PT: alpha 0.6, lambda 1.0
% 2) PT: alpha 0.6, lambda 2.0
% 3) PT: alpha 0.6, lambda 3.0
% 4) PT: alpha 0.8, lambda 1.0
% 5) PT: alpha 0.8, lambda 2.0
% 6) PT: alpha 0.8, lambda 3.0
% 7) PT: alpha 0.9, lambda 1.0
% 8) PT: alpha 0.9, lambda 2.0
% 9) PT: alpha 0.9, lambda 3.0
% 10) Iso: eta -0.5
% 11) Iso: eta -0.2
% 12) Iso: eta 0.0
% 13) Iso: eta 0.2
% 14) Iso: eta 0.5
% 15) Iso: eta 0.8
% 16) Iso: eta 1.0
% 17) Iso: eta 1.2
% 18) Iso: eta 1.5
% 19) time
% 20) time
% 21) time
% 22) time
% 23) time
% 24) time
% 25) time
% 26) time
% 27) time


%% Load data, define parameters for the different utility functions
[startDir,~] = fileparts(mfilename('fullpath'));
load(fullfile(startDir,'data','allData.mat'));
alphas = [0.6,0.75,0.9];
lambdas = [1,2,3];
etas = linspace(-0.5,1.5,9);
betas = [exp(-1),exp(5)];
p = 0.5; %input for isoelastic and expected utility functions which allow
%estimation of utilities with probabilities different from 0.5, which is 
%not needed in this paradigm

subjList = [1];
Choice_multi_synth = {};
Choice_add_synth = {};

counter = 0;

%% Compute choices for prospect theory with different alphas and lambdas
%the function "computeChoicesModelRecovery.m" returns choice probabilities 
%and synthetic choices probabilistically realized on the basis of the 
%choice probabilities, for both sessions
for a = 1:numel(alphas)
    for g = 1:numel(lambdas)
        counter = counter + 1;
        alpha = alphas(a);
        lambda = lambdas(g);

        fprintf('%d) PT: alpha %.1f, lambda %.1f\n',counter,alpha,lambda);

        [CP_left_PT_add,CP_left_PT_multi,choice_left_PT_add,choice_left_PT_multi] = ...
            computeChoicesModelRecovery('PT',[alpha,lambda,betas,p],subjList,...
            {LinU_Gam1_1_add,LinU_Gam1_2_add,LinU_Gam2_1_add,LinU_Gam2_2_add,...
            LinU_Gam1_1_multi,LinU_Gam1_2_multi,LinU_Gam2_1_multi,LinU_Gam2_2_multi});

        Choice_add_synth{end+1} = choice_left_PT_add{:};
        Choice_multi_synth{end+1} = choice_left_PT_multi{:};
    end
end

%% Compute choices for isoelastic utility with different etas
for e = 1:numel(etas)
    counter = counter + 1;
    eta = etas(e);

    fprintf('%d) Iso: eta %.1f\n',counter,eta);

    [CP_left_iso_add,CP_left_iso_multi,choice_left_iso_add,choice_left_iso_multi] = ...
        computeChoicesModelRecovery('isoUt',[eta,betas,p],subjList,...
        {LinU_Gam1_1_add,LinU_Gam1_2_add,LinU_Gam2_1_add,LinU_Gam2_2_add,...
        LinU_Gam1_1_multi,LinU_Gam1_2_multi,LinU_Gam2_1_multi,LinU_Gam2_2_multi,...
        Wealth_add(1),Wealth_multi(1)});

    Choice_add_synth{end+1} = choice_left_iso_add{:};
    Choice_multi_synth{end+1} = choice_left_iso_multi{:};
end

%% Compute choices for time criterion model
% This uses 9 times choices based on linear utility in the additive session
% and logarithmic utility in the multiplicative session. Choices are
% somewhat variable between the 9 synthetic agents due to probabilistic
% realisation of choice
betas = repmat(betas,9,1);
for b = 1:size(betas,1)
    counter = counter + 1;
    
    fprintf('%d) time\n',counter);
    
    [CP_left_time_add,CP_left_time_multi,choice_left_time_add,choice_left_time_multi] = ...
        computeChoicesModelRecovery('time',[betas(b,:),p],subjList,{delta_EU_Lin_add,delta_EU_Log_multi});
    
    Choice_add_synth{end+1} = choice_left_time_add{:};
    Choice_multi_synth{end+1} = choice_left_time_multi{:};
end

%% Save synthetic data file 
%Data file with 27 different choices, but otherwise 27 times the same
%information from subject 1 (gambles presented to subject 1 and derived
%differences in utility)
for i = 1:numel(Choice_add_synth)
    LinU_Gam1_1_add{i} = LinU_Gam1_1_add{1};
    LinU_Gam1_2_add{i} = LinU_Gam1_2_add{1};
    LinU_Gam2_1_add{i} = LinU_Gam1_1_add{1};
    LinU_Gam2_2_add{i} = LinU_Gam2_2_add{1};
    LinU_Gam1_1_multi{i} = LinU_Gam1_1_multi{1};
    LinU_Gam1_2_multi{i} = LinU_Gam1_2_multi{1};
    LinU_Gam2_1_multi{i} = LinU_Gam2_1_multi{1};
    LinU_Gam2_2_multi{i} = LinU_Gam2_2_multi{1};
    LogU_Gam1_1_add{i} = LogU_Gam1_1_add{1};
    LogU_Gam1_2_add{i} = LogU_Gam1_2_add{1};
    LogU_Gam2_1_add{i} = LogU_Gam1_1_add{1};
    LogU_Gam2_2_add{i} = LogU_Gam2_2_add{1};
    LogU_Gam1_1_multi{i} = LogU_Gam1_1_multi{1};
    LogU_Gam1_2_multi{i} = LogU_Gam1_2_multi{1};
    LogU_Gam2_1_multi{i} = LogU_Gam2_1_multi{1};
    LogU_Gam2_2_multi{i} = LogU_Gam2_2_multi{1};
    delta_EU_Lin_add{i} = delta_EU_Lin_add{1};
    delta_EU_Log_add{i} = delta_EU_Log_add{1};
    delta_EU_Lin_multi{i} = delta_EU_Lin_multi{1};
    delta_EU_Log_multi{i} = delta_EU_Log_multi{1};
    NoBrainerChoiceCorrect_add{i} = NoBrainerChoiceCorrect_add{1};
    NoBrainerChoiceCorrect_multi{i} = NoBrainerChoiceCorrect_multi{1};  
    Wealth_add(i) = Wealth_add(1);
    Wealth_multi(i) = Wealth_multi(1);
end
Choice_add = Choice_add_synth;
Choice_multi = Choice_multi_synth;

save(fullfile(startDir,'data','allData_synth_modelRecov.mat'),...
    'LinU_Gam1_1_add', 'LinU_Gam1_2_add', 'LinU_Gam2_1_add', 'LinU_Gam2_2_add',...
    'LinU_Gam1_1_multi', 'LinU_Gam1_2_multi', 'LinU_Gam2_1_multi', 'LinU_Gam2_2_multi',...
    'LogU_Gam1_1_add', 'LogU_Gam1_2_add', 'LogU_Gam2_1_add', 'LogU_Gam2_2_add',...
    'LogU_Gam1_1_multi', 'LogU_Gam1_2_multi', 'LogU_Gam2_1_multi', 'LogU_Gam2_2_multi',...
    'delta_EU_Lin_add', 'delta_EU_Log_add', 'delta_EU_Lin_multi', 'delta_EU_Log_multi',...
    'NoBrainerChoiceCorrect_add', 'NoBrainerChoiceCorrect_multi', ...
    'Wealth_add', 'Wealth_multi', 'Choice_add', 'Choice_multi');