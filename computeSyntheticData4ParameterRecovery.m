%% Parameter Recovery script
%This scripts "presents" synthetic agents with the gambles presented to the
%subjects. The agents operate with the isoelastic utility function with 
%different parameter values. The script generates choice probabilities 
%based on the function with different parameter values, and then 
%probabilistically realises left-right choices. Can then be used to be fed 
%into HLM script to see whether parameters can be recovered.
%The script computeSyntheticData4ParameterRecovery.m is used to create 
%synthetic agents with only one type of utility function (isoelastic), but
%where we loop through different combinations of eta parameter values for 
%the two different sessions.
%The script computeSyntheticData4ModelRecovery.m is used to create 
%synthetic agents with different utility functions.

%Script now generates 9 x 19 synthetic agents, based on the gambles
%presented to our 19 subjects. The 19 agents have slightly different etas 
%drawn from normal distributions, e.g. 19 etas around 0 for the additive 
%and 0.5 for the multiplicative session, varying around these values with 
%a given SD defined in the script

%% Setting variables determining how script is run
[startDir,~] = fileparts(mfilename('fullpath'));
saveFiles = 1; %if saving of results during debugging not wanted, set to 0
fixedBeta = 1; %if we want to maximize discriminability between parameters,
%we set beta to be high (= 500) at a fixed level for all combinations of
%etas. If fixedBeta is set to 0, the script instead estimates a beta for 
%each eta from linear regression of the estimated beta values from all
%subjects on their estimated eta values 

%% Load relevant variablbes from allData.mat
load(fullfile(startDir,'data','allData.mat'),...
    'LinU_Gam1_1_add', 'LinU_Gam1_2_add', 'LinU_Gam2_1_add', 'LinU_Gam2_2_add',...
    'LinU_Gam1_1_multi', 'LinU_Gam1_2_multi', 'LinU_Gam2_1_multi', 'LinU_Gam2_2_multi',...
    'LogU_Gam1_1_add', 'LogU_Gam1_2_add', 'LogU_Gam2_1_add', 'LogU_Gam2_2_add',...
    'LogU_Gam1_1_multi', 'LogU_Gam1_2_multi', 'LogU_Gam2_1_multi', 'LogU_Gam2_2_multi',...
    'delta_EU_Lin_add', 'delta_EU_Log_add', 'delta_EU_Lin_multi', 'delta_EU_Log_multi',...
    'NoBrainerChoiceCorrect_add', 'NoBrainerChoiceCorrect_multi', ...
    'Wealth_add', 'Wealth_multi');

nSubjects = numel(LinU_Gam1_1_add);
subjList = 1:nSubjects;

%% Retrieving randomly distributed etas from normal distribution
etasAdd = -0.5:0.5:1.5;
etasMulti = -0.5:0.5:1.5;
etas = zeros(numel(etasAdd)*numel(etasMulti),numel(subjList),2); %9 groups 
%of 19 subjects with different combinations of 2 etas (add and multi)

SD = 0.1;
counter = 1;
for a = 1:numel(etasAdd)
    for m = 1:numel(etasMulti)
        etas(counter,:,1) = computeEtaDistribution(etasAdd(a),SD,nSubjects);
        etas(counter,:,2) = computeEtaDistribution(etasMulti(m),SD,nSubjects);
        counter = counter +1;
    end
end

%% Calculate beta values to pair with eta
%Etas and betas (in log) have a close-to-linear relationship (as can be 
%seen when plotting the estimated etas and betas from our subjects (mapEtas
%and log(mapBetas)) against each other). Thus, we here run a linear
%regression, predicting betas with the given eta value, as estimated for
%each subject. Then we take the randomly sampled eta values and calculate
%the corresponding beta values from the parameter values from the
%regression

switch fixedBeta
    case 0
        load(fullfile(startDir,'data','data_mapBetas.mat'));
        load(fullfile(startDir,'data','data_mapEtas.mat'));

        mapEtaValsAdd = mapEtaVals(:,1); %maximum a-posteriori estimates for the 
        %etas for all subjects in the additive session
        mapEtaValsAdd = mapEtaValsAdd(mapEtaValsAdd~=0); %remove subject with 0 values
        mapEtaValsMulti = mapEtaVals(:,2); %maximum a-posteriori estimates for the 
        %etas for all subjects in the multiplicative session
        mapEtaValsMulti = mapEtaValsMulti(mapEtaValsMulti~=0);%remove subject with 0 values

        mapBetaValsAdd = mapBetaVals(:,1);%maximum a-posteriori estimates for the 
        %betas for all subjects in the additive session
        mapBetaValsAdd = mapBetaValsAdd(mapBetaValsAdd~=0);%remove subject with 0 values
        mapBetaValsMulti = mapBetaVals(:,2);%maximum a-posteriori estimates for the 
        %betas for all subjects in the multiplicative session
        mapBetaValsMulti = mapBetaValsMulti(mapBetaValsMulti~=0);%remove subject with 0 values

        xAdd = [ones(numel(mapEtaValsAdd),1) mapEtaValsAdd]; %add vector of 1's to
        %eta vector in order to get intercept and slope from regression
        yAdd = log(mapBetaValsAdd); %transform beta values into log space
        bAdd = xAdd\yAdd; %run linear regression

        xMulti = [ones(numel(mapEtaValsMulti),1) mapEtaValsMulti];
        yMulti = log(mapBetaValsMulti);
        bMulti = xMulti\yMulti;

        betasAdd = exp(bAdd(1) + (etas(:,:,1) * bAdd(2))); %calculate beta values 
        %from randomly sampled eta values, applying parameters from regression (and
        %tranforming beta values back into non log space)
        betasMulti = exp(bMulti(1) + (etas(:,:,2) * bMulti(2)));
        betas = cat(3,betasAdd,betasMulti);
    case 1
        betas = ones(size(etas))*500;
end

if saveFiles
    save(fullfile(startDir,'data','data_etasForParamRecov.mat'),'etas');
    save(fullfile(startDir,'data','data_betasForParamRecov.mat'),'betas');
end

close all;

%% Plot distribution of randomly drawn etas and their betas

h1 = figure(1);
for i = 1:numel(etasAdd)*numel(etasMulti)
    subplot(numel(etasAdd),numel(etasMulti),i)
    histogram(etas(i,:,1),'FaceColor','b','FaceAlpha',.5)
    hold on 
    histogram(etas(i,:,2),'FaceColor','y','FaceAlpha',.5)
    text(-1,14,sprintf('%d',i));
    if i == 9
        legend('Additive','Multiplicative')
    end
    xlim([-0.5,1.5]);
    ylim([0,16]);
end
title('Etas')
set(h1,'color','w');
set(h1,'Position',[0 0 1000 1000]);

if saveFiles
    print(fullfile(startDir,'figs','EtasForParamRecov'),'-dpng');
end

h2 = figure(2);
for i = 1:numel(etasAdd)*numel(etasMulti)
    subplot(5,5,i)
    histogram(betas(i,:,1),'FaceColor','b','FaceAlpha',.5)
    hold on 
    histogram(betas(i,:,2),'FaceColor','y','FaceAlpha',.5)
    text(-1,14,sprintf('%d',i));
    if i == 9
        legend('Additive','Multiplicative')
    end
    ylim([0,16]);
end
title('Betas')
set(h2,'color','w');
set(h2,'Position',[20 0 1000 1000]);

if saveFiles
    print(fullfile(startDir,'figs','BetasForParamRecov'),'-dpng');
end

%% For debugging: Plot betas against etas
addEtasPlot = etas(:,:,1);
addBetasPlot = betas(:,:,1);
multEtasPlot = etas(:,:,2);
multBetasPlot = betas(:,:,2);

figure(4)
plot(reshape(addEtasPlot,prod(size(addEtasPlot)),1),reshape(log(addBetasPlot),prod(size(addBetasPlot)),1),'.')
hold on
plot(reshape(multEtasPlot,prod(size(multEtasPlot)),1),reshape(log(multBetasPlot),prod(size(multBetasPlot)),1),'.')
legend('add' ,'mult');
ylabel('betas (beta = b0 + b1*eta)');
xlabel('etas (randomly drawn)');
title('Synthetic agents');
set(figure(4),'color','w');
set(figure(4),'Position',[40 0 1000 1000]);

figure(5);
plot(mapEtaVals(:,1),log(mapBetaVals(:,1)),'.')
hold on
plot(mapEtaVals(:,2),log(mapBetaVals(:,2)),'.')
legend('add' ,'mult');
ylabel('betas (beta = b0 + b1*eta)');
xlabel('etas (randomly drawn)');
title('Actual subjects');
set(figure(5),'color','w');
set(figure(5),'Position',[60 0 1000 1000]);

%% Compute choices for isoelastic utility with different parameter values
%the function "computeChoicesModelRecovery" returns choice probabilities 
%and synthetic choices probabilistically realized on the basis of the 
%choice probabilities, for both sessions
counter = 1;
p = 0.5; %input for isoelastic utility function which allows estimation of 
%utilities with probabilities different from 0.5, which is not needed in
%this paradigm
for a = 1:numel(etasAdd)
    for m = 1:numel(etasMulti)
        Choice_add_synth = {}; %Initialise variables
        Choice_multi_synth = {}; %Initialise variables
        CP_add_synth = {}; %Initialise variables
        CP_multi_synth = {}; %Initialise variables
        fprintf('%d) Iso eta add: %.1f, eta multi: %.1f\n',counter,etasAdd(a),etasMulti(m));

        %Output variables: 312 choice probabilities for each trial for 
        %multi and add session and 312 probabilistically realized choices 
        %for both sessions
        %Function input: (1) utility function, 
        %(2) parameters (eta multi, eta add, beta multi, beta add, p (always 0.5)),
        %(3) subjects (here: 1),
        %(4) outcomes: for multi session, we want to use absolute change in
        %wealth (e.g. 1000 (wealth) * 1.5 = 1500 -> change in wealth 500)
        %which is precalculated in variables LinU_Gamx_x_multi
        %for add session, that change is the same as the growth increment
        %itself (i.e. Gam1_1_add and LinU_Gam1_1_add are the same)
        for i = 1:size(etas,2) %for n subjects 
            etaAdd = etas(counter,i,1); %one of the 19 etas which were drawn from normal distribution
            etaMulti = etas(counter,i,2);
            betaAdd = betas(counter,i,1);
            betaMulti = betas(counter,i,2);
            
            if i == 1 %whether to create plot of utility against change in wealth
                plotUt = 1;
            else
                plotUt = 0;
            end
            
            [CP_left_iso_add,CP_left_iso_multi,...
            choice_left_iso_add,choice_left_iso_multi] = ...
            computeChoicesParameterRecovery([etaAdd,etaMulti,betaAdd,betaMulti,p],...
            {LinU_Gam1_1_add{i},LinU_Gam1_2_add{i},LinU_Gam2_1_add{i},LinU_Gam2_2_add{i},...
            LinU_Gam1_1_multi{i},LinU_Gam1_2_multi{i},LinU_Gam2_1_multi{i},LinU_Gam2_2_multi{i},...
            Wealth_add(i),Wealth_multi(i)},...
            plotUt,counter);
        
            %Add output of above function "computeChoices" to cell array
            Choice_add_synth{end+1} = choice_left_iso_add{:};
            Choice_multi_synth{end+1} = choice_left_iso_multi{:};
            CP_add_synth{end+1} = CP_left_iso_add{:};
            CP_multi_synth{end+1} = CP_left_iso_multi{:};
            
            %Plot for debugging, showing choice probabilities
            if i == 1
                figure(3);
                subplot(5,5,counter);
                plot(CP_left_iso_add{:},'.');
                hold on;
                plot(CP_left_iso_multi{:},'.');
                ylim([0,1]);
                legend('add','multi');
                title(sprintf('%d) add: %.3f,  multi: %.3f',counter,etaAdd,etaMulti));
                set(figure(3),'color','w');
                set(figure(3),'Position',[80 0 1000 1000]);
            end
        end
        
        %Having created choice data for all subjects with this combination 
        %of eta parameters, now save synthetic data files with all the 
        %relevant original variables for all subjects, but with choices 
        %replaced with the ones from the synthetic agents realized on the 
        %basis of the isoelastic utility function with these eta values
        Choice_add = Choice_add_synth;
        Choice_multi = Choice_multi_synth;
        
        if saveFiles
            save(fullfile(startDir,'data',sprintf('allData_synth_paramRecov_eta%s-%s_beta500.mat',...
            num2str(etasAdd(a)*10,'%02d\n'),num2str(etasMulti(m)*10,'%02d\n'))),...
            'LinU_Gam1_1_add', 'LinU_Gam1_2_add', 'LinU_Gam2_1_add', 'LinU_Gam2_2_add',...
            'LinU_Gam1_1_multi', 'LinU_Gam1_2_multi', 'LinU_Gam2_1_multi', 'LinU_Gam2_2_multi',...
            'LogU_Gam1_1_add', 'LogU_Gam1_2_add', 'LogU_Gam2_1_add', 'LogU_Gam2_2_add',...
            'LogU_Gam1_1_multi', 'LogU_Gam1_2_multi', 'LogU_Gam2_1_multi', 'LogU_Gam2_2_multi',...
            'delta_EU_Lin_add', 'delta_EU_Log_add', 'delta_EU_Lin_multi', 'delta_EU_Log_multi',...
            'NoBrainerChoiceCorrect_add', 'NoBrainerChoiceCorrect_multi', ...
            'Wealth_add', 'Wealth_multi', 'Choice_add', 'Choice_multi');  
        end
        
        counter = counter + 1;
    end
end


