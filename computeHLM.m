function computeHLM(runModelNum,nBurnin,nSamples,nThin,nChains,subjList,whichJAGS,runPlots,synthMode,doParallel)

%% Hiercharchical Latent Mixture (HLM) model
% This is a general script for running several types of hierarchical
% bayesian model via JAGS. It can run hiearchical latent mixture models in
% which different utility models can be compared via inference on the model
% indicator variables, and it can run without latent mixtures of models for
% instance in order to estimate parameters of a given utility model. It
% takes as input the following:

% runModelNum - a number which selects which JAGS model to run.
% nBurnin - a number which specifies how many burn in samples to run.
% nThin - a number specifying the thinnning number
% nChains - number of chains
% subjList - list of subject numbers to include
% whichJAGS - sets which copy of matjags to run
% runPlots - sets whether to plot data/illustrations(1) or to suppress (0)
% synthMode - sets whether to run on real data (1), synthetic data for model recovery (2) or parameter recovery (3)

%% Set paths
[startDir,~] = fileparts(mfilename('fullpath'));%specify your starting directory here (where this script runs from)
cd(startDir);%move to starting directory
jagsDir=[startDir,'/JAGS'];
addpath(fullfile(pwd,'/matjags'));%set path to include matjags folder
addpath(fullfile(pwd,'/data'));%set path to include data folder
figDir=[startDir,'/figs'];%specify where you want to save figures within starting directory

%% Choose & load data
switch synthMode
    case 1, dataSource ='allData'; %Real experimental data
    case 2, dataSource ='allData_synth_modelRecov';%Synthetic data for model recovery from 3 different models
    case 4, dataSource ='allData_synth_paramRecov_eta-5--5_beta500';%Synthetic data for model recovery with an eta of -0.5 for add and -0.5 for mult and beta of 500
    case 5, dataSource ='allData_synth_paramRecov_eta-5-00_beta500';%Ditto but with the parameters specified in the title, and so on.
    case 6, dataSource ='allData_synth_paramRecov_eta-5-05_beta500';
    case 7, dataSource ='allData_synth_paramRecov_eta-5-10_beta500';
    case 8, dataSource ='allData_synth_paramRecov_eta-5-15_beta500';
    case 9, dataSource ='allData_synth_paramRecov_eta00--5_beta500';
    case 10,dataSource ='allData_synth_paramRecov_eta00-00_beta500';
    case 11,dataSource ='allData_synth_paramRecov_eta00-05_beta500';
    case 12,dataSource ='allData_synth_paramRecov_eta00-10_beta500';
    case 13,dataSource ='allData_synth_paramRecov_eta00-15_beta500';
    case 14,dataSource ='allData_synth_paramRecov_eta05--5_beta500';
    case 15,dataSource ='allData_synth_paramRecov_eta05-00_beta500';
    case 16,dataSource ='allData_synth_paramRecov_eta05-05_beta500';
    case 17,dataSource ='allData_synth_paramRecov_eta05-10_beta500';
    case 18,dataSource ='allData_synth_paramRecov_eta05-15_beta500';
    case 19,dataSource ='allData_synth_paramRecov_eta10--5_beta500';
    case 20,dataSource ='allData_synth_paramRecov_eta10-00_beta500';
    case 21,dataSource ='allData_synth_paramRecov_eta10-05_beta500';
    case 22,dataSource ='allData_synth_paramRecov_eta10-10_beta500';
    case 23,dataSource ='allData_synth_paramRecov_eta10-15_beta500';
    case 24,dataSource ='allData_synth_paramRecov_eta15--5_beta500';
    case 25,dataSource ='allData_synth_paramRecov_eta15-00_beta500';
    case 26,dataSource ='allData_synth_paramRecov_eta15-05_beta500';
    case 27,dataSource ='allData_synth_paramRecov_eta15-10_beta500';
    case 28,dataSource ='allData_synth_paramRecov_eta15-15_beta500';
end
load(dataSource)

%% Set model specific variables
%Set model name, number of utility models, and prior over model indicator
%variable z
switch runModelNum
    case 1 %parameter estimation of eta
        modelName = 'JAGS_Iso_Subjectwise_Conditionwise'; priorName='NoZPrior';
        
    case 2 %model selection between 3 models with parameter expansion of model indicator
        modelName = 'JAGS_StrongModels_Subjectwise_ExpandZ'; priorName='FlatZprior'; pz=repmat(1/12,1,12);%sets flat prior over models [time, pt, iso]
   
end

%% Set key variables
maxTrials=314;%highest number of trials
nTrials=312;%most sessions have 312 trials, only one session was 314 therefore set to 312
nConditions=2;%number of conditions, keep this to 2
doDIC=1;%compute Deviance information criteria? This is the hierarchical equivalent of an AIC, the lower the better
nSubjects=length(subjList);%number of subjects

%% Set bounds of hyperpriors
%hard code the upper and lower bounds of hyperpriors, typically uniformly
%distributed in log space. These values will be imported to JAGS.

%beta - prior on log since cannot be less than 0; note same bounds used for independent priors on all utility models
muLogBetaL=-2.3;muLogBetaU=3.4;muLogBetaM=(muLogBetaL+muLogBetaU)/2; %bounds on mean of distribution log beta
sigmaLogBetaL=0.01;sigmaLogBetaU=sqrt(((muLogBetaU-muLogBetaL)^2)/12);sigmaLogBetaM=(sigmaLogBetaL+sigmaLogBetaU)/2;%bounds on the std of distribution of log beta

%Alpha - prior on log since cannot be less than 0
muLogAlphaL=-2.3;muLogAlphaU=0;muLogAlphaM=(muLogAlphaL+muLogAlphaU)/2;%bounds on mean of distribution of log Alpha
sigmaLogAlphaL=0.01;sigmaLogAlphaU=sqrt(((muLogAlphaU-muLogAlphaL)^2)/12);sigmaLogAlphaM=(sigmaLogAlphaL+sigmaLogAlphaU)/2; %bounds on std of distribution of log Alpha

%Lambda - prior on log since cannot be less than 0
muLogLambdaL=0;muLogLambdaU=1.6;muLogLambdaM=(muLogLambdaL+muLogLambdaU)/2;%bounds on the mean of distribution of log Lambda
sigmaLogLambdaL=0.01;sigmaLogLambdaU=sqrt(((muLogLambdaU-muLogLambdaL)^2)/12);sigmaLogLambdaM=(sigmaLogLambdaL+sigmaLogLambdaU)/2;%bounds on the std of log Lambda

%eta
muEtaL=-2.5;muEtaU=2.5;muEtaM=(muEtaL+muEtaU)/2;%bounds on mean of distribution of eta
sigmaEtaL=0.01;sigmaEtaU=sqrt(((muEtaU-muEtaL)^2)/12);sigmaEtaM=(sigmaEtaL+sigmaEtaU)/2;%bounds on std of eta

%% Visualise priors at the edges and centre of hyperprior space
if runPlots==1
    figName='Hyperprior limits';
    nVals=1000;x=linspace(-7,8,nVals);%space of log and linear units to visualise
    figure;nrows=4;ncols=2;
    
    %betas calculate normals and normalise
    y1=pdf('norm',x,muLogBetaL,sigmaLogBetaL);y1=y1/max(y1);y2=pdf('norm',x,muLogBetaM,sigmaLogBetaL);y2=y2/max(y2);
    y3=pdf('norm',x,muLogBetaU,sigmaLogBetaL);y3=y3/max(y3);y4=pdf('norm',x,muLogBetaL,sigmaLogBetaU);y4=y4/max(y4);
    y5=pdf('norm',x,muLogBetaM,sigmaLogBetaU);y5=y5/max(y5);y6=pdf('norm',x,muLogBetaU,sigmaLogBetaU);y6=y6/max(y6);
    y7=exp(y1);y7=y7/max(y7);y8=exp(y2);y8=y8/max(y8);y9=exp(y3);y9=y9/max(y9);y10=exp(y4);y10=y10/max(y10);y11=exp(y5);y11=y11/max(y11);y12=exp(y6);y12=y12/max(y12);
    subplot(nrows,ncols,1),plot(x,[y1;y2;y3;y4;y5;y6]),xlim([-7,8]);ntitle('log beta')
    subplot(nrows,ncols,2),plot(exp(x),[y7;y8;y9;y10;y11;y12]),xlim([0,60]);ntitle('beta' )
    
    %Alphas calculate normals and normalise
    y1=pdf('norm',x,muLogAlphaL,sigmaLogAlphaL);y1=y1/max(y1);y2=pdf('norm',x,muLogAlphaM,sigmaLogAlphaL);y2=y2/max(y2);
    y3=pdf('norm',x,muLogAlphaU,sigmaLogAlphaL);y3=y3/max(y3);y4=pdf('norm',x,muLogAlphaL,sigmaLogAlphaU);y4=y4/max(y4);
    y5=pdf('norm',x,muLogAlphaM,sigmaLogAlphaU);y5=y5/max(y5);y6=pdf('norm',x,muLogAlphaU,sigmaLogAlphaU);y6=y6/max(y6);
    y7=exp(y1);y7=y7/max(y7);y8=exp(y2);y8=y8/max(y8);y9=exp(y3);y9=y9/max(y9);y10=exp(y4);y10=y10/max(y10);y11=exp(y5);y11=y11/max(y11);y12=exp(y6);y12=y12/max(y12);
    subplot(nrows,ncols,3),plot(x,[y1;y2;y3;y4;y5;y6]),xlim([-4.5,2]);ntitle('log Alpha' )
    subplot(nrows,ncols,4),plot(exp(x),[y7;y8;y9;y10;y11;y12]),xlim([0,1.5]);ntitle('Alpha' )
    
    %Lambdas calculate normals and normalise
    y1=pdf('norm',x,muLogLambdaL,sigmaLogLambdaL);y1=y1/max(y1);y2=pdf('norm',x,muLogLambdaM,sigmaLogLambdaL);y2=y2/max(y2);
    y3=pdf('norm',x,muLogLambdaU,sigmaLogLambdaL);y3=y3/max(y3);y4=pdf('norm',x,muLogLambdaL,sigmaLogLambdaU);y4=y4/max(y4);
    y5=pdf('norm',x,muLogLambdaM,sigmaLogLambdaU);y5=y5/max(y5);y6=pdf('norm',x,muLogLambdaU,sigmaLogLambdaU);y6=y6/max(y6);
    y7=exp(y1);y7=y7/max(y7);y8=exp(y2);y8=y8/max(y8);y9=exp(y3);y9=y9/max(y9);y10=exp(y4);y10=y10/max(y10);y11=exp(y5);y11=y11/max(y11);y12=exp(y6);y12=y12/max(y12);
    subplot(nrows,ncols,5),plot(x,[y1;y2;y3;y4;y5;y6]),xlim([-1.5,3]);ntitle('log Lambda')
    subplot(nrows,ncols,6),plot(exp(x),[y7;y8;y9;y10;y11;y12]),xlim([0,7]);ntitle('Lambda'),figtitle(figName);
    
    %etas
    y1=pdf('norm',x,muEtaL,sigmaEtaL);y1=y1/max(y1);y2=pdf('norm',x,muEtaM,sigmaEtaL);y2=y2/max(y2);
    y3=pdf('norm',x,muEtaU,sigmaEtaL);y3=y3/max(y3);y4=pdf('norm',x,muEtaL,sigmaEtaU);y4=y4/max(y4);
    y5=pdf('norm',x,muEtaM,sigmaEtaU);y5=y5/max(y5);y6=pdf('norm',x,muEtaU,sigmaEtaU);y6=y6/max(y6);
    y7=exp(y1);y7=y7/max(y7);y8=exp(y2);y8=y8/max(y8);y9=exp(y3);y9=y9/max(y9);y10=exp(y4);y10=y10/max(y10);y11=exp(y5);y11=y11/max(y11);y12=exp(y6);y12=y12/max(y12);
    subplot(nrows,ncols,8),plot(x,[y1;y2;y3;y4;y5;y6]),xlim([-7,7]);ntitle('Eta')
    
    saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');
end

%% Print information for user
disp('**************');
disp(['running model#_',num2str(runModelNum),':'])
disp([modelName,'_started:_',datestr(clock)])
disp(['MCMC number_',num2str(whichJAGS)])
disp(['running on_',dataSource])
disp(['with_',priorName])
disp('**************');

%% Initialise matrices
%initialise matrices with nan values of size subjects x conditions x trials
choice = nan(nSubjects,nConditions,maxTrials); %initialise choice data matrix 
dx1 = nan(nSubjects,nConditions,maxTrials); dx2 = dx1; dx3 = dx1; dx4=dx1;%initialise changes in wealth
out1= dx1;out2=dx1;out3=dx1;out4=dx1;%initialise outcomes
deuLin= nan(nSubjects,nConditions,maxTrials);deuLog= nan(nSubjects,nConditions,maxTrials);%initialise changes in linear utility / log utility

%% Compile choice & gamble data
% Jags cannot deal with partial observations, so we need to specify gamble info for all nodes. This doesn't change anything.
for i = 1:nSubjects
    for c = 1:nConditions
        switch c %condition
            case 1% add
                trialInds=1:length(Choice_add{subjList(i)});%generate indices for each trial
                choice(i,c,trialInds)=Choice_add{subjList(i)}(trialInds);%assign to temporary variables
                dx1(i,c,trialInds)=LinU_Gam1_1_add{subjList(i)}(trialInds);%assign changes in wealth dx for outcome 1
                dx2(i,c,trialInds)=LinU_Gam1_2_add{subjList(i)}(trialInds);%same for outcome 2 etc.
                dx3(i,c,trialInds)=LinU_Gam2_1_add{subjList(i)}(trialInds);
                dx4(i,c,trialInds)=LinU_Gam2_2_add{subjList(i)}(trialInds);
                out1(i,c,trialInds)=dx1(i,c,trialInds);out2(i,c,trialInds)=dx2(i,c,trialInds);%specify as outcomes 1 to 4
                out3(i,c,trialInds)=dx3(i,c,trialInds);out4(i,c,trialInds)=dx4(i,c,trialInds);
                deuLin(i,c,trialInds)=delta_EU_Lin_add{subjList(i)}(trialInds);%specify changes in expected utility for each gamble for linear utility
                deuLog(i,c,trialInds)=delta_EU_Log_add{subjList(i)}(trialInds);%specify changes in expected utility for each gamble for log utility
            case 2% multi
                trialInds=1:length(Choice_multi{subjList(i)});
                choice(i,c,trialInds)=Choice_multi{subjList(i)}(trialInds);
                dx1(i,c,trialInds)=LinU_Gam1_1_multi{subjList(i)}(trialInds);
                dx2(i,c,trialInds)=LinU_Gam1_2_multi{subjList(i)}(trialInds);
                dx3(i,c,trialInds)=LinU_Gam2_1_multi{subjList(i)}(trialInds);
                dx4(i,c,trialInds)=LinU_Gam2_2_multi{subjList(i)}(trialInds);
                out1(i,c,trialInds)=LogU_Gam1_1_multi{subjList(i)}(trialInds);out2(i,c,trialInds)=LogU_Gam1_2_multi{subjList(i)}(trialInds);
                out3(i,c,trialInds)=LogU_Gam2_1_multi{subjList(i)}(trialInds);out4(i,c,trialInds)=LogU_Gam2_2_multi{subjList(i)}(trialInds);
                deuLin(i,c,trialInds)=delta_EU_Lin_multi{subjList(i)}(trialInds);
                deuLog(i,c,trialInds)=delta_EU_Log_multi{subjList(i)}(trialInds);
        end
    end
end

%% Compute Time average growth for choices
% since outcomes are specified in terms of linear utility for additive
% condition and log utility for multiplicative condition, then the average
% for each gamble corresponds to time average additive growth per trial, or
% time average multiplicative growth per trial.
chosen1=choice.*out1;chosen2=choice.*out2;chosen3=abs(choice-1).*out3;chosen4=abs(choice-1).*out4;%outcome of the chosen gamble
%chosenX is a nSubjects x nConditions x nTrials matrix, 0 if not chosen,
%growth factor if chosen, nan if trial not presented to subject for each of
%the four stimuli
chosenTop=chosen1+chosen3;chosenBot=chosen2+chosen4;%top or bottom outcome of the chosen gamble
taChosen=(chosenTop+chosenBot)./2;taChosen=nanmean(taChosen,3);%time average of chosen gamble,
%mean over all trials for every subject for both conditions (nSubjects x
%nConditions)
csvwrite([startDir,'/data/data_TimeAveragesChosenGambles'],taChosen);%save in data directory for further analysis later


%% Compute trajectories of wealth according to choices for each agent
%randomise as coinflip between top and bottom
%multiply the binaries with top,and the not-binary with bottom and sum the
%two to get the actual outcome. partial sum for additive condition, and
%product for mult.
% randi([0 1], n,m)
% size(ChosenTop)

%% Compute no brainer choice frequencies for statewise dominance
if runPlots
    figure;
    subplot(3,1,1);y=[];for i =1:19;y=[y,mean(NoBrainerChoiceCorrect_multi{i})];end;bar(y);xlabel('Subjects');ylabel('Choice proportion add')
    subplot(3,1,2),x=[];for i =1:19,x=[x,mean(NoBrainerChoiceCorrect_add{i})];end,bar(x);ylabel('Choice proportion mult')
    subplot(3,1,3),z=(x+y)/2;bar(z);ylabel('Choice proportion both')
    figtitle('No brainer choice proportions')
    csvwrite('data_proportionCorrectNoBrainers',[x',y',z']);
end

%% Plot isoelastic utility exemplars
if runPlots
    figure;
    figName='isoelastic utility functions';
    figtitle(figName);
    dx=-500:100:1000;%change in wealth due to outcome
    w0=1005;%initial wealth
    w1=w0+dx;%new wealth at next time step
    u0Lin=w0;u1Lin=w1;duLin=u1Lin-u0Lin;
    u0Log=log(w0);u1Log=log(w1);duLog=u1Log-u0Log;
    eta=-1;%risk seeking
    u0Iso=((w0.^(1-eta))-1)/(1-eta);u1Iso=((w1.^(1-eta))-1)/(1-eta);duIso=u1Iso-u0Iso;
    eta=2;%sub log risk aversion
    u0Iso2=((w0.^(1-eta))-1)/(1-eta);u1Iso2=((w1.^(1-eta))-1)/(1-eta);duIso2=u1Iso2-u0Iso2;
    plot(dx,duIso/max(duIso),'g',dx,duLin/max(duLin),'b',dx,duLog/max(duLog),'r',dx,duIso2/max(duIso2),'m')
    xlabel('change in wealth due to outcome');
    ylabel('utility');
    legend({'eta -1','eta 0 (LinU)', 'eta 1 (LogU)', 'eta 2'},'Location','SouthEast');
    axis('square');
    saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');
end

%% Compile wealth data
wealths=[Wealth_add;Wealth_multi];
wealths=wealths(:,subjList);%limits to subjects in subjlist

%% Truncate data
%truncate down to 312 trials
truncInds=1:nTrials;
choice=choice(:,:,truncInds);
dx1=dx1(:,:,truncInds);dx2=dx2(:,:,truncInds);dx3=dx3(:,:,truncInds);dx4=dx4(:,:,truncInds);
deuLin=deuLin(:,:,truncInds);deuLog=deuLog(:,:,truncInds);

%% Add gamble data for missing trials
%sub 2 multi had only 299 trials, therefore we add random gambles to pad out to 312. this
%allows jags to work since doesn't work for partial observation. this does not affect
% parameter estimation. nans in the choice data are allowed as long as all covariates are not nan.
if nConditions==2 && any(subjList==2)
    dx1(2,2,300:312)=dx1(2,2,1:13);dx2(2,2,300:312)=dx2(2,2,1:13);
    dx3(2,2,300:312)=dx3(2,2,1:13);dx4(2,2,300:312)=dx4(2,2,1:13);
    deuLin(2,2,300:312)=deuLin(2,2,1:13);
    deuLog(2,2,300:312)=deuLog(2,2,1:13);
end

%% Nan check
disp([num2str(length(find(isnan(choice)))),'_nans in choice data']);%nans in choice data do not matter
disp([num2str(length(find(isnan(dx1)))),'_nans in gambles 1 matrix'])% nans in gamble matrices do, since model will not run
disp([num2str(length(find(isnan(dx2)))),'_nans in gambles 2 matrix'])
disp([num2str(length(find(isnan(dx3)))),'_nans in gambles 3 matrix'])
disp([num2str(length(find(isnan(dx4)))),'_nans in gambles 4 matrix'])
disp([num2str(length(find(isnan(deuLin)))),'_nans in deu_lin'])
disp([num2str(length(find(isnan(deuLog)))),'_nans in deu_log'])

%% Visualise changes in wealth
if runPlots
    figure;subplot(2,1,1),histogram([(dx1(:,1,:)+dx2(:,1,:))/2,(dx3(:,1,:)+dx4(:,1,:))/2],'Normalization','pdf')
    hold on,histogram([(dx1(:,2,:)+dx2(:,2,:))/2,(dx3(:,2,:)+dx4(:,2,:))/2],'Normalization','pdf');
    ntitle('EV of dx - blue: add, red: mult')
    subplot(2,1,2)
    tmp1=(out1(:,1,:)+out2(:,1,:))/2;tmp2=(out3(:,1,:)+out4(:,1,:))/2;data1=[tmp1(:);tmp2(:)];
    tmp1=(out1(:,2,:)+out2(:,2,:))/2;tmp2=(out3(:,2,:)+out4(:,2,:))/2;data2=[tmp1(:);tmp2(:)];
    nBins1=15;nBins2=15;
    hAx1 = gca;% create one axis or get current axes
    posAx1 = get(hAx1, 'Position');% get position of axis
    hAx2 = axes('Position', posAx1);% create an overlapping axis at the same location
    h1=histogram(hAx1,data1,nBins1,'Normalization','probability');% histogram for first data vector
    h2=histogram(hAx2,data2,nBins2,'Normalization','probability','FaceColor','r');
    set(hAx2,'Color','none');% make second axis transparent
    ylabel(hAx1,'Probability');% ylabel for histogram 1
    set(hAx2,'XAxisLocation','top');
    xlabel(hAx2,'Time average multiplicative growth');xlabel(hAx1,'Time average additive growth');set(hAx1,'XColor','r');set(hAx2,'XColor','b')
end

%% Compute Time averages for eta beta space
if runPlots
    nEtas=100;nBetas=100;TimeAvAdd=[];TimeAvMult=[];%number of etas and betas to simulate; initialise
    etaParam=linspace(muEtaL,muEtaU,nEtas);betaParam=linspace(exp(muLogBetaL),exp(muLogBetaU),nBetas);
    for e=1:nEtas
        for b=1:nBetas
            [TimeAvAdd(e,b),TimeAvMult(e,b)]=computeEtaBeta2TimeAv(etaParam(e),betaParam(b),wealths,dx1,dx2,dx3,dx4);
        end
    end
    [etaVals,betaVals]=meshgrid(betaParam,etaParam);
    figure,subplot(2,2,1),surf(etaVals,betaVals,TimeAvMult),title('mult');ylabel('beta'),xlabel('eta'),zlabel('Time av');
    subplot(2,2,2),surf(etaVals,betaVals,TimeAvAdd),title('add');ylabel('beta'),xlabel('eta'),zlabel('Time av');
    subplot(2,2,3),yyaxis right,plot(etaParam,TimeAvMult(:,1),'r'),ylabel('Time Average Multiplicative Growth Rate')
    yyaxis left,plot(etaParam,TimeAvAdd(:,1),'b'),xlabel('eta'),ylabel('Time Average Additive Growth Rate')
    line([0, 0], ylim, 'LineWidth', 1.5, 'Color', 'b','LineStyle','--'),line([1, 1], ylim, 'LineWidth', 1.5, 'Color', 'r','LineStyle','--')
    subplot(2,2,4),yyaxis right,plot(etaParam,(TimeAvMult(:,nBetas)),'r'),ylabel('Time Average Multiplicative Growth Rate')
    yyaxis left,plot(etaParam,(TimeAvAdd(:,nBetas)),'b'),xlabel('eta'),ylabel('Time Average Additive Growth Rate')
    line([0, 0], ylim, 'LineWidth', 1.5, 'Color', 'b','LineStyle','--'),line([1, 1], ylim, 'LineWidth', 1.5, 'Color', 'r','LineStyle','--')
    save('data_etaBeta2Time','TimeAvAdd','TimeAvMult','-v7.3') 
end

%% Configure data structure for graphical model & parameters to monitor
%everything you want jags to use
switch runModelNum
    case {1} %Parameter estimation of eta
        dataStruct = struct(...
            'wealths',wealths,'nSubjects', nSubjects,'nConditions',nConditions,...
            'nTrials',nTrials,'dx1',dx1,'dx2',dx2,'dx3',dx3,'dx4',dx4,'y',choice,...
            'muLogBetaL',muLogBetaL,'muLogBetaU',muLogBetaU,'sigmaLogBetaL',sigmaLogBetaL,'sigmaLogBetaU',sigmaLogBetaU,...
            'muEtaL',muEtaL,'muEtaU',muEtaU,'sigmaEtaL',sigmaEtaL,'sigmaEtaU',sigmaEtaU);
        
    case {2}
        dataStruct = struct(...
            'wealths',wealths,'nSubjects', nSubjects,'nConditions',nConditions,...
            'nTrials',nTrials,'dx1',dx1,'dx2',dx2,'dx3',dx3,'dx4',dx4,...
            'y',choice,'muLogBetaL',muLogBetaL,...
            'muLogBetaU',muLogBetaU,'sigmaLogBetaL',sigmaLogBetaL,'sigmaLogBetaU',sigmaLogBetaU,...
            'muLogAlphaL',muLogAlphaL,'muLogAlphaU',muLogAlphaU,'sigmaLogAlphaL',sigmaLogAlphaL,...
            'sigmaLogAlphaU',sigmaLogAlphaU,'muLogLambdaL',muLogLambdaL,'muLogLambdaU',muLogLambdaU,...
            'sigmaLogLambdaL',sigmaLogLambdaL,'sigmaLogLambdaU',sigmaLogLambdaU,'muEtaL',muEtaL,...
            'muEtaU',muEtaU,'sigmaEtaL',sigmaEtaL,'sigmaEtaU',sigmaEtaU,...
            'pz',pz);
end

for i = 1:nChains
    switch runModelNum
        
        case {1}  %Parameter estimation
            monitorParameters = {'beta_iso','eta','mu_eta','sigma_eta','mu_log_beta_iso','sigma_log_beta_iso'};
            S=struct; init0(i)=S; %sets initial values as empty so randomly seeded
            
        case {2} %Model selection for subjectwise strong models
            
            monitorParameters = {...
                'eta_iso','eta_tw','alphaGain','alphaLoss','lambda',...%utility params
                'beta_tw','beta_pt','beta_iso',...%betas
                'z','px_z1','px_z2','delta_z1','sum_z'};%model indicator
            S=struct; init0(i)=S;   %sets initial values as empty so randomly seeded
            
    end
end

%% Run JAGS sampling via matJAGS
tic;fprintf( 'Running JAGS ...\n' ); % start clock to time % display

[samples, stats] = matjags( ...
    dataStruct, ...                           % Observed data
    fullfile(jagsDir, [modelName '.txt']), ...% File that contains model definition
    init0, ...                                % Initial values for latent variables
    whichJAGS,...                             % Specifies which copy of JAGS to run on
    'doparallel' , doParallel, ...            % Parallelization flag
    'nchains', nChains,...                    % Number of MCMC chains
    'nburnin', nBurnin,...                    % Number of burnin steps
    'nsamples', nSamples, ...                 % Number of samples to extract
    'thin', nThin, ...                        % Thinning parameter
    'dic', doDIC, ...                         % Do the DIC?
    'monitorparams', monitorParameters, ...   % List of latent variables to monitor
    'savejagsoutput' , 1 , ...                % Save command line output produced by JAGS?
    'verbosity' , 1 , ...                     % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 1 ,...                        % clean up of temporary files?
    'rndseed',1);                             % Randomise seed; 0=no; 1=yes

toc % end clock

%% Save stats and samples
disp('saving samples and stats...')
save(['samples_stats/' modelName,'_',priorName,'_',dataSource,'_burn_',num2str(nBurnin),'_samps_',num2str(nSamples),'_chains_',num2str(nChains),'_',datestr(clock)],'stats','samples','-v7.3')

%% Print readouts
disp('stats:'),disp(stats)%print out structure of stats output
disp('samples:'),disp(samples);%print out structure of samples output
try
    rhats=fields(stats.Rhat);
    for lp = 1: length(rhats)
        disp(['stats.Rhat.',rhats{lp}]);
        eval(strcat('stats.Rhat.',rhats{lp}))
    end
catch
    disp('no field for stats.Rhat')
end
