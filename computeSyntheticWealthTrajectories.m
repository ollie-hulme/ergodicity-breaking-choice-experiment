%% computeSyntheticWealthTrajectories
% This script computes wealth trajectories for synthetic agents playing the ergodicity game
% repeatedly over several timescales, hours to weeks to years.
%
% Danish Research Centre for Magnetic Resonance, www.drcmr.dk
% Oliver Hulme,2018-2019

%% Clean up & set mod
clc;clear;close all; tic;%clean up and set timer

%% Set modes and specify figures to plot
testmode=0; %flags whether to run in test mode or not for debugging
plotFigs=8:13;%determines which figures to plot and save

%% Set Prospect Theory variables
nExponents=3;nLambdas=3;nPTagents=nExponents*nLambdas; %#ok<NASGU> %number of parameters & agents
exponentGain=linspace(0.3,0.9,nExponents);lossLambda=linspace(1,3,nLambdas);[exponentGain,lossLambda]=meshgrid(exponentGain,lossLambda); %initialise parameters
exponentGain=exponentGain(:);exponentLoss=exponentGain;%for purposes of this simulation since redundant in models with lambda
lossLambda=lossLambda(:);

%% Set time variables
s_pertrial=9.5;% average length of gamble trial in active session
t_perhr=round((60*60)/s_pertrial); t_perday=t_perhr*24;t_perweek=t_perday*7; t_peryr=t_perday*365; %trials per unit time
nOutcomes=8;nUniqueGamblesPairs=144;
cycles_peryear=ceil(t_peryr/nUniqueGamblesPairs);%number of complete cycles of unique gambles in 1 yr
nYears=1;nCycles=cycles_peryear*nYears;nTrials=t_perweek;%

%% Set other variables
nCond=2;%number of conditions to simulate
nAgents=nPTagents+3;%n PT settings, + 3 (linear and log and time criterion)
wealth=1000;%sets wealth to the endowment wealth, choices should be independent of wealth
synAgent(1).add=[];synAgent(1).mult=[];%set up struct to contain all synthetic agents

%% Randomly permute trials without replacement
disp('permuting trials:');tic
randTrials = arrayfun(@(x)randperm(nUniqueGamblesPairs),(1:nCycles)','UniformOutput',0); %#ok<NASGU>
randTrials = cell2mat(randTrials);toc; % generates a list of permutations of unique gamble pairs trials without replacement over n cycles
if testmode==0
    nTrials=t_peryr*nYears;%only when not in test mode does number of Trials correspond to a year
end
randTrials=randTrials(1:nTrials);%shortens number of trials according to nTrials
disp('***')

%% Flip coins
%if outside of agent loop then same coin outcome sequence for all
coinFlip=round(rand(1,nTrials));

%% Column names for big matrix of variables
colOut_1=1:2;% outcomes for gamble 1
colOut_2=3:4;% ditto gamble 2
colDelU_1=5:6;% change in utility for each outcome of gamble 1
colDelU_2=7:8;% ditto gamble 2
colExDelU_1=9;% expected change in utility gamble 1
colExDelU_2=10;% ditto gamble 2
colDelExDelU=11;% difference in expected change in utility between gamble 1 and gamble 2 (pos when 1 > 2)
colChosenGamb=12;% the gamble chosen by the agent
colChosenOut=13:14;% the outcomes of the chosen gamble
colCoin=15;% the coin outcome heads or tails
colRealiseOut=16;% the outcome realised by the coin toss
colCumWealth=17;% the wealth cumulated at that trial
colCP=18;%the choice probabilities
nCols=18;%number of columns specified

%% Outcome space
middleValue=ceil((nOutcomes+1)/2);%calculate middle fractal number e.g. 5 if 9 fractals
outcomes=exp(linspace(log(0.4467),log(2.2387),nOutcomes+1));% calculate growth factors outcomes
outcomes=outcomes(outcomes~=outcomes(middleValue));%keep all except the middle growth factor
[Xm,Ym]=meshgrid(outcomes(outcomes<0.999),outcomes(outcomes>1.001));
gambles(:,:,1)= [Xm(:),Ym(:)];%set multiplicative outcome space
outcomes=linspace(-428,428,nOutcomes+1);%overwrite variable with additive outcomes, growth increments
outcomes=outcomes(outcomes~=outcomes(middleValue));%keep all except the middle increment
[Xa,Ya]=meshgrid(outcomes(outcomes<-0.001),outcomes(outcomes>0.001));
gambles(:,:,2)= [Xa(:),Ya(:)];%set additive outcome space

%% Set gamble space
nGambles=(length(outcomes)/2)^2;%its divided by two since these are mixed gambles only
nGamblePairs=nGambles^2;%this is all pairwise combinations of outcomes
[W,Z]=meshgrid(1:nGambles,1:nGambles);%16 by 16 gamble space
gamble_pairs(:,:,1)=[gambles(W(:),:,1),gambles(Z(:),:,1)];%1 indicates for multiplicative condition
gamble_pairs(:,:,2)=[gambles(W(:),:,2),gambles(Z(:),:,2)];%1 indicates for additive condition

%% Filter non-unique pairs of gambles
disp('filtering unique gambles:');tic; %start timer
x=nan(1,nGamblePairs);%pre allocate temporary variable x
for lp=1:nGamblePairs
    x(lp)=size(unique(gamble_pairs(lp,:,1)),2);%find gamble pairs with no shared outcomes
end
uniqueIndices=find(x==4);uniqueGambles=gamble_pairs(uniqueIndices,:,:);%space of gamble pairs with unique outcomes, this has two sheets, one for each dynmamic
toc;disp('***');%finish timer

%% Compute time and ensemble averages of growth factors
%key: time vs. ens-emble, GF growth factor, GI growth increment, m multi, a additive.
w=[10^3,10^6,10^9];%three wealths of different orders of magnitude
timeGFm=repmat((Xm.*Ym).^0.5,1,1,3);%geometric mean (same as the exp (E (log growth factors)))% Mult: wealth invariant, thus doesnt change as t tend to inf.
timeGFa=(computeIncrement2Factor(Xa,w).*computeIncrement2Factor(Ya,w)).^0.5;%geomentric mean %  Add: wealth variant, asymptotes to 1 as t tends to inf.
ensGFm=repmat((Xm+Ym)/2,1,1,3);%ensemble mean of growth factors, wealth invariant  % Mult: wealth invariant, thus doesnt change as t tends to inf.
ensGFa=(computeIncrement2Factor(Xa,w)+computeIncrement2Factor(Ya,w))/2;%wealth variant   % Add: wealth variant, asymptotes to 1 as t tends to inf.
timeGIm=(computeFactor2Increment(Xm,w)+computeFactor2Increment(Ym,w))/2;%    % Mult: wealth variant, asymptotes to zero for losses and inf for gains
timeGIa=repmat((Xa+Ya)/2,1,1,3);%wealth invariant  % Add: wealth invariant, thus doest change as t tends to inf
ensGIm=(computeFactor2Increment(Xm,w)+computeFactor2Increment(Ym,w))/2;%wealth variantMult: wealth variant. Asymptotes to zero for losses, and inf for gains
ensGIa=repmat((Xa+Ya)/2,1,1,3);%wealth invariant
timeEGm=repmat((log(Xm)+log(Ym))/2,1,1,3);%wealth invariant
timeEGa=(log(computeIncrement2Factor(Xa,w))+log(computeIncrement2Factor(Ya,w)))/2;%wealth variant %  Add: Wealth variant, asymptotes to zero as t tend to inf
ensEGm=repmat((log(Xm)+log(Ym))/2,1,1,3);%wealth invariant
ensEGa=(log(computeIncrement2Factor(Xa,w))+log(computeIncrement2Factor(Ya,w)))/2; % Add: wealth variant, for additive dynamics it assymptotes to zero.

%% Fig. 1. Plot all time and ensemble averages
if find(plotFigs==1)>0
    figure,
    Ct=1;%counter
    
    %m
    subplot(8,6,(0*6)+Ct),bar3(timeGFm(:,:,1));zlim([0, 1.6]);set(gca,"yticklabel",[1.2, 1.5,1.8,2.2 ],"xticklabel",[0.4,0.5,0.7,0.8]);ntitle('timeGFm lo w')
    subplot(8,6,(1*6)+Ct),bar3(timeGFm(:,:,2));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('timeGFm med w')
    subplot(8,6,(2*6)+Ct),bar3(timeGFm(:,:,3));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('timeGFm hi w')
    subplot(8,6,(3*6)+Ct),bar3(timeGFm(:,:,1));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('timeGFm long t limit')%time average of mult is wealth invariant as time tends to inf
    
    %a
    subplot(8,6,(4*6)+Ct),bar3(timeGFa(:,:,1));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('timeGFa');
    subplot(8,6,(5*6)+Ct),bar3(timeGFa(:,:,2));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(6*6)+Ct),bar3(timeGFa(:,:,3));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(7*6)+Ct),bar3(ones(4));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);%time average of add asymptotes to 1 as time tends to inf
    Ct=Ct+1;
    
    %m
    subplot(8,6,(0*6)+Ct),bar3(ensGFm(:,:,1));zlim([0, 1.6]);set(gca,"yticklabel",[1.2, 1.5,1.8,2.2 ],"xticklabel",[0.4,0.5,0.7,0.8]);ntitle('ensGFm');
    subplot(8,6,(1*6)+Ct),bar3(ensGFm(:,:,2));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(2*6)+Ct),bar3(ensGFm(:,:,3));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(3*6)+Ct),bar3(ensGFm(:,:,1));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);%time average of mult is wealth invariant as t tends to inf
    
    %a
    subplot(8,6,(4*6)+Ct),bar3(ensGFa(:,:,1));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('ensGFa');
    subplot(8,6,(5*6)+Ct),bar3(ensGFa(:,:,2));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(6*6)+Ct),bar3(ensGFa(:,:,3));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(7*6)+Ct),bar3(ones(4));zlim([0, 1.6]);set(gca,"yticklabel",[],"xticklabel",[]);%time average of add asymptotes to 1 as time tends to inf
    Ct=Ct+1;
    
    %m
    subplot(8,6,(0*6)+Ct),bar3(timeGIm(:,:,1));set(gca,"yticklabel",[1.2, 1.5,1.8,2.2 ],"xticklabel",[0.4,0.5,0.7,0.8]);ntitle('timeGIm');
    subplot(8,6,(1*6)+Ct),bar3(timeGIm(:,:,2));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(2*6)+Ct),bar3(timeGIm(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(3*6)+Ct),bar3(timeGIm(:,:,1)>0);zlim([0,0.91]);set(gca,"yticklabel",[],"xticklabel",[]);%goes to zero or inf -this is a hack to crop the bars to represent infinity
    
    %a
    subplot(8,6,(4*6)+Ct),bar3(timeGIa(:,:,1));set(gca,"yticklabel",[],"xticklabel",[]);ntitle('timeGIa');
    subplot(8,6,(5*6)+Ct),bar3(timeGIa(:,:,2));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(6*6)+Ct),bar3(timeGIa(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(7*6)+Ct),bar3(timeGIa(:,:,1));set(gca,"yticklabel",[],"xticklabel",[]);
    Ct=Ct+1;
    
    %m
    subplot(8,6,(0*6)+Ct),bar3(ensGIm(:,:,1));set(gca,"yticklabel",[1.2, 1.5,1.8,2.2 ],"xticklabel",[0.4,0.5,0.7,0.8]);ntitle('ensGIm');
    subplot(8,6,(1*6)+Ct),bar3(ensGIm(:,:,2));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(2*6)+Ct),bar3(ensGIm(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(3*6)+Ct),bar3(ensGIm(:,:,1)>0);zlim([0,0.91]);set(gca,"yticklabel",[],"xticklabel",[]);%time average- goes to zero or inf -this is a hack to crop the bars to represent infinity
    
    %a
    subplot(8,6,(4*6)+Ct),bar3(ensGIa(:,:,1));set(gca,"yticklabel",[],"xticklabel",[]);ntitle('ensGIa');
    subplot(8,6,(5*6)+Ct),bar3(ensGIa(:,:,2));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(6*6)+Ct),bar3(ensGIa(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(7*6)+Ct),bar3(ensGIa(:,:,1));set(gca,"yticklabel",[],"xticklabel",[]);%time average
    Ct=Ct+1;
    
    %m
    subplot(8,6,(0*6)+Ct),bar3(timeEGm(:,:,1));set(gca,"yticklabel",[1.2, 1.5,1.8,2.2 ],"xticklabel",[0.4,0.5,0.7,0.8]);ntitle('ensEGm');%this is the geometric mean of the growth increments
    subplot(8,6,(1*6)+Ct),bar3(timeEGm(:,:,2));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(2*6)+Ct),bar3(timeEGm(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(3*6)+Ct),bar3(timeEGm(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);%time average
    
    %a
    subplot(8,6,(4*6)+Ct),bar3(timeEGa(:,:,1));set(gca,"yticklabel",[],"xticklabel",[]);ntitle('timeEGa');%this is the geometric mean of the growth increments
    subplot(8,6,(5*6)+Ct),bar3(timeEGa(:,:,2));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(6*6)+Ct),bar3(timeEGa(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(7*6)+Ct),bar3(zeros(size(timeEGa(:,:,3))));set(gca,"yticklabel",[],"xticklabel",[]);%time average
    Ct=Ct+1;
    
    %m
    subplot(8,6,(0*6)+Ct),bar3(ensEGm(:,:,1));set(gca,"yticklabel",[1.2, 1.5,1.8,2.2 ],"xticklabel",[0.4,0.5,0.7,0.8]);ntitle('ensEGm');%this is the geometric mean of the growth increments
    subplot(8,6,(1*6)+Ct),bar3(ensEGm(:,:,2));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(2*6)+Ct),bar3(ensEGm(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(3*6)+Ct),bar3(ensEGm(:,:,1));set(gca,"yticklabel",[],"xticklabel",[]);%time average
    
    %a
    subplot(8,6,(4*6)++Ct),bar3(ensEGa(:,:,1));set(gca,"yticklabel",[],"xticklabel",[]);ntitle('ensEGa');%this is the geometric mean of the growth increments
    subplot(8,6,(5*6)+Ct),bar3(ensEGa(:,:,2));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(6*6)+Ct),bar3(ensEGa(:,:,3));set(gca,"yticklabel",[],"xticklabel",[]);
    subplot(8,6,(7*6)+Ct),bar3(zeros(size(timeEGa(:,:,3))));set(gca,"yticklabel",[],"xticklabel",[]);%time average
    
    figtitle('Growth rates and utilities over gamble space');%title for whole figure, has to be called after subplot
end

%% Fig. 2  Plot time averages only
if find(plotFigs==2)>0
    figure
    % time average additive growth
    subplot(4,1,1), bar3(timeGIa(:,:,1));zlim([-200, 200]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('time av additive growth for add');
    subplot(4,1,2), bar3((timeGIm(:,:,1)>0)+0.01);zlim([0,1]);zlim([-1, 1]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('time av additive growth for mult');
    % time average exponential growth
    subplot(4,1,3), bar3(zeros(size(timeEGa(:,:,3))));zlim([-0.2, 0.2]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('time av exponential growth for add');%time average
    subplot(4,1,4), bar3(timeEGm(:,:,3));zlim([-0.31, 0.31]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('time av exponential growth for mult');
    figtitle('Time average growths');
    saveas(gcf,'timeAverageGrowthsBar','epsc');
end

%% Fig. 3. Expected change in linear utility
% This is the same as calculating an ensemble average of wealth
% change. As such it is identical to ensemble and time average
% growth increments. Only for additive dynamics does time average
% linear utility equal ensemble average, and is wealth invariant.
% For multiplicative dynamics, linear utility asymptotes to 0 for
% losses and inf for gains
if find(plotFigs==3)>0
    figure; Ct=1;
    evDelLinUm=ensGIm;%same as ensemble average of growth increments
    evDelLinUa=ensGIa;%same as ensemble average of growth increments
    subplot(8,4,(0*4)+Ct),bar3(evDelLinUa(:,:,1));ntitle('Linear utility lo w add')
    subplot(8,4,(1*4)+Ct),bar3(evDelLinUa(:,:,2));ntitle('Linear utility med w add')
    subplot(8,4,(2*4)+Ct),bar3(evDelLinUa(:,:,3));ntitle('Linear utility hi w add')
    subplot(8,4,(3*4)+Ct),bar3(evDelLinUa(:,:,1));ntitle('Linear utility long time limit add')%limit of time average
    subplot(8,4,(4*4)+Ct),bar3(evDelLinUm(:,:,1));ntitle('Linear utility lo w mult')
    subplot(8,4,(5*4)+Ct),bar3(evDelLinUm(:,:,2));ntitle('Linear utility med w mult')
    subplot(8,4,(6*4)+Ct),bar3(evDelLinUm(:,:,3));ntitle('Linear utility high w mult')
    subplot(8,4,(7*4)+Ct),bar3(evDelLinUm(:,:,1)>0);zlim([0,0.91]);set(gca,"yticklabel",[],"xticklabel",[]);ntitle('Linear utility long time limit mult')%time average- goes to zero or inf -this is a hack to crop the bars to represent infinity
    Ct=Ct+1;
    
    %% Expected change in log utility
    % This is the same as the time and ensemble average exponential growth.
    % Only for multiplicative dynamics does time average log utility equal its ensemble average, and is wealth invariant.
    % For additive dynamics, log utility asymptotes to 0 for gains, and losses
    evDelLogUm=ensEGm;%same as time average exponential growth rate
    evDelLogUa=ensEGa;%same as ensemble average exponential growth rate
    subplot(8,4,(0*4)+Ct),bar3(evDelLogUa(:,:,1));ntitle('Log utility')%this is the geometric mean of the growth increments
    subplot(8,4,(1*4)+Ct),bar3(evDelLogUa(:,:,2));
    subplot(8,4,(2*4)+Ct),bar3(evDelLogUa(:,:,3));
    subplot(8,4,(3*4)+Ct),bar3(zeros(4,4));%limit
    subplot(8,4,(4*4)+Ct),bar3(evDelLogUm(:,:,1));%this is the geometric mean of the growth increments
    subplot(8,4,(5*4)+Ct),bar3(evDelLogUm(:,:,2));
    subplot(8,4,(6*4)+Ct),bar3(evDelLogUm(:,:,3));
    subplot(8,4,(7*4)+Ct),bar3(evDelLogUm(:,:,3));%limit
    Ct=Ct+1;
    
    %% Expected change in PT utility
    % PT is wealth invariant for additive amounts since it doesnt
    % depend on wealth. For mult dynamics the wealth tends to either
    % zero or inf, so PT likewise goes to zero or inf.
    w=[10^3,10^9,10^12];%three wealths of different orders of magnitude
    %             u1m=computeUtility_PT(computeFactor2Increment(Xm,w),0.88,0.88,2.25,w);
    %             u2m=computeUtility_PT(computeFactor2Increment(Ym,w),0.88,0.88,2.25,w);
    u1m=computeExpectedUtility({computeFactor2Increment(Xm,w)},0.88,0.88,2.25,0.5);
    u2m=computeExpectedUtility({computeFactor2Increment(Ym,w)},0.88,0.88,2.25,0.5);
    evDelPTUm=(u1m+u2m)/2;
    %             u1a=computeUtility_PT(repmat(Xa,1,1,3),0.88,0.88,2.25,w);
    %             u2a=computeUtility_PT(repmat(Ya,1,1,3),0.88,0.88,2.25,w);
    u1a=computeExpectedUtility({repmat(Xa,1,1,3)},0.88,0.88,2.25,0.5);
    u2a=computeExpectedUtility({repmat(Ya,1,1,3)},0.88,0.88,2.25,0.5);
    evDelPTUa=(u1a+u2a)/2;
    subplot(8,4,(0*4)+Ct),bar3(evDelPTUa(:,:,1));ntitle('Prospect theory')
    subplot(8,4,(1*4)+Ct),bar3(evDelPTUa(:,:,2));
    subplot(8,4,(2*4)+Ct),bar3(evDelPTUa(:,:,3));
    subplot(8,4,(3*4)+Ct),bar3(evDelPTUa(:,:,1));%limit
    subplot(8,4,(4*4)+Ct),bar3(evDelPTUm(:,:,1));
    subplot(8,4,(5*4)+Ct),bar3(evDelPTUm(:,:,2));
    subplot(8,4,(6*4)+Ct),bar3(evDelPTUm(:,:,3));
    subplot(8,4,(7*4)+Ct),bar3(evDelLinUm(:,:,1)>0);zlim([0,0.91]);set(gca,"yticklabel",[],"xticklabel",[]);%time average- goes to zero or inf -this is a hack to crop the bars to represent infinity
    Ct=Ct+1;
    
    %% Expected change in time criterion utility
    %linear utility for additive, logarithmic utility for
    %multiplicative
    subplot(8,4,(0*4)+Ct),bar3(evDelLinUa(:,:,1));ntitle('Time Criterion')
    subplot(8,4,(1*4)+Ct),bar3(evDelLinUa(:,:,2));
    subplot(8,4,(2*4)+Ct),bar3(evDelLinUa(:,:,3));
    subplot(8,4,(3*4)+Ct),bar3(evDelLinUa(:,:,1));%limit of time average
    subplot(8,4,(4*4)+Ct),bar3(evDelLogUm(:,:,1));%this is the geometric mean of the growth increments
    subplot(8,4,(5*4)+Ct),bar3(evDelLogUm(:,:,2));
    subplot(8,4,(6*4)+Ct),bar3(evDelLogUm(:,:,3));
    subplot(8,4,(7*4)+Ct),bar3(evDelLogUm(:,:,3));%limit
    figtitle('Utilities at different wealths and long time limits')
    saveas(gcf,'timeAverageUtilitiesBar','epsc');
end

%% Fig. 4. Plot limits
if find(plotFigs==4)>0
    figure;ct=1;
    subplot(2,4,ct),bar3(evDelLinUa(:,:,1));%limit linear a
    ct=ct+1;
    subplot(2,4,ct),bar3(zeros(4,4));%limit log a
    ct=ct+1;
    subplot(2,4,ct),bar3(evDelPTUa(:,:,1));%limit pt a
    ct=ct+1;
    subplot(2,4,ct),bar3(evDelLinUa(:,:,1));%limit time crit a
    ct=ct+1;
    subplot(2,4,ct),bar3(evDelLinUm(:,:,1)>0);zlim([0,0.91]);set(gca,"yticklabel",[],"xticklabel",[]);%limit linear m
    ct=ct+1;
    subplot(2,4,ct),bar3(evDelLogUm(:,:,3));%limit log m
    ct=ct+1;
    subplot(2,4,ct),bar3(evDelLinUm(:,:,1)>0);zlim([0,0.91]);set(gca,"yticklabel",[],"xticklabel",[]);%limit pt m
    ct=ct+1;
    subplot(2,4,ct),bar3(evDelLogUm(:,:,3));%limit time crit m
    ct=ct+1;
    saveas(gcf,'timeAverageUtilities_timeLimits_Bar','epsc');
end

%% Utility agent loop
%Here we loop over different synthetic agents with different
%utility functions
for ulp=1:nAgents;disp(['agent_',num2str(ulp),':started']);tic %utility function loop%
    
    % Condition loop
    for condlp=1:nCond;disp(['condition loop',num2str(condlp),':started']);tic %condition loop, 1=mult, 2=add
        
        gambles_year=uniqueGambles(randTrials,:,condlp);% Compile outcomes and utilities over a year
        m=nan(nTrials,nCols);
        m(:,1:4)=uniqueGambles(randTrials,:,condlp);
        
        %initialise wealth
        if condlp==1
            wealth=log(1000);%initialise starting wealth for each agent
        elseif condlp==2
            wealth=1000;
        end
        
        % Trial loop
        for trlp=1:nTrials
            
            % Compute utilities according to dynamic
            switch ulp
                case 1 %linear utility
                    if condlp==1
                        m(trlp,[colDelU_1,colDelU_2])=computeFactor2Increment(m(trlp,[colOut_1,colOut_2]),wealth);%the change in utilities for each of the 4 outcomes
                    elseif condlp==2
                        m(trlp,[colDelU_1,colDelU_2])=m(trlp,[colOut_1,colOut_2]);%the change in utilities for each of the 4 outcomes
                    end
                case 2 %log utility
                    if condlp==1
                        m(trlp,[colDelU_1,colDelU_2])= log(m(trlp,[colOut_1,colOut_2]));
                    elseif condlp==2
                        m(trlp,[colDelU_1,colDelU_2])=log(computeIncrement2Factor(m(trlp,[colOut_1,colOut_2]),wealth));%the change in utilities for each of the 4 outcomes
                    end
                case {3,4,5,6,7,8,9,10,11}%PT
                    if condlp==1
                        incr=computeFactor2Increment(m(trlp,[colOut_1,colOut_2]),wealth);%compute increments dependend on wealth
                        %                                 [m(trlp,[colDelU_1,colDelU_2])]=computeUtility_PT(incr,exponentGain(ulp-2),exponentGain(ulp-2),lossLambda(ulp-2),wealth);
                        [m(trlp,[colDelU_1,colDelU_2])]=computeExpectedUtility({incr},exponentGain(ulp-2),exponentGain(ulp-2),lossLambda(ulp-2),0.5);
                    elseif condlp==2
                        %                                 [m(trlp,[colDelU_1,colDelU_2])]=computeUtility_PT(m(trlp,[colOut_1,colOut_2]),exponentGain(ulp-2),exponentGain(ulp-2),lossLambda(ulp-2),wealth);
                        [m(trlp,[colDelU_1,colDelU_2])]=computeExpectedUtility({m(trlp,[colOut_1,colOut_2])},exponentGain(ulp-2),exponentGain(ulp-2),lossLambda(ulp-2),0.5);
                    end
                case 12 %time criterion
                    if condlp==1%if mult then log
                        m(trlp,[colDelU_1,colDelU_2])= log(m(trlp,[colOut_1,colOut_2]));
                    elseif condlp==2%if add then linear
                        m(trlp,[colDelU_1,colDelU_2])=m(trlp,[colOut_1,colOut_2]);%the change in utilities for each of the 4 outcomes
                    end
            end %switch ulp
            
            %expected utilities of each gamble
            [m(trlp,colExDelU_1),m(trlp,colExDelU_2)]=computeExpectationValue(m(trlp,[colDelU_1,colDelU_2]));%the expectation value of change of utilities for each of 2 gambles(EDU)
            
            %differences in expected utilities between gambles
            m(trlp,colDelExDelU)=m(trlp,colExDelU_1)-m(trlp,colExDelU_2); %difference in expectation value of change in utility, between gamble 1 and 2, and is pos where gamble1 > gamble2
            
            %                     %choose gambles
            %                     beta=20;%sets sensitivity parameter to some arbitrary value
            %                     m(trlp,colCP)=1/(1+exp(-((m(trlp,colDelExDelU))*beta)));
            
            m(trlp,colChosenGamb)=sign(m(trlp,colDelExDelU));%+1 if gamble 1 chosen (and thus higher EDU), -1 if gamble 2
            if m(trlp,colChosenGamb)==1
                m(trlp,colChosenOut)=m(trlp,colOut_1);%chose gamble 1 when higher utility
            elseif m(trlp,colChosenGamb)==-1
                m(trlp,colChosenOut)=m(trlp,colOut_2);%and likewise for gamble 2
            elseif m(trlp,colChosenGamb)==0 %randomly choose if utilities are identical
                randChoice=round(rand(1));
                if randChoice==1
                    m(trlp,colChosenOut)=m(trlp,colOut_1);
                else
                    m(trlp,colChosenOut)=m(trlp,colOut_2);
                end
            end
            
            %Record coin flip outcomes (same for all agents)
            m(trlp,colCoin)=coinFlip(trlp);%determines which of two outcomes is realised
            
            %Realises outcomes for chosen gamble
            if m(trlp,colCoin)==1
                m(trlp,colRealiseOut)=m(trlp,colChosenOut(1));
            else
                m(trlp,colRealiseOut)=m(trlp,colChosenOut(2));
            end
            
            %Apply to wealth
            if condlp==1%mult
                wealth=wealth+log(m(trlp,colRealiseOut));%otherwise numbers exceed numerical limits
            elseif condlp==2%add
                wealth=wealth+m(trlp,colRealiseOut);
            end
            m(trlp,colCumWealth)=wealth;
        end% trial loop
        
        %% Save big matrix into condition structure
        if condlp==1
            condition.mult=m;
        elseif condlp==2
            condition.add=m;
        end
        m=[];%clear matrix
        disp(['condition loop',num2str(condlp),':ended']);toc %condition loop
    end
    
    %% Save data into agent structure
    synAgent(ulp)=condition; %#ok<SAGROW>
    disp(['agent_',num2str(ulp),':ended']);toc %utility loop
    disp('***')
end

%% Configure utility functions
n_dx=1000;
w=1000;%endowed wealth
mxW=2000;%w*2.2387;%maximum wealth after 1 outcome
mnW=0;%w*0.4467;%min wealth after 1 outcome
W=linspace(mnW,mxW,n_dx);%all possible end wealths
dX=W-w;%all possible delta wealths
duLin=dX/max(dX);%normalise
duLog=log(W)-log(w);
duLog=duLog/max(duLog);%normalise
duPro=nan(nPTagents,n_dx);
for ptlp=1:nPTagents
    duPro(ptlp,:)=computeExpectedUtility({dX},exponentGain(ptlp),exponentLoss(ptlp),lossLambda(ptlp),0.5);
    duPro(ptlp,:)=duPro(ptlp,:)/max(duPro(ptlp,:));%normalise
end
du=[duLin;duLog;duPro];%matrix containing all three utility functions

%% Concatenate additive and multiplicative

%         concatWealth=[repmat(wealth,1,3);nan(nTrials-1,3)];
%         for cyclp=0:(cycles_peryear-1)
%             a=concatWealth(~isnan(concatWealth(:,1)));
%             trials=(cyclp*288)+1:(cyclp+1)*288;%trials to pull out for concatenation into wealth trajectory of additive and multiplicative
%             cumWealthAdd=cumsum([concatWealth(trials(1),:);Sim.add(trials,36:38)]);%cumulative sum for growth increments during additive session
%             cumWealthMult=cumprod(Sim.mult(trials,36:38)).*cumWealthAdd(trials(1),:);%cumulative product for growth factors during multiplicative session
%             cumWealth=[cumWealthAdd;cumWealthMult];%concatenate cumulative wealths for both sessions
%             concatWealth(trials(1):trials(1)+576,:)=cumWealth;%aggregate into larger vector for 1 year trajectory
%         end

%% Plot pre-processing
timeVec=(1:nTrials)*9.5;%set up vector of times for each trial in real experimental time in seconds
cumWealthAdd=nan(nAgents,nTrials);cumWealthMult=cumWealthAdd;
for aglp=1:nAgents
    cumWealthAdd(aglp,:)=synAgent(aglp).add(:,colCumWealth);
    cumWealthMult(aglp,:)=synAgent(aglp).mult(:,colCumWealth);
end

%% Fig. 5. Plot utility functions of different agents
if find(plotFigs==5)>0
    figure, plot(dX,du(1,:),'r',dX,du(2,:),'b',dX,du(3,:),'m','LineWidth',1.5);ylim([-1,1]);xlim([-1000,1000]);axis square;
    saveas(gcf,'utilityFunctions3models','epsc');
end
if find(plotFigs==6)>0 | find(plotFigs==7)>0
    figure, h=plot(dX,du(3:11,:),'LineWidth',1.5);ylim([-1,1]);xlim([-1000,1000]);axis square
    saveas(gcf,'utilityFunctionsPTmodels','epsc');
    c=get(h,'Color');
    colorInds=reshape([c{1:9,:}],3,9)';
    sizeInds=repmat(120,1,9);
    figure,scatter(exponentGain',lossLambda',sizeInds,colorInds,'s','filled'),xlim([0,1]),ylim([0,3.25]);axis square;
    saveas(gcf,'parameterSpacePT','epsc');
end

%% Fig. 6. Plot wealth for 3 model classes
if find(plotFigs==6)>0
    figure % compare synthetic agents
    threeModels=[1,2,3];%linear, log and the eighth PT which is close to the canonical [0.88 2.21] of K&T original paper
    try %#ok<TRYNC>
        subplot(4,2,1); h5=plot(timeVec(1:t_perhr),  cumWealthAdd(threeModels,1:t_perhr));
        c=get(h5,'Color');
        subplot(4,2,3), h6=plot(timeVec(1:t_perday), cumWealthAdd(threeModels,1:t_perday));
        subplot(4,2,5), h7=plot(timeVec(1:t_perweek),cumWealthAdd(threeModels,1:t_perweek));
        subplot(4,2,7),h8=plot(timeVec(1:t_peryr),  cumWealthAdd(threeModels,1:t_peryr));
    end
    try
        subplot(4,2,2);  h1=plot(timeVec(1:t_perhr), cumWealthMult(threeModels,1:t_perhr));legend
        c=get(h1,'Color');
        subplot(4,2,4); h2=plot(timeVec(1:t_perday), cumWealthMult(threeModels,1:t_perday));
        subplot(4,2,6); h3=plot(timeVec(1:t_perweek),cumWealthMult(threeModels,1:t_perweek));
        subplot(4,2,8); h4=plot(timeVec(1:t_peryr),  cumWealthMult(threeModels,1:t_peryr));
    catch
        disp('not enough trials computed for all plots')
    end
    saveas(gcf,['model class wealth trajectories_',datestr(clock)],'epsc')
end

%% Fig. 7. Plot wealth for all agents over all time scales
if find(plotFigs==7)>0
    figure % compare synthetic agents
    try
        subplot(4,4,2); h1=plot(timeVec(1:t_perhr),  cumWealthMult(:,1:t_perhr));legend;c=get(h1,'Color');
        subplot(4,4,6); h2=plot(timeVec(1:t_perday), cumWealthMult(:,1:t_perday));subplot(4,4,10); h3=plot(timeVec(1:t_perweek),cumWealthMult(:,1:t_perweek));
        subplot(4,4,14); h4=plot(timeVec(1:t_peryr), cumWealthMult(:,1:t_peryr));
    catch,disp('not enough trials computed for all plots')
    end
    
    try
        subplot(4,4,1); h5=plot(timeVec(1:t_perhr),  cumWealthAdd(:,1:t_perhr));c=get(h5,'Color');
        subplot(4,4,5), h6=plot(timeVec(1:t_perday), cumWealthAdd(:,1:t_perday));subplot(4,4,9), h7=plot(timeVec(1:t_perweek),cumWealthAdd(:,1:t_perweek));
        subplot(4,4,13),h8=plot(timeVec(1:t_peryr),  cumWealthAdd(:,1:t_peryr));
    end
    saveas(gcf,['wealth trajectories_',datestr(clock)],'epsc')
end

%% Fig. 8. Plot wealth over week for multiplicative for all agents
if find(plotFigs==8)>0
    figure
    h9=plot(timeVec(1:t_perweek),cumWealthMult(1:11,1:t_perweek));hold on
    plot(timeVec(1:t_perweek),cumWealthMult(12,1:t_perweek),'LineWidth',2)
    saveas(gcf,['wealth trajectories_week_mult_',datestr(clock)],'epsc')
end

%% Fig. 9. Plot wealth over week for additive for all agents
if find(plotFigs==9)>0
    figure
    plot(timeVec(1:t_perweek),cumWealthAdd(1:11,1:t_perweek)),hold on
    plot(timeVec(1:t_perweek),cumWealthAdd(12,1:t_perweek),'LineWidth',3)
    saveas(gcf,['wealth trajectories_week_add_',datestr(clock)],'epsc')
end

%% Fig. 10. Plot wealth over yr for multiplicative for all agents
if find(plotFigs==10)>0
    figure
    plot(timeVec(1:t_peryr),cumWealthMult(1:11,1:t_peryr)),hold on
    plot(timeVec(1:t_peryr),cumWealthMult(12,1:t_peryr),'LineWidth',3)
    saveas(gcf,['wealth trajectories_yr_mult_',datestr(clock)],'epsc')
end

%% Fig. 11. Plot wealth over yr for additive for all agents
if find(plotFigs==11)>0
    figure
    plot(timeVec(1:t_peryr),cumWealthAdd(1:11,1:t_peryr)),hold on
    plot(timeVec(1:t_peryr),cumWealthAdd(12,1:t_peryr),'LineWidth',3)
    saveas(gcf,['wealth trajectories_yr_add_',datestr(clock)],'epsc')
end

%% Fig. 12. Plot wealth over day for multiplicative for all agents
if find(plotFigs==12)>0
    figure
    plot(timeVec(1:t_perday),cumWealthMult(1:11,1:t_perday)),hold on
    plot(timeVec(1:t_perday),cumWealthMult(12,1:t_perday),'LineWidth',3)
    saveas(gcf,['wealth trajectories_day_mult_',datestr(clock)],'epsc')
end

%% Fig. 13. Plot wealth over day for additive for all agents
if find(plotFigs==13)>0
    figure
    plot(timeVec(1:t_perday),cumWealthAdd(1:11,1:t_perday)),hold on
    plot(timeVec(1:t_perday),cumWealthAdd(12,1:t_perday),'LineWidth',3)
    saveas(gcf,['wealth trajectories_day_add_',datestr(clock)],'epsc')
end

%% Fig. 14. Rank agents and scatter
if find(plotFigs==14)>0
    [~,p] = sort(cumWealthMult(:,nTrials),'descend');
    rankMult = 1:nAgents;
    rankMult(p) = rankMult;
    [~,p] = sort(cumWealthAdd(:,nTrials),'descend');
    rankAdd = 1:nAgents;
    rankAdd(p) = rankAdd;
    colorInds=reshape([c{1:12,:}],3,12)';
    figure
    scatter(rankMult,rankAdd,10,colorInds,'s','filled');legend
    saveas(gcf,['agentRanks',datestr(clock)],'epsc');
end

%% Fig. 13. Plot finite time average growth rates
if find(plotFigs==13)>0
    finTimeAvAddGrowth=(cumWealthAdd(:,end)-wealth)/nTrials;
    finTimeAvMultGrowth=(log(cumWealthMult(:,end)/wealth))/nTrials;
    finTimeAvAddGrowthNorm=finTimeAvAddGrowth/(max(finTimeAvAddGrowth));
    finTimeAvMultGrowthNorm=finTimeAvMultGrowth/(max(finTimeAvMultGrowth));
    generalFinTimeAvGrowth=finTimeAvAddGrowthNorm+finTimeAvMultGrowthNorm;
    figure,subplot(3,1,1),bar(finTimeAvAddGrowth)
    subplot(3,1,2),bar(finTimeAvMultGrowth)
    subplot(3,1,3),bar([generalFinTimeAvGrowth;finTimeAvAddGrowthNorm(1)+finTimeAvMultGrowthNorm(2)])
end

%% Compute duration of script
disp('entire script:ended');toc