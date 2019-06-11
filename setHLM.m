function setHLM(runModelNum,synthMode,whichJAGS,whichQuals,runPlots,doParallel)
%% setHLM
% setHLM sets up multiple HLM models to run sequentially according to inputs
% This function takes the following inputs:
% runModelNum - which models to run
% synthMode   - which data to use: 1 real participants, 2 synthetic agents for model recovery, 3 synthetic agents for param recovery
% whichJAGS   - which copy of matjags to run on. this allows parallel jobs to run as long as they use different matjags
% whichQuals  - sets the order of qualities to run
% runPlots    - plot results, illustrations, schematics
% 
% There are three qualities for several variables, each selected by whichQuals
% qualities  are 'bronze','silver','gold'
% gold is highest quality but takes longest, bronzest lowest but fastest
% etc.

%% Specifies qualities to be selected from
numRuns      = length(whichQuals);%how many separate instances of an MCMC to run
nBurnin      = [1e2,1e3,1e4,2e4,4e4];%from 100 to 40k
nSamples     = [5e1,5e2,5e3,1e4,2e4];%from 50 to 20k
nChains      = [4,4,4,4,4];%
nThin        = 10;%thinnning factor, 1 = no thinning, 2=every 2nd etc.

%% Specifies subjects
subjList{1}  = [1:4,6:19];%the subjects from the experiment (excludes 5 who didnt learn the stimuli)
subjList{2}  = 1:27;%synthetic agents for model recovery
subjList{3}  = 1:9;%synthetic agents for parameter recover, 9 equally spaced agents forming a 3x3 grid in eta space
subjList{4}  = 1:19;subjList{5}=1:19;subjList{6}=1:19;subjList{7}=1:19;subjList{8}=1:19;%synthetic agents for parameter recover
subjList{9}  = 1:19;subjList{10}=1:19;subjList{11}=1:19;subjList{12}=1:19;%synthetic agents for parameter recover
subjList{13}  = 1:19;subjList{14}=1:19;subjList{15}=1:19;subjList{16}=1:19;%synthetic agents for parameter recover
subjList{17}  = 1:19;subjList{18}=1:19;subjList{19}=1:19;subjList{20}=1:19;%synthetic agents for parameter recover
subjList{21}  = 1:19;subjList{22}=1:19;subjList{23}=1:19;subjList{24}=1:19;%synthetic agents for parameter recover
subjList{25}  = 1:19;subjList{26}=1:19;subjList{27}=1:19;subjList{28}=1:19;%synthetic agents for parameter recover

%% Runs HLMs sequentiallt
for i=1:numRuns
    computeHLM(runModelNum,nBurnin(whichQuals(i)),nSamples(whichQuals(i)),nThin,nChains(whichQuals(i)),subjList{synthMode},whichJAGS,runPlots,synthMode,doParallel)
end