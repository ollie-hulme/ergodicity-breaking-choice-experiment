% runHLM provides inputs to setHLM for a given run on the cluster

% It provides the following inputs when calling setHLM:

% runModelNum - which models to run
% synthMode   - which data to use: 1 real participants, 2 synthetic agents for model recovery, 3 synthetic agents for param recovery
% whichJAGS   - which copy of matjags to run on. this allows parallel jobs to run as long as they use different matjags
% whichQuals  - sets the order of qualities to run
% runPlots    - plot results, illustrations, schematics

% The idea is that this is written into by the user, then called by a
% cluster job via the terminal:

% runHLM1.m is called by runHLM_job1.sh
% runHLM2.m is called by runHLM_job2.sh, and so on.

%% Add to path
cd .. ;%move to base directory 
addpath(genpath(pwd));%adds base directory and subfolders to path, important for running shell scripts from terminal

%% Specify variables
runModelNum=6;
synthMode=1;
whichJAGS=1;
whichQuals=3:5;
runPlots=0;
doParallel=0;

%% Call setHLM
setHLM(runModelNum,synthMode,whichJAGS,whichQuals,runPlots,doParallel)