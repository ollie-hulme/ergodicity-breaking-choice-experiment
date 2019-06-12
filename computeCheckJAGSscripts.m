%% computeCheckJAGSscripts
% computeCheckJAGSscripts is used for debugging JAGS scripts. JAGS scripts
% cannot be debugged in a sequential manner as in matlab, thus this script 
% does sanity check computations whether variables computed by the JAGS 
% script are meaningful

%% Load jags output 
clc
cd 'samples_stats'
disp('load the jags output')
uiopen %opens interface to load relevant samples

%% Check PT algebra
%         lamb1[i,c,t]       = ifelse(dx1[i,c,t]>0, 1, lambda[i,c]) #set lambda to 1 for positive outcomes, and to lambda param for negative
%         lamb2[i,c,t]       = ifelse(dx2[i,c,t]>0, 1, lambda[i,c])
%         lamb3[i,c,t]       = ifelse(dx3[i,c,t]>0, 1, lambda[i,c])
%         lamb4[i,c,t]       = ifelse(dx4[i,c,t]>0, 1, lambda[i,c])
%         eadx1[i,c,t]       = pow(adx1[i,c,t],alpha[i,c]) #exponentiate absolute value of outcome by alpha       
%         eadx2[i,c,t]       = pow(adx2[i,c,t],alpha[i,c])
%         eadx3[i,c,t]       = pow(adx3[i,c,t],alpha[i,c])
%         eadx4[i,c,t]       = pow(adx4[i,c,t],alpha[i,c])
%         pu1[i,c,t]         = lamb1[i,c,t] *  eadx1[i,c,t] #multiply by lambda variable
%         pu2[i,c,t]         = lamb2[i,c,t] *  eadx2[i,c,t]
%         pu3[i,c,t]         = lamb3[i,c,t] *  eadx3[i,c,t]
%         pu4[i,c,t]         = lamb4[i,c,t] *  eadx4[i,c,t]
%         epug1[i,c,t]       = (pu1[i,c,t]+pu2[i,c,t])/2  #calculate mean utility for gamble
%         epug2[i,c,t]       = (pu3[i,c,t]+pu4[i,c,t])/2           
%         deupt[i,c,t]       = epug1[i,c,t]-epug2[i,c,t] #difference in expected util        
%         sdeupt[i,c,t]      = -1 * beta_pt[i,c] * deupt[i,c,t] # sensitivity-scaled difference in eu
%         tmppt[i,c,t]       = (1)/(1+(exp(sdeupt[i,c,t]))) # choice probability
%         theta[i,c,t,2]     = max(0.0001,min(0.9999,tmppt[i,c,t])) # ensure 0 < cp < 1

%% Set the address of the sample
ch=1;%chain
sm=3;%samples
sb=12;%subjects
cn=1;%condition
t=196;%trial

%% Lambda
disp(['lambda=',num2str(samples.lambda(ch,sm,sb,cn))])
disp(['positive?_',num2str(samples.dx1(ch,sm,sb,cn,t)>0)])
disp(['sets lambda to_',num2str(samples.lamb1(ch,sm,sb,cn,t))])
disp(['positive?_',num2str(samples.dx2(ch,sm,sb,cn,t)>0)])
disp(['sets lambda to_',num2str(samples.lamb2(ch,sm,sb,cn,t))])
disp(['positive?_',num2str(samples.dx3(ch,sm,sb,cn,t)>0)])
disp(['sets lambda to_',num2str(samples.lamb3(ch,sm,sb,cn,t))])
disp(['positive?_',num2str(samples.dx4(ch,sm,sb,cn,t)>0)])
disp(['sets lambda to_',num2str(samples.lamb4(ch,sm,sb,cn,t))])
disp('*****')

%% Exponentiated by alpha absolute dx
disp(['alpha=',num2str(samples.alpha(ch,sm,sb,cn))])
disp(['eadx1=',num2str(samples.eadx1(ch,sm,sb,cn,t))])
disp(['adx1=',num2str(samples.adx1(ch,sm,sb,cn,t))])
disp(['eadx1 calculated as=',num2str(samples.adx1(ch,sm,sb,cn,t)^samples.alpha(ch,sm,sb,cn))])
disp('****')
disp(['eadx2=',num2str(samples.eadx2(ch,sm,sb,cn,t))])
disp(['adx2=',num2str(samples.adx2(ch,sm,sb,cn,t))])
disp(['eadx2 calculated as=',num2str(samples.adx2(ch,sm,sb,cn,t)^samples.alpha(ch,sm,sb,cn))])
disp('****')
disp(['eadx3=',num2str(samples.eadx3(ch,sm,sb,cn,t))])
disp(['adx3=',num2str(samples.adx3(ch,sm,sb,cn,t))])
disp(['eadx3 calculated as=',num2str(samples.adx3(ch,sm,sb,cn,t)^samples.alpha(ch,sm,sb,cn))])
disp('****')
disp(['eadx4=',num2str(samples.eadx4(ch,sm,sb,cn,t))])
disp(['adx4=',num2str(samples.adx4(ch,sm,sb,cn,t))])
disp(['eadx4 calculated as=',num2str(samples.adx4(ch,sm,sb,cn,t)^samples.alpha(ch,sm,sb,cn))])

%% Prospect utility
disp('****')
disp(['pu1=',num2str(samples.pu1(ch,sm,sb,cn,t))])
disp(['pu1 calcd as=',num2str(samples.lamb1(ch,sm,sb,cn,t)*samples.eadx1(ch,sm,sb,cn,t))])
% predicted_pu1=computeUtility_PT(samples.dx1(ch,sm,sb,cn,t),samples.alpha(ch,sm,sb,cn),samples.alpha(ch,sm,sb,cn),samples.lambda(ch,sm,sb,cn),1000)
predicted_pu1=computeExpectedUtility({samples.dx1(ch,sm,sb,cn,t)},samples.alpha(ch,sm,sb,cn),samples.alpha(ch,sm,sb,cn),samples.lambda(ch,sm,sb,cn),0.5)
disp('****')
disp(['pu2=',num2str(samples.pu2(ch,sm,sb,cn,t))])
disp(['pu2 calcd as=',num2str(samples.lamb2(ch,sm,sb,cn,t)*samples.eadx2(ch,sm,sb,cn,t))])
% predicted_pu2=computeUtility_PT(samples.dx2(ch,sm,sb,cn,t),samples.alpha(ch,sm,sb,cn),samples.alpha(ch,sm,sb,cn),samples.lambda(ch,sm,sb,cn),1000)
predicted_pu2=computeExpectedUtility({samples.dx2(ch,sm,sb,cn,t)},samples.alpha(ch,sm,sb,cn),samples.alpha(ch,sm,sb,cn),samples.lambda(ch,sm,sb,cn),0.5)
disp('****')
disp(['pu3=',num2str(samples.pu3(ch,sm,sb,cn,t))])
disp(['pu3 calcd as=',num2str(samples.lamb3(ch,sm,sb,cn,t)*samples.eadx3(ch,sm,sb,cn,t))])
% predicted_pu3=computeUtility_PT(samples.dx3(ch,sm,sb,cn,t),samples.alpha(ch,sm,sb,cn),samples.alpha(ch,sm,sb,cn),samples.lambda(ch,sm,sb,cn),1000)
predicted_pu3=computeExpectedUtility({samples.dx3(ch,sm,sb,cn,t)},samples.alpha(ch,sm,sb,cn),samples.alpha(ch,sm,sb,cn),samples.lambda(ch,sm,sb,cn),0.5)
disp('****')
disp(['pu4=',num2str(samples.pu4(ch,sm,sb,cn,t))])
disp(['pu4 calcd as=',num2str(samples.lamb4(ch,sm,sb,cn,t)*samples.eadx4(ch,sm,sb,cn,t))])
% predicted_pu4=computeUtility_PT(samples.dx4(ch,sm,sb,cn,t),samples.alpha(ch,sm,sb,cn),samples.alpha(ch,sm,sb,cn),samples.lambda(ch,sm,sb,cn),1000)
predicted_pu4=computeExpectedUtility({samples.dx4(ch,sm,sb,cn,t)},samples.alpha(ch,sm,sb,cn),samples.alpha(ch,sm,sb,cn),samples.lambda(ch,sm,sb,cn),0.5)

%% expected prospect utility of gamble and its delta
disp('****')
disp(['epug1=',num2str(samples.epug1(ch,sm,sb,cn,t))]);
disp(['epug1 calcd as=',num2str((samples.pu1(ch,sm,sb,cn,t)+samples.pu2(ch,sm,sb,cn,t))/2)]);
disp(['epug2=',num2str(samples.epug2(ch,sm,sb,cn,t))]);
disp(['epug2 calcd as=',num2str((samples.pu3(ch,sm,sb,cn,t)+samples.pu4(ch,sm,sb,cn,t))/2)]);
disp(['deupt=',num2str(samples.deupt(ch,sm,sb,cn,t))]);
disp(['deupt calcd as=',num2str((samples.epug1(ch,sm,sb,cn,t)-samples.epug2(ch,sm,sb,cn,t)))]);

%% sensitivity scaled difference in expected utility
disp(['sdeupt=',num2str(samples.sdeupt(ch,sm,sb,cn,t))]);
disp(['sdeupt calcd as=',num2str(-1 * samples.beta_pt(ch,sm,sb,cn) * samples.deupt(ch,sm,sb,cn,t))]);%            
disp(['tmppt=',num2str(samples.tmppt(ch,sm,sb,cn,t))]);
disp(['tmppt calcd as=',num2str((1)/(1+exp(samples.sdeupt(ch,sm,sb,cn,t))))])