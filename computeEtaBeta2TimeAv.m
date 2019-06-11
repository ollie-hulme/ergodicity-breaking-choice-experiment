function [timeAvAdd,timeAvMult] = computeEtaBeta2TimeAv (eta,beta,wealthInit,dx1,dx2,dx3,dx4)
%% computeEtaBeta2TimeAv
% computeEtaBeta2TimeAv computes the time average growth over the
% gambles given an eta and beta value. 

%% change in log wealth and wealths
wealthInit=repmat(wealthInit',1,1,size(dx1,3));

dlogx1=log(wealthInit+dx1)-log(wealthInit);dlogx2=log(wealthInit+dx2)-log(wealthInit);
dlogx3=log(wealthInit+dx3)-log(wealthInit);dlogx4=log(wealthInit+dx4)-log(wealthInit);

wealth1=dx1+wealthInit; wealth2=dx2+wealthInit;
wealth3=dx3+wealthInit;wealth4=dx4+wealthInit;%calculate updated wealths

%% compute utilities
uInit=((wealthInit.^(1-eta))-1)./(1-eta);%utility of initial wealth
u1=   ((wealth1.^(1-eta))-1)./(1-eta); u2=((wealth2.^(1-eta))-1)./(1-eta);%utilities of wealth following outcomes
u3=   ((wealth3.^(1-eta))-1)./(1-eta); u4=((wealth4.^(1-eta))-1)./(1-eta);

%% changes in utilities and their expectation value
du1=u1-uInit;du2=u2-uInit;du3=u3-uInit;du4=u4-uInit;%changes in utility for each outcome
edug1=(du1+du2)./2;edug2=(du3+du4)./2;%expected changes in utility
deu=edug1-edug2;%difference in expected changes in utility

%% logistic function
sdeu = -1.* beta.*deu;%scale by sensitivity parameter and make negative  
theta= 1 / (1+(exp(sdeu)));%choice probabilities

%% time averages of each gamble
g1a=(dx1(:,1,:)+dx2(:,1,:))/2;%add
g2a=(dx3(:,1,:)+dx4(:,1,:))/2;
g1m=(dlogx1(:,2,:)+dlogx2(:,2,:))/2;%mult
g2m=(dlogx3(:,2,:)+dlogx4(:,2,:))/2;

%% double check this below by random sampling single gambles and calculating time avs.
timeAvAdd=(theta(:,1,:).*g1a) + ((1-theta(:,1,:).*g2a));%time average additive growth computed from dx's 
timeAvMult=(theta(:,2,:).*g1m) + ((1-theta(:,2,:).*g2m));%time average multiplicative growth
timeAvAdd=mean(timeAvAdd,3); %mean over all trials
timeAvMult=mean(timeAvMult,3);
timeAvAdd=timeAvAdd(1,:); %take first subject only
timeAvMult=timeAvMult(1,:);

end