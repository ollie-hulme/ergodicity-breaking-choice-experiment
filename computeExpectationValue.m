function [eu1, eu2]=computeExpectationValue(utilities)
%% computeExpectationValue
%takes a list of utilities and computes a mean for each pair

eu1=mean(utilities(:,1:2),2);
eu2=mean(utilities(:,3:4),2);
