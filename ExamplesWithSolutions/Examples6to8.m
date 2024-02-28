%% The script to run examples 6,7 and 8 of Arts etal. (2024)
clear all
clc

%% Example 6 Condition Based Maintenance:
Cu = 6000
Cp = 500
L = 10
lambda = 0.311850311840798 % from excel sheet MLE
pmf = poisspdf(0:L,lambda);
epsilon = 10^(-6);
[ M , gstar , P0, P1] = CBMdp(pmf,Cu,Cp,epsilon) %all hard work is in this function
CostAgeBased = 49.29195464 % from excel sheet optimization
SavingCBMvsUBM = (CostAgeBased-gstar)/CostAgeBased


%% Example 7 estimate hyperparameters with MLE

% Data from paper (x_{i,\max} and t_{i,\max})
tmax = [24.8; 27.2; 44.2; 13.5; 58.0];
xmax = [10; 10; 10; 10; 10];
% Set up the log-likelihood
LLfun = @(params) -sum(log( nbinpdf(xmax,params(1),params(2)./(params(2)+tmax))));
% Optimize the log-likelihood
param0 = [5,2] %some initial guess for the solver
paramhat = fmincon(LLfun,param0,[],[],[],[],[0.001,0.001],[Inf,Inf]);
OLL = -LLfun(paramhat) %optimized log-likelihood

%% Example 8 POMDP
Cu = 6000;
Cp = 500;
L = 10;
a = 7.409679083029343; % from code above (=paramhat(1))
b = 21.691893050859001; % from code above (=paramhat(2))
[CostRateBayesian, Tbayesian] = cbmBayesPoisson(Cu,Cp,L,a,b) %heavy lifting happens in this function
SavingLearning = (gstar - CostRateBayesian)/CostRateBayesian