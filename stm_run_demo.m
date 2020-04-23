function [] = stm_run_demo()
% STM: y = xbeta + Zu + epsilon, u: NT by 1 spatio-temporal AR1(phi) kronecker CAR(gamma)
rng('default');  rng(8);

% input data: Y, X, Z, W
% Y is the N*T*J column vector of response variable: stacked first by
% location, then by time, then by group replicates
% X is the corresponding matrix of covariates
% Z is the corespoonding vector with spatio-temporal index (1 to N*T) 
% W is N*N spatial adjancency matrix (binary, 1 if two locations are adjacent)
load('stm_run_demo_data.mat', 'Y','X','Z','W')

% set environmental variables (ev)
ev.verbose = 0;  
ev.saveAsMat = false;   ev.saveAsMatFileName = 'stmOut.mat';  % save the MCMC output
ev.nonspatial = 0;  ev.nontemporal = 0;  ev.nonrandom = 0;  % run the full spatio-temporal mixed model
ev.betaprior = 1; %1=flat prior on the data domain; 2=spike&slab prior on wavelet transformed data
% ev.niter = 6e3;  ev.burnin = 5e3;   %MCMC setting  
ev.niter = 60;  ev.burnin = 50;   %short run for this demo
ev.computeDIC = 1;

% fit STM
out = stm(Y, X, Z, W, ev); 

% posterior summary of parameters
summaryPara = prctile(out.matPara,[2.5 25 50 75 97.5]);  %Percentiles
format short;   disp(summaryPara)

% goodness-of-fit measure for Bayesian model
meanEs = mean(out.Es, 1); 
E1 = meanEs(1); E2 = meanEs(2); DIC = -4*E1 + 2*E2;  %DIC measure
% format long g;  disp([-2*E1, -2*E1+2*E2, DIC])
 fprintf('Deviance = %3.2f, pD = %3.2f, DIC = %3.2f\n\n', [-2*E1, -2*E1+2*E2, DIC])

end
