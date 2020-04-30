function [] = sfc_demo()
% Demo for SFC: Spatio-Functional Clustering with wavelet smoothing

% Y: T by N response variable 
% W: N by N spatial adjacency matrix
load('sfc_demo_data.mat','Y','W')

% set environmenal variables (ev)
ev.simu = true; %If true, will simulate data using initial parameters as the true values
ev.verbose = true; 
ev.niter = 5e4;	ev.burnin = 25e3;  ev.thin = 10;  ev.preRun = 1e4; 
ev.shrinkage = true; %wavelet shrinkage? If true, some gamma will be 0 to achieve dimension reduction
ev.updateAlpha = true; 
ev.nrMin = 1; %minimal number of regions within a cluster, related to the max number of clusters (d)
ev.usepenalty = false; %give penalty to large number of clusters (d)? 
ev.useExp = false; %use exponential penalty?
ev.fixcluster = false;  ev.fixd = false;
ev.seed = 20;   %random seeds: e.g. 20

% set priors
% mean_sigma2 = 10; var_sigma2 = 1e4;  %Prior mean and variance of Inverse Gamma for noise level \sigma^2
mean_lambda = 5; var_lambda = 1e4;  %Prior mean and variance of Inverse Gamma for signal-noise ratio \lambda
ev.a_pi = 1;  ev.b_pi = 1;  %shape parameter of Beta prior for shrinkage level p
% ev.alpha_sig2 = 2+mean_sigma2^2/var_sigma2;   ev.invbeta_sig2 = mean_sigma2*(ev.alpha_sig2-1);
ev.alpha_sig2 = 0;  ev.invbeta_sig2 = 0; %Jeffreys prior
ev.alpha_lambda = 2+mean_lambda^2/var_lambda;   ev.invbeta_lambda = mean_lambda*(ev.alpha_lambda-1);

%************************* fit SFC
out = sfc(Y, W, ev); 

% summarize SFC fit
nsample = length(out.pa);  ds = nan(1,nsample); 
for i = 1:nsample; ds(i) = out.pa{i}.d; end
tabulate(ds)
tab = tabulate(ds);   dmax = tab(tab(:,2)==max(tab(:,2)),1); 

% obtain central cluster
N = numel(out.pa{1}.labs);   nsample = length(out.pa);   
PCount = zeros(N);
for i = 1:nsample; for r = 1:out.pa{i}.d; ind = find(out.pa{i}.labs==r); PCount(ind, ind) = PCount(ind,ind)+1; end; end
PCount = 1-PCount./nsample;
pdistvec = zeros(1, N*(N-1)/2); k = 1;
for i = 1:(N-1); for j = (i+1):N;  pdistvec(k) = PCount(i,j); k = k+1;    end; end
clusterRe = linkage(pdistvec, 'ward');
ClusterIndex = cluster(clusterRe, 'maxclust', dmax);
labsCentral = zeros(1, N); for r = 1:dmax;  labsCentral(ClusterIndex==r) = r; end

if ev.simu  % for simulated data, true clustering and parameters are known. Validate. 
    %----------- cluster measure
    % use other methods such as Kmeans
    [AR,RI,MI,HI] = RandIndex(out.pa_true.labs, labsCentral);
    % disp( [AR,RI,MI,HI] )  %adjust rand, rand, disagreement, agreement
    fit.cl = [AR,RI,MI,HI];
    eval = evalclusters(out.pa_true.data', 'kmeans', 'silhouette', 'klist', 1:10);
    rng('default');  rng(8);
    [labsKmeans, ctrs] = kmeans(out.pa_true.data',  eval.OptimalK); %Kmeans on original data domain
    [AR,RI,MI,HI] = RandIndex(out.pa_true.labs, labsKmeans');
    fit.cl = [fit.cl; [AR,RI,MI,HI] ]; 
    
    %----------- clsuter mean curve estimation
    labs0 = out.pa_true.labs; 
    dtrue = max(labs0);  T = numel(out.pa_true.beta)/dtrue; 
    meanCurveTrue = out.pa_true.Wav' * reshape(out.pa_true.beta, [T,dtrue]); 
    meanCurveKmeans = nan(T, dtrue); 
    tmp = ctrs(labsKmeans,:); 
    for r = 1:dtrue; meanCurveKmeans(:,r) = mean(tmp(labs0==r,:), 1)'; end
    meanCurveSFC = zeros(T, dtrue);
    for r = 1:dtrue
        for i = 1:length(out.pa)
            betas = out.pa_true.Wav' * reshape(out.pa{i}.beta, [T, out.pa{i}.d]);
            tmp = betas(:, out.pa{i}.labs(labs0==r)); 
            meanCurveSFC(:,r) = meanCurveSFC(:,r) + mean(tmp, 2);  %average over sites in the cluster
        end
    end
    meanCurveSFC = meanCurveSFC/length(out.pa); %average over MCMC iterations
    fit.rmse = sqrt([mean((meanCurveTrue-meanCurveSFC).^2, 1); ...
        mean((meanCurveTrue-meanCurveKmeans).^2, 1)]);  %average over time (function) for each cluster
    
    % save for figures
    fit.meanCurveTrue = meanCurveTrue; 
    fit.meanCurveSFC = meanCurveSFC;
    
    fprintf('SFC:  AR=%.3f,  RI=%.3f, MI=%.3f, HI=%.3f, RMSE1=%.3f, RMSE2=%.3f, RMSE3=%.3f, RMSE4=%.3f\n',...
        [fit.cl(1,:), fit.rmse(1,:)])
    fprintf('Kmean:  AR=%.3f,  RI=%.3f, MI=%.3f, HI=%.3f, RMSE1=%.3f, RMSE2=%.3f, RMSE3=%.3f, RMSE4=%.3f\n',...
        [fit.cl(2,:), fit.rmse(2,:)])
    
    save('out.mat', 'fit')
end

end


