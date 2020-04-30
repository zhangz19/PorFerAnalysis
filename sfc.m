function [out] = sfc(Y, W, ev)
% SFC: Spatio-Functional Clustering with wavelet smoothing
% input Y: T by N response variable 
% input W: N by N spatial adjacency matrix
% input ev: environmental variables for general setting
% (C)  Zhen Zhang (2020)      zhangquake@gmail.com

[ev.T,  ev.N] = size(Y);  
ev.p = 1;  % this version of SFC works for p=1 (intercept only).
ev.NT = ev.N*ev.T;   ev.pT = ev.p*ev.T; 
ev.J = floor(log2(ev.T));  % for wavelet: number of levels of wavelet resolution. usually assume ev.T = power(2, ev.J)
ev.indJ = zeros(1, ev.T);  t = 1; for j = 0:ev.J; for k = 1:round(2^(j-1));   ev.indJ(t) = j;   t = t+1;   end; end;   
% construct wavelet transformation matrix
ev.X = eye(ev.T); Wav = nan(ev.T);  for i1 = 1:ev.T;  Wav(:,i1) = wavedec(ev.X(:,i1), ev.J, 'db1');  end;  ev.Wav = Wav; %Wav'; 
ev.Ywav = Wav*Y;
ev.indJ = repmat(ev.indJ, [1, ev.p]);
ev.D = adjance2distance(W); % distance matrix
ev.Nstar = floor(ev.N/ev.nrMin);  %Nstar: max number of clusters
if ev.updateAlpha
    kappagap = 1e-3;     kappavec = kappagap:kappagap:(1-kappagap);
    ev.kappaconst = log(kappavec) - log(1-(1-kappavec).^ev.Nstar);
    ev.kappavec = kappavec;     ev.kappagap = kappagap;
end
ev.alpha_lambda = ev.alpha_lambda*ones(1, 1+ev.J);  ev.invbeta_lambda = ev.invbeta_lambda*ones(1, 1+ev.J);
ev.a_pi = ev.a_pi*ones(1, 1+ev.J);  ev.b_pi = ev.b_pi*ones(1, 1+ev.J); 

%------------------ simulate data if ev.simu == true
%will overwrite the observed data ev.Ywav (on wavelet domain)
pa_true = []; if ev.simu; [ev.Ywav, pa_true] = simulateData(ev); end 

rng('default');  rng(ev.seed);  % set random seeds

%------------------ set initial values
initialD = 8;   %set 0 if initial search over d=1:10
pa = initialCluster(ev, initialD);  %if second argument = 0, explore a range of numbers of clusters (d=1:10)
nr = histc(pa.labs, 1:pa.d); 
pa.beta = nan(ev.pT, pa.d); pa.lambda = nan(pa.d, 1+ev.J);  SSE = 0;  pa.pis = ones(pa.d, 1+ev.J); 
for r = 1:pa.d
    pa.beta(:,r) = mean(ev.Ywav(:,pa.labs==r),2);
    tmp = 1;  if nr(r)>1; tmp = var(ev.Ywav(:,pa.labs==r),[],2)'; end   %for lambda, initial larger, more data driven
    % tmp = 10*tmp;  %one can initilize with diffuse prior (large lambda)
    pa.lambda(r,:)  = accumarray((ev.indJ+1)', tmp', [], @mean)'; 
    SSE = SSE + sum(sum( (ev.Ywav(:,pa.labs==r) - repmat(pa.beta(:,r), [1,nr(r)])).^2 ));
end
pa.beta = reshape(pa.beta, [ev.T*pa.d,1]);
pa.sigma2 = SSE/ev.NT;
pa.lambda = pa.lambda/pa.sigma2;
pa.kappa = 1e-10;  %code will use kappa for \alpha in the paper. As alpha will be used for recording log acceptance rate

%------------------ run the MCMC, store the results
nsample = (ev.niter-ev.burnin)/ev.thin;
out.pa_true = pa_true; %If ev.simu==true, this is the true parameter
out.pa = cell(1, nsample); out.ML = nan(1, nsample); %save results including marginal likelihood
steps = zeros(1,ev.niter); accepts = zeros(1,ev.niter); alphas = zeros(1,ev.niter); betarates = zeros(1,ev.niter);
xds = zeros(1,ev.niter);  %record number of clusters
preRun_lambda = nan(ev.preRun, ev.J+1);    preRun_pis = nan(ev.preRun, ev.J+1); 
MLconst = gammaln(ev.alpha_sig2+ev.NT/2) - ev.NT/2*log(pi);

[ML, pa] = getMarginalLoglike(pa, ev, false);
tmp = getFullLoglike(pa, ev);
FL = tmp.FL;  MSE = tmp.MSE;
fprintf('Initial:  ML=%.3f,  FL=%.3f,  MSE=%.3f, d=%d\n', [ML+MLconst, FL, MSE, pa.d])

tic
iter_save = 1; 
for i = 1:ev.niter
    %if ev.verbose == 1;  fprintf('%6d', pa.d);  if(~mod(i,20));  fprintf('\n');  end;  end
    if ~ ev.fixcluster && i > ev.preRun
        u0 = rand(1);  if ev.fixd; u0 = unifrnd(0.8, 1); end  %if fix d: shiftStep or switchStep (1/2 chance)
        if u0 <= 0.4 && pa.d > 1 && ~ ev.fixd
            [pa, ~] = mergeStep(pa, ML, ev);  steps(i) = 2; accepts(i) = pa.merge; alphas(i) = pa.alpha;
        elseif ((u0 > 0.4 && u0 <= 0.8 && pa.d < ev.Nstar) || pa.d == 1) && ~ ev.fixd
            [pa, ~] = splitStep(pa, ML, ev);  steps(i) = 1; accepts(i) = pa.growth; alphas(i) = pa.alpha;
        elseif u0 > 0.8 && u0 <= 0.9
            [pa, ~] = shiftStep(pa, ML, ev, W);  steps(i) = 3; accepts(i) = pa.shift; alphas(i) = pa.alpha;
        elseif u0 > 0.9
            [pa, ~] = switchStep(pa, ML, ev);  steps(i) = 4; accepts(i) = pa.switch; alphas(i) = pa.alpha;
        end
    end
    
    [pa] = updateSigma2(pa, ev);
    [pa] = updateBeta(pa, ev);
    [pa] = updateLambdaP(pa, ev);
    if ev.updateAlpha;   [pa] = updateKappa(pa, ev);   end
    
    [ML, pa] = getMarginalLoglike(pa, ev, false);
    tmp = getFullLoglike(pa, ev);
    FL = tmp.FL;  MSE = tmp.MSE; 
    if ~ mod(i, 1e3)
        fprintf('iter=%d,  ML=%.3f,  FL=%.3f,  MSE=%.3f, d=%d, sigma2=%.3f\n', [i, ML+MLconst, FL, MSE, pa.d,pa.sigma2])
    end
    
    if i <= ev.preRun
        preRun_lambda(i,:) = mean(pa.lambda,1);    preRun_pis(i,:) = mean(pa.pis,1); 
    end  %average lambda over clusters
    if i == ev.preRun  %update hyperparameters
        for j = 1:size(preRun_lambda,2)
            tmp = gamfit(1./preRun_lambda(ceil(ev.preRun/2):end,j)); %tmp(2): scale (in the denominator)
            ev.alpha_lambda(j) = tmp(1); ev.invbeta_lambda(j) = 1/tmp(2);  %invbeta_lambda: in the numerator
            tmp = preRun_pis(ceil(ev.preRun/2):end,j); 
            if ev.shrinkage && numel(unique(tmp'))~=1
                tmp = betafit(tmp); 
                ev.a_pi(j) = tmp(1); ev.b_pi(j) = tmp(2); 
            end
        end
    end

    betarates(i) = mean(pa.beta == 0);   xds(i) = pa.d;
    if i > ev.burnin &&  ~ mod(i-ev.burnin, ev.thin)
        out.pa{iter_save} = pa;    out.ML(iter_save) = ML; 
        iter_save = iter_save + 1; 
    end
end

out.runtime = toc/60;  %CPUtime in minutes
if ev.verbose
    fprintf('summary of moves and acceptance rate of samples for use:\n')
    tstep = steps((ev.burnin+1):ev.thin:end);
    tacc = accepts((ev.burnin+1):ev.thin:end);
    tabulate(tstep)
    disp([mean(tacc(tstep==1)), mean(tacc(tstep==2)), mean(tacc(tstep==3)), mean(tacc(tstep==4))])
    fprintf('\n%d iterations are done with elapsed time %.2f minutes.\n', ev.niter, out.runtime)
end
end


%****************** util function for SFC
function [pa, ML1] = shiftStep(pa, ML0, ev, W)
% propose whole beta under new cluster
pa1 = pa; distIndex = 1:ev.N; indic = NaN(1, pa.d);
for r = 1:pa.d
    group = distIndex(W(distIndex == pa.center(r),:)==1);
    for s = 1:pa.d;  group(group==pa.center(s)) = [];  end
    gsize = length(group);   if gsize == 0;  indic(r) = 0;  else indic(r) = 1;  end
end
subCenter = pa.center(indic==1); nsize = length(subCenter); % number of centers with >=1 non-center neighbors
target = randsample(1:length(subCenter),1); target = subCenter(target);
group = find(W(target,:)==1); % neighbors of selected center
for r = 1:pa.d;  group(group==pa.center(r)) = [];   end
gsize = length(group); % number of non-center neighbors
if gsize == 1;  newCenter = group;  else  newCenter = randsample(group,1);  end
pa1.center(pa.center==target) = newCenter;  pa1.labs = ClusterGen(pa1.center, ev); %membership may change
nr = histc(pa1.labs, 1:pa1.d);
if min(nr) < ev.nrMin
    pa.alpha = -99;  pa.shift = 0; ML1 = ML0;
else
    [ML1,pa1] = getMarginalLoglike(pa1, ev, true);
    indic = NaN(1, pa1.d);
    for r = 1:pa1.d
        group = distIndex(W(distIndex==pa1.center(r),:)==1);
        for j = 1:pa1.d;   group(group==pa1.center(j)) = [];   end
        gsize = length(group);
        if gsize == 0
            indic(r) = 0;
        else indic(r) = 1;
        end
    end
    if sum(indic)==0
        alpha = 0; shift = 0;
    else
        subCenter = pa1.center(indic==1); Revnsize = length(subCenter); target = newCenter;
        group = find(W(target,:)==1);
        for r = 1:pa.d
            group(group==pa1.center(r)) = [];
        end
        Revgsize = length(group); PRatio = log(nsize) + log(gsize) - log(Revnsize) - log(Revgsize);
        logratio = ML1 - ML0 + PRatio; alpha = min(0,logratio);
        u = log(rand(1)); shift = 0;
        if u <= alpha
            pa = pa1; shift = 1;
        else
            ML1 = ML0;
        end
    end
    pa.alpha = alpha; pa.shift = shift;
end
end


function [pa, ML1] = switchStep(pa, ML0, ev)
perm = 1:pa.d; ind = randsample(perm,2); perm(ind) = ind([2,1]);
pa1 = pa; pa1.center = pa.center(perm);
pa1.labs = ClusterGen(pa1.center, ev);
nr = histc(pa1.labs, 1:pa1.d);
if(min(nr) < ev.nrMin)
    pa.alpha = -99; pa.switch = 0; ML1 = ML0;
else
    [ML1,pa1] = getMarginalLoglike(pa1, ev, true);
    alpha = min(0, ML1 - ML0);
    u = log(rand(1)); switchs = 0;
    if u <= alpha
        pa = pa1; switchs = 1;
    else
        ML1 = ML0;
    end
    pa.alpha = alpha;
    pa.switch = switchs;
end
end


function [pa, ML1] = splitStep(pa, ML0, ev)
pa1 = pa;   pa1.d = pa.d + 1;
newCenter = randsample(setdiff(1:ev.N, pa.center),1);
pa1.center = nan(1, pa1.d); ind = randsample(pa1.d, 1);  indt = (ind-1)*ev.pT + (1:ev.pT);
pa1.center(ind) = newCenter;   i0 = setdiff(1:pa1.d, ind);
pa1.center(i0) = pa.center;  pa1.labs = ClusterGen(pa1.center, ev); nr = histc(pa1.labs, 1:pa1.d);
pa1.beta = zeros(pa1.d*ev.pT,1);  pa1.beta(setdiff(1:size(pa1.beta,1), indt)) = pa.beta;
pa1.lambda = zeros(pa1.d, ev.J+1);  pa1.pis = ones(pa1.d, ev.J+1); 
k = 1;
for i = 1:pa1.d
    if any(pa.center==pa1.center(i))
        pa1.lambda(i,:) = pa.lambda(k,:);  pa1.pis(i,:) = pa.pis(k,:);  k = k+1;
    end
end
if min(nr) < ev.nrMin 
    pa.alpha = -99; pa.growth = 0; ML1 = ML0;
else
    [ML1, pa1] = getMarginalLoglike(pa1, ev, true);
    if ev.usepenalty
        logratio = ML1 - ML0 - 0.5*ev.pT*log(ev.NT);
    elseif ev.useExp
        logratio = ML1 - ML0 - tan(pi/2*pa1.kappa);
    else
        logratio = ML1 - ML0 + log(1-pa1.kappa);
    end
    alpha = min(logratio,0);  u = log(rand(1)); growth = 0;
    if u <= alpha
        pa = pa1; growth = 1;
    else
        ML1 = ML0;
    end
    pa.alpha = alpha;  pa.growth = growth;
end
end


function [pa, ML1] = mergeStep(pa, ML0, ev)
DelCenter = randsample(pa.center, 1);
pa1 = pa; ind = find(pa.center==DelCenter);
indt = (ind-1)*ev.pT + (1:ev.pT); pa1.center(ind) = []; pa1.labs = ClusterGen(pa1.center, ev);
pa1.pis(ind,:) = [];  pa1.lambda(ind,:) = [];
pa1.d = length(pa1.center); pa1.beta = pa.beta(setdiff(1:size(pa.beta,1), indt));
nr = histc(pa1.labs, 1:pa1.d);
if(min(nr) < ev.nrMin)
    pa.alpha = -99; pa.merge = 0; ML1 = ML0;
else
    [ML1, pa1] = getMarginalLoglike(pa1, ev, true);
    if ev.usepenalty
        logratio = ML1 - ML0 + 0.5*ev.pT*log(ev.NT);
    elseif ev.useExp
        logratio = ML1 - ML0 + tan(pi/2*pa1.kappa);
    else
        logratio = ML1 - ML0 - log(1-pa1.kappa);
    end
    alpha = min(logratio,0);
    u = log(rand(1)); merge = 0;
    if u <= alpha;  pa = pa1; merge = 1; else  ML1 = ML0;  end
    pa.alpha = alpha;
    pa.merge = merge;
end
end


function [pa] = updateBeta(pa, ev)
% This version works only for p = 1
for r = 1:pa.d
    lambda = pa.lambda(r, ev.indJ+1);
    indr = find(pa.labs==r); nr = length(indr); rT = (r-1)*ev.pT + (1:ev.pT); betavec = pa.beta(rT);
    if ev.shrinkage
        for i = 1:ev.p
            t = 1;
            for j = 0:ev.J
                for k = 1:round(2^(j-1))
                    nott = find((1:ev.T)~=t);  it = (i-1)*ev.T+t;
                    if i == 1
                        Mu = sum(sum( repmat(ev.X(:,t),[1,nr]).*( ev.Ywav(:,indr) - ...
                            repmat(ev.X(:,nott)*betavec(nott),[1,nr]) )./pa.sigma2 ));
                        Sigma = pa.sigma2/( nr + 1/lambda(it) );
                        Mu = Sigma*Mu;
                    end
                    logBF = 0.5*Mu^2/Sigma - 0.5*log(1+nr*lambda(it));
                    Odds = exp(logBF)*pa.pis(r,j+1)/(1-pa.pis(r,j+1));
                    alpha = 1;  if ~isinf(Odds);  alpha = Odds/(1+Odds);  end
                    u = rand(1);  if u <= alpha;  betavec(it) = normrnd(Mu, sqrt(Sigma));  else  betavec(it) = 0;  end
                    t = t+1;
                end
            end
        end
    else  % if no shrinkage, just block update from multivariate normal.
        Sigma = pa.sigma2./( nr + 1./lambda );   Mu = Sigma .* sum(ev.Ywav(:,indr), 2)' / pa.sigma2; 
        betavec =  normrnd(Mu, sqrt(Sigma)); 
    end
    pa.beta(rT) = betavec;
end
end


function [pa] = updateSigma2(pa, ev)
res = zeros(ev.T, ev.N); betapart = 0; p0 = 0;
for r = 1:pa.d
    lambda = pa.lambda(r, ev.indJ+1); rT = (r-1)*ev.pT + (1:ev.pT); betavec = pa.beta(rT);
    p0 = p0 + sum(betavec~=0);
    res(:,pa.labs==r) = ev.Ywav(:,pa.labs==r) - repmat(ev.X*pa.beta((r-1)*ev.pT + (1:ev.pT)),[1,sum(pa.labs==r)]);
    betapart = betapart + sum(betavec.^2./lambda');
end
a_ast = 0.5*(ev.NT + p0) + ev.alpha_sig2;
b_ast = (0.5*sum(sum((res).^2)) + 0.5*betapart + ev.invbeta_sig2)^-1;
pa.sigma2 = 1./gamrnd(a_ast, b_ast);
end


function [x] = updateLambdaP(x, ev)
% update x.lambda(rijk) == x.lambda(rij), i ==1
for r = 1:x.d
    beta = x.beta((r-1)*ev.pT+(1:ev.pT));
    for i = 1:ev.p
        for j = 0:ev.J
            vec = beta(ev.indJ == j);
            x.lambda(r,j+1) = 1./gamrnd(0.5*sum(vec~=0)+ev.alpha_lambda(j+1), ...
                (0.5*sum(vec.^2)/x.sigma2+ev.invbeta_lambda(j+1))^(-1));
            if j > 0 %fix pi at j==0 to be 1
                vec = (vec~=0);
                x.pis(r,j+1) = betarnd(ev.a_pi(j)+sum(vec), ev.b_pi(j)+sum(1-vec));
            end
        end
    end
end
end


function [x] = updateKappa(x, ev)
% update hyperparameters kappa for number of clusters
% u = rand(1);
% x.kappa = 1-(1-u)^(1/(x.d+1));
if ev.useExp
    loglike = ev.kappaconst - x.d*tan(pi/2*ev.kappavec);
else
    loglike = ev.kappaconst + (x.d-1)*log(1-ev.kappavec);
end
MaxLogLike = max(loglike);
P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
U0 = rand(1);
cump = cumsum([0 P(1:(end-1))]);
i0 = sum(U0 > cump);
x.kappa = ev.kappavec(1);
if i0 > 1;   x.kappa = ev.kappavec(i0-1) + ev.kappagap/P(i0)*(U0-cump(i0));   end
end


function [ML1, pa] = getMarginalLoglike(pa, ev, proposeBeta)
L1 = zeros(1,pa.d); L2 = zeros(1,pa.d); rs = kron((1:pa.d)',ones(ev.pT,1)); 
n1 = 0.5*ev.NT+ev.alpha_sig2; ind1s = cell(1,pa.d); nr = histc(pa.labs, 1:pa.d);
for r = 1:pa.d
    indr = find(pa.labs==r);
    betavec = pa.beta(rs==r); ind1 = find(betavec~=0);
    if proposeBeta
        if all(betavec==0)
            gammavec = ones(ev.pT, 1);
            for i = 1:ev.p
                t = 1;
                for j = 0:ev.J
                    pa.pis(r,j+1) = 1;
                    pa.lambda(r,j+1) = 1./gamrnd(ev.alpha_lambda(j+1), 1/ev.invbeta_lambda(j+1));
                    if j > 1;   pa.pis(r,j+1) = betarnd(ev.a_pi(j), ev.b_pi(j));   end
                    if ev.shrinkage;  gammavec(ev.indJ==j) = binornd(1, pa.pis(r,j+1), [1,round(2^(j-1))]);  end
                    t = t+1;
                end
            end
            ind1 = find(gammavec~=0);
        else
            ind1 = find(betavec~=0);
        end
        ind1s{r} = ind1;
    end
    lambda = pa.lambda(r, ev.indJ+1);
    L1(r) = sum(sum(ev.Ywav(:,indr).^2)) - sum(1./(nr(r)+1./lambda(ind1)').*(sum(ev.Ywav(ind1,indr),2)).^2);
    L2(r) = sum(log(nr(r)*lambda(ind1)+1));
end
L1 = sum(L1); L2 = sum(L2);
if proposeBeta
    pa.sigma2 = 1/gamrnd(n1, 1/(ev.invbeta_sig2+L1/2));
    pa.beta = zeros(pa.d*ev.pT,1);
    for r = 1:pa.d
        lambda = pa.lambda(r, ev.indJ+1);
        betavec = zeros(ev.pT,1); ind1 = ind1s{r};
        Mu = sum(ev.X(:,ind1)'*ev.Ywav(:,pa.labs==r), 2);
        Sigma = 1./(nr(r) + 1./lambda(ind1)'); Mu = Sigma.*Mu;
        betavec(ind1) = sqrt(Sigma).*normrnd(0,sqrt(pa.sigma2), [length(ind1),1]) + Mu;
        pa.beta(rs==r) = betavec;
        % plot(Ywav(:,x1.labs==r), 'k.'); hold on; plot(1:pT, X(:,ind1)*betavec(ind1)); hold off
    end
end
% ML1 = -n1*log(1+0.5*L1/ev.invbeta_sig2) - 0.5*L2;
ML1 = -n1*log(L1) - 0.5*L2;
end


function [out] = getFullLoglike(pa, ev)
% get the full log-likelihood (FL)
FL = 0;  SSE = 0; 
nr = histc(pa.labs, 1:pa.d);
rs = kron((1:pa.d)', ones(ev.pT,1));
for r = 1:pa.d
    res = ev.Ywav(:,pa.labs==r) - repmat(ev.X*pa.beta(rs==r), [1, nr(r)]);
    FL = FL - 0.5*nr(r)*ev.T*log(2*pi*pa.sigma2) - 0.5/pa.sigma2*sum(sum(res.^2));
    SSE = SSE + sum(sum(res.^2)); 
end
out.FL = FL;   out.MSE = SSE/ev.NT; 
end


function [labs] = ClusterGen(CenterIndex, ev)
labs = NaN(1, ev.N);
labs(CenterIndex) = 1:length(CenterIndex);
for i = setdiff(1:ev.N, CenterIndex)
    DistMi = ev.D(i, CenterIndex);
    ind = find(abs(DistMi - min(DistMi)) < 1e-10);
    labs(i) = ind(1); % smallest index position
end
end


function [paInitial] = initialCluster(ev, K)
exploreD = 1:10; if K>0; exploreD = K; end
nD = numel(exploreD); paAll = cell(1,nD);  WSS = zeros(1,nD);  %WSS: with-cluster sum of suqared
for k = 1:nD
    pa.d = exploreD(k);
    % pa.center = randsample(ev.N, pa.d);
    [~, ctrs] = kmeans(ev.Ywav',  pa.d);
    for r = 1:pa.d
        pa.center(r) = 1; cdist =  mean((ev.Ywav(:,1)-ctrs(r,:)').^2);
        for i = 1:ev.N
            cdist1 =  mean((ev.Ywav(:,i)-ctrs(r,:)').^2);
            if cdist1 < cdist && (r==1 || (r>1 && all(pa.center(1:(r-1))~=i)))
                pa.center(r) = i;  cdist = cdist1;
            end
        end
    end
    pa.labs = ClusterGen(pa.center, ev);    nr = histc(pa.labs, 1:pa.d);
    for r = 1:pa.d
        Yr = ev.Ywav(:,pa.labs==r); 
        WSS(k) = WSS(k) + sum(sum(  (Yr - repmat(mean(Yr, 2), [1, nr(r)])).^2  )); 
    end
    paAll{k} = pa; 
end
% paInitial = paAll{WSS==min(WSS)};  %choose the smallest WSS
inc = WSS(2:end) - WSS(1:end-1);   i0 = find(inc>0); if isempty(i0); i0 = nD; end;
paInitial = paAll{i0};   %choose the WSS (before first jump-up)
end


function [Ywav, pa] = simulateData(ev)
rng('default');  rng(8);  % set random seed
trueD = 4; 
pa = initialCluster(ev, trueD);
 nr = histc(pa.labs, 1:pa.d);
Ywav = nan(ev.T, ev.N);   pa.sigma2 = 0.6;  pa.beta = nan(ev.T, pa.d);  adj = 5*[-.2, .3, 0, -.4];
for r = 1:pa.d
    indr = find(pa.labs==r);
    pa.beta(:, r) = mean(ev.Ywav(:, indr), 2) + adj(r); 
    for i = 1:nr(r)
        Ywav(:,indr(i)) = pa.beta(:,r) + sqrt(pa.sigma2)*randn(ev.T, 1);  %for intercept only (p=1)
    end
end
plot(ev.Wav'*pa.beta); 
pa.beta = reshape(pa.beta, [numel(pa.beta), 1]); 
pa.Wav = ev.Wav;  
pa.data = ev.Wav'*Ywav; 
end


function [D] = adjance2distance(W)
% compute the distance matrix D from spatial adjancency matrix W
% here the distance of two sites is the number of sites in between
n = size(W,1);  D = zeros(n);
for i = 1:(n-1)
    for j = (i+1):n
        inds = find(W(i,:)~=0); k = 1;
        while all(inds~=j) && k < n
            inds = find(sum(W(inds,:),1)~=0);
            k = k+1;
        end
        D(i,j) = k; D(j,i) = k;
    end
end
end





