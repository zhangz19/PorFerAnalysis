function [out] = stm(Y, X, Z, W, ev)
%y = xbeta + Zu + epsilon, u: NT by 1 spatio-temporal AR1(phi) kron CAR(gamma)

[pa, ev] = initPara(Y, X, Z, W, ev);  %pa stores parameters;  E stores updated environmental variables (won't further update) 

% MCMC running
out.matPara = zeros(ev.nsample, ev.p+4); 
out.Es = [];   Us = zeros(ev.nsample, ev.NT);
if ev.computeDIC == 1;  out.Es = zeros(ev.nsample, 2);  end  % 2 components of logLik for DIC4

V = ev.M - W.*pa.gamma; %LV = chol(V,'lower'); LV = LV';
V1 = V\eye(ev.N); LV = chol(V1,'lower'); LV = LV\eye(ev.N);

tic
for iter = 1:ev.niter
    if ev.verbose == 1
        %fprintf('%6d', b) %fprintf('%3d%3.2f', [b, x.phi])
        fprintf('%3.2f %3.2f %3.2f %3.2f %3.2f\n', [pa.phi, pa.gamma, pa.sigma2, pa.tau2, mean(pa.u)])
        % if ~mod(b,20);  fprintf('\n');  end
    end
    [pa] = updateBeta(pa, Y, X, Z, ev);
    [pa,res] = updateSigma2(pa, Y, X, Z, ev);
    if ev.nonrandom ~= 1
        [pa,delta,res] = updateTau2(pa, LV, res, ev);
        if ev.nontemporal ~= 1
            [pa] = updatePhi(pa, LV, ev);
        end
        if ev.nonspatial ~= 1
            [pa,V,LV] = updateGamma(pa, W, ev);
        end
        [pa] = updateU(pa, V, delta, res, ev);
    end
    if(iter > ev.burnin)
        out.matPara((iter-ev.burnin),:) = [pa.beta', pa.sigma2, pa.tau2, pa.phi, pa.gamma];
        Us((iter-ev.burnin),:) = pa.u';
        if ev.computeDIC == 1
            out.Es(iter-ev.burnin,:) = getDIC(pa, Y, X, Z, W, ev, LV);
        end
    end
end

CPUtime = toc; out.CPUtime = CPUtime/60;  % in minute
fprintf('\n%d iterations are done with elapsed time %.2f minutes.\n', ev.niter, out.CPUtime)
if ev.saveAsMat; save(ev.saveAsMatFileName, 'out'); end
end

function [pa, ev] = initPara(Y, X, Z, W, ev)
ev.M = diag(max(1, sum(W, 1)));
ev.N = size(W, 1);  ev.NT = numel(unique(Z));  ev.T = ev.NT/ev.N;  
ev.J = numel(Y)/ev.NT;
ev.J0 = floor(log2(ev.T));  %for wavelet decomposition, round sample size to a power of 2
ev.nsample = ev.niter - ev.burnin;
switch ev.betaprior 
    case 1 % flat prior
        ev.X = [kron(eye(ev.J), ones(ev.NT,1)), kron(eye(ev.J), kron((1:ev.T)',ones(ev.N,1)))];
        ev.p = size(ev.X,2);
        invXtX = (ev.X'*ev.X)\eye(ev.p);  ev.facMu = invXtX*X';  ev.L0 = chol(invXtX, 'lower');
    case 2 % construct wavelet transformation matrix
        ev.p = 1; %no functional covariates in this version
        tmp = eye(ev.T); Wav = tmp;  for i1 = 1:T;  Wav(:,i1) = wavedec(tmp(:,i1),J0,'db1');  end
        ev.X = Wav';
end
ev.pT = ev.p*ev.T;

% set priors
meansigma2 = 0.01; varsigma2 = 10^2;
ev.alphasig = 2+meansigma2^2/varsigma2;  ev.invbetasig = meansigma2*(ev.alphasig-1);

meantau2 = 0.01; vartau2 = 10^2;
ev.alphatau = 2+meantau2^2/vartau2;  ev.invbetatau = meantau2*(ev.alphatau-1);

invM = inv(ev.M);  eigs = eig(sqrt(invM)*W*sqrt(invM));
ev.lgamma = max(1/min(eigs),-1);  ev.ugamma = 1/max(eigs);
ev.gap_gamma = 1e-2;  ev.gammas = (ev.lgamma+ev.gap_gamma):ev.gap_gamma:(ev.ugamma-ev.gap_gamma); 
ev.len = numel(ev.gammas);
ev.lgamma0 = zeros(1, ev.len);
for i = 1:ev.len; ev.lgamma0(i) = 0.5*ev.T*sum(log(eig(ev.M-ev.gammas(i)*W))); end

ev.gap_phi = 1e-2;  ev.phis = (-1+ev.gap_phi):ev.gap_phi:(1-ev.gap_phi);
ev.lphi0 = - 0.5*(ev.T-1)*ev.N*log(1-ev.phis.^2);

% set initial values
pa.gamma = 0.7337; %0
pa.phi = 0;
pa.tau2 = 0.0583; %1e-10;
pa.u = normrnd(0,sqrt(pa.tau2),[ev.NT,1]); %zeros(NT,1);
switch ev.betaprior
    case 1
        pa.beta = ev.facMu*Y;   pa.sigma2 = mean((Y - X*pa.beta).^2);
    case 2
        pa.beta = zeros(ev.pT,1);  pa.sigma2 = 0;
        for i = 1:ev.N
            ywave = wavedec(Y(:,i),ev.J0,'db1');
            pa.beta = pa.beta + ywave;
            pa.sigma2 = pa.sigma2 + sum((Y(:,i)-ywave).^2);
        end
        pa.beta = pa.beta./ev.N;  pa.sigma2 = pa.sigma2./ev.NT;
end
ev.eta2 = 100+zeros(ev.p, ev.T);
% pia = 6; pib = 1; pis = betarnd(pia, pib, [p J+1]); 
ev.pis = 0.5 + zeros(ev.p, ev.J0+1);

if ev.nonspatial == 1;  pa.gamma = 0; ev.M = eye(ev.N);  end
if ev.nontemporal == 1;  pa.phi = 0;  end
if ev.nonrandom == 1;  pa.u = zeros(ev.NT,1);  pa.tau2 = 0;  end
end

function [pa] = updateBeta(pa, Y, X, Z, ev)
switch ev.betaprior
    case 1
        Mu = ev.facMu*(Y-pa.u(Z));
        tmp = normrnd(0,1,[ev.p,1]);
        pa.beta = sqrt(pa.sigma2).*(ev.L0*tmp) + Mu;
    case 2
        for i = 1:ev.p
            t = 1;
            for j = 0:ev.J
                for k = 1:round(2^(j-1))
                    ind0 = false(1,ev.T); ind0(t) = true;
                    inds = false(1, ev.pT); inds(t + (i-1)*ev.T) = true;
                    
                    if i == 1
                        Mu = sum(sum( repmat(X(:,ind0),[1,ev.N]).*...
                            ( Y - pa.u - repmat(X(:,~ind0)*pa.beta(~ind0),[1,ev.N]) )./pa.sigma2 ));
                        Sigma = 1/( ev.N/pa.sigma2 + 1/ev.eta2(i,t) );
                        Mu = Sigma*Mu;
                    end
                    
                    logBF = 0.5*Mu^2/Sigma - 0.5*log(ev.eta2(i,t)) + 0.5*log(Sigma);
                    Odds = exp(logBF)*ev.pis(i,j+1)/(1-ev.pis(i,j+1));
                    alpha = 1;
                    if ~isinf(Odds);  alpha = Odds/(1+Odds);   end
                    u = rand(1);
                    if u <= alpha
                        pa.beta(inds) = normrnd(Mu, sqrt(Sigma));
                    else
                        pa.beta(inds) = 0;
                    end
                    t = t+1;
                end
            end
        end
end
end

function [pa, res] = updateSigma2(pa, Y, X, Z, ev) 
res = Y - X*pa.beta;
azeros = 0.5*numel(res) + ev.alphasig;
bzeros = (0.5*sum((res-pa.u(Z)).^2) + ev.invbetasig)^-1;
pa.sigma2 = 1./gamrnd(azeros, bzeros);
end

function [pa, delta, res] = updateTau2(pa, LV, res, ev)
U = reshape(pa.u, [ev.N, ev.T]); U = U';
uMu = fastKronMulti(LV, pa.phi, U);
bzeros = (0.5*uMu + ev.invbetatau)^-1;
azeros = 0.5*ev.NT + ev.alphatau;
pa.tau2 = 1./gamrnd(azeros, bzeros);
delta = pa.tau2/pa.sigma2;
res = sum(reshape(res, [ev.NT, ev.J]), 2).*delta;
delta = delta*ev.J;
res = reshape(res, [ev.N, ev.T]);  res = res';
end

function [pa] = updatePhi(pa, LV, ev)
U = reshape(pa.u, [ev.N, ev.T]);  U = U';
loglike = (-0.5/pa.tau2).*fastKronMulti(LV, ev.phis, U) + ev.lphi0;
MaxLogLike = max(loglike);
P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
U0 = rand(1);
cump = cumsum([0 P(1:(end-1))]);
i0 = sum(U0 > cump);
pa.phi = ev.phis(1);
if i0 > 1;  pa.phi = ev.phis(i0-1) + ev.gap_phi/P(i0)*(U0-cump(i0));  end
end

function [pa, V, LV] = updateGamma(pa, W, ev)
U = reshape(pa.u, [ev.N, ev.T]);  U = U';
fac1 = 0;
for t = 1:(ev.T-1)
    if t == 1
        a = 1;
    else
        a = 1+pa.phi^2;
    end
    fac1 = fac1 + a.*(U(t,:)*W*U(t,:)') - 2*pa.phi.*(U(t,:)*W*U(t+1,:)');
end
fac1 = fac1 + U(ev.T,:)*W*U(ev.T,:)';
loglike = ev.lgamma0 + (0.5/pa.tau2)*fac1/(1-pa.phi^2)*ev.gammas;
MaxLogLike = max(loglike);
P = exp(loglike-MaxLogLike)/sum(exp(loglike-MaxLogLike));
U0 = rand(1);
cump = cumsum([0 P(1:(end-1))]);
i0 = sum(U0 > cump);
pa.gamma = ev.gammas(1);
if i0 > 1;  pa.gamma = ev.gammas(i0-1) + ev.gap_gamma/P(i0)*(U0-cump(i0));  end
V = ev.M - W.*pa.gamma; %LV = chol(V,'lower'); LV = LV';
V1 = V\eye(ev.N); LV = chol(V1,'lower'); LV = LV\eye(ev.N);
end

function [pa] = updateU(pa, V, delta, res, ev)
ur = zeros(ev.N, ev.T); tmpmu = zeros(ev.N, ev.T); mu = zeros(ev.N, ev.T);
h = (1-pa.phi^2)^-1;
L0 = cell(1, ev.T); TL = cell(1, ev.T);
ur0 = normrnd(0, sqrt(pa.tau2), [ev.N, ev.T]);
% forward step
for t = 1:ev.T
    if t == 1
        Star = h*V + delta.*eye(ev.N);
        TL{t} = chol(Star, 'lower'); L0{t} = TL{t}\eye(ev.N);
        ur(:,t) = (ur0(:,t)'*L0{t})';
        TL{t} = -pa.phi*h*V*L0{t}'; TinvStar = TL{t}*L0{t};
        tmpmu(:,t) = L0{t}*res(t,:)';
    else
        tmpmu(:,t) = res(t,:)' - TL{t-1}*tmpmu(:,t-1); % TL: t-1
        if t>1 && t<ev.T
            Star = (1+pa.phi^2)*h*V + delta.*eye(ev.N) + pa.phi*h*TinvStar*V; % TinvStar: t-1
        else
            Star = h*V + delta.*eye(ev.N) + pa.phi*h*TinvStar*V;
        end
        TL{t} = chol(Star, 'lower'); L0{t} = TL{t}\eye(ev.N);
        tmpmu(:,t) = L0{t}*tmpmu(:,t);
        ur(:,t) = (ur0(:,t)'*L0{t})';
        TL{t} = -pa.phi*h*V*L0{t}'; TinvStar = TL{t}*L0{t}; % update TL, TinvStar to t
    end
end

% backward step
Ls = cell(1,ev.T-1);
for t0 = 1:ev.T
    t = ev.T+1-t0;  % from T to 1, the column
    if t == ev.T
        mu(:,t) = L0{t}'*tmpmu(:,t);  % L0: T
    else
        for t2 = (t+1):ev.T
            if t2 == t+1 % last value: t+2, L0 info from t+1
                Ls{t} = L0{t+1};
            end
            Ls{t2-1} = -Ls{t2-1}*TL{t}*L0{t};
            ur(:,t) = ur(:,t) + (ur0(:,t2)'*Ls{t2-1})';
        end
        mu(:,t) = L0{t}'*(tmpmu(:,t) - TL{t}'*mu(:,t+1)); % info at t
    end
end
ur = ur + mu;
pa.u = reshape(ur, [ev.NT,1]);
end

function [uMu] = fastKronMulti(LD, phi, ur)
% Evaluate uMu = ur'{kron(invD, invB)}ur
% note ur is T by nr matrix
d = length(phi); % phi is 1 by d vector
[T,nr] = size(ur); uMu = zeros(1,d);
h = 1./sqrt(1-phi.^2);
hphi = phi.*h; 

if length(size(LD)) < 3
    for s1 = 1:nr
        delta = zeros(T,d);
        for s2 = 1:s1
            if d > 1
                urs = ur(:,s2);
                urnew = repmat(urs, [1 d]);
                urnew(2:T,:) = kron(h, urs(2:T)) - kron(hphi, urs(1:(T-1)));
                delta = delta + urnew.*LD(s1,s2);
            else
                urs = ur(:,s2); urnew = urs;
                urnew(2:T) = (urs(2:T) - phi*urs(1:(T-1)))*h;
                delta = delta + urnew.*LD(s1,s2);
            end
        end
        uMu = uMu + sum((delta).^2, 1); %delta'*delta
    end
else
    for s1 = 1:nr
        delta = zeros(T,size(LD,3));
        for s2 = 1:s1
            urs = ur(:,s2); urnew = urs;
            urnew(2:T) = (urs(2:T) - phi*urs(1:(T-1)))*h;
            delta = delta + kron(reshape(LD(s1,s2,:),[1 size(LD,3)]), urnew);
        end
        uMu = uMu + sum((delta).^2, 1);
    end
end
end

function [Es] = getDIC(pa0, Y, X, Z, W, ev, LV)
nreplicate = 50; pa = pa0; Es = zeros(1,2);  
U = reshape(pa.u, [ev.N, ev.T]);  U = U';
res = Y - X*pa.beta - pa.u(Z);
Q = numel(res); 
logLikPart1 = -0.5*sum((res).^2)/pa.sigma2 - 0.5*Q*log(2*pi*pa.sigma2);
if ev.nonrandom == 1
    logLikPart2 = 0;
else
    lgamma = 0.5*ev.T*sum(log(eig(ev.M-pa.gamma*W)));
    lphi = - 0.5*(ev.T-1)*ev.N*log(1-pa.phi^2);
    logLikPart2 = (-0.5/pa.tau2).*fastKronMulti(LV, pa.phi, U) + lphi + lgamma - 0.5*ev.NT*log(2*pi*pa.tau2);
end
Es(1) = logLikPart1+logLikPart2;
betabar = zeros(ev.p,1); sigma2bar = 0; tau2bar = 0; phibar = 0; gammabar = 0;
for i = 1:nreplicate
    [pa] = updateBeta(pa, Y, X, Z, ev);  betabar = betabar + pa.beta;
    [pa,res] = updateSigma2(pa, Y, X, Z, ev);   sigma2bar = sigma2bar + pa.sigma2;
    if ev.nonrandom ~= 1 
        [pa, ~, ~] = updateTau2(pa, LV, res, ev);   tau2bar = tau2bar + pa.tau2;
        if ev.nontemporal ~= 1
            [pa] = updatePhi(pa,  LV, ev);  phibar = phibar + pa.phi;
        end
        if ev.nonspatial ~= 1
            [pa, ~, LV] = updateGamma(pa, W, ev);  gammabar = gammabar + pa.gamma;
        end
    end
end
betabar = betabar/nreplicate; sigma2bar = sigma2bar/nreplicate; tau2bar = tau2bar/nreplicate;
phibar = phibar/nreplicate; gammabar = gammabar/nreplicate;
res = Y - X*betabar - pa.u(Z);
logLikPart1 = (-0.5/sigma2bar)*sum((res).^2) - 0.5*Q*log(2*pi*sigma2bar);
if ev.nonrandom == 1
    logLikPart2 = 0;
else
    lgamma = 0.5*ev.T*sum(log(eig(ev.M-gammabar*W)));
    lphi = - 0.5*(ev.T-1)*ev.N*log(1-phibar^2);
    V = ev.M - W.*gammabar; %LV = chol(V,'lower'); LV = LV';
    V1 = V\eye(ev.N);  LV1 = chol(V1,'lower');  LV1 = LV1\eye(ev.N);
    logLikPart2 = (-0.5/tau2bar)*fastKronMulti(LV1, phibar, U) + lphi + lgamma - 0.5*ev.NT*log(2*pi*tau2bar);
end
Es(2) = logLikPart1 + logLikPart2;
end
