%% load data
nsrcx             = 102;
nsrcy             = 102;
nrecx             = 101;
nrecy             = 101;
load ../freq_40.mat
true = out;
out = out(1:2:end,1:2:end,:);
%% define restriction
mode.interp = 1;
if mode.interp==1
    perc  = 0.8;
    index = randperm(nsrcx*nsrcy);
    indm  = index(1:floor(nsrcx*nsrcy*perc));
    out(:,:,indm) = 0;
else
    indm = []; 
end
%% add noise
mode.noise = 1;
if mode.noise==1
    perc = 0.01;
    indn = setdiff(1:nsrcx*nsrcy,indm);
    out(:,:,indn(1:floor(length(indn)*perc))) = 1e3*randn(nrecx,nrecy,length(indn(1:floor(length(indn)*perc))))+1i*randn(nrecx,nrecy,length(indn(1:floor(length(indn)*perc))));
end
out               = reshape(permute(reshape(out,nrecx,nrecy,nsrcx,nsrcy),[1 3 2 4]),nrecx*nsrcx,nrecy*nsrcy);
out               = out(1:nsrcx*nrecx,1:nsrcy*nrecy);
b                 = out;
ind               = find(out==0);
b(ind)            = 0;
params.numr       = nsrcx*nrecx;%202*101
params.numc       = nsrcy*nrecy;%202*101
params.nr         = 100;
params.ind        = ind;
params.mode       = 1;
params.ls         = 1;
params.funForward = @NLfunforwardCS;
opts              = spgSetParms('project', @TraceNorm_project_hassan, ...
                                'primal_norm', @TraceNorm_primal, ...
                                'dual_norm', @TraceNorm_dual, ...
                                'proxy', 1, ...
                                'ignorePErr', 1, ...
                                'iterations',250,...
                                'verbosity',2,...
                                'weights', []);
opts.funPenalty  = @funLS;
Linit            = randn(params.numr,params.nr)+1i*randn(params.numr,params.nr);
Rinit            = randn(params.numc,params.nr)+1i*randn(params.numc,params.nr);
xinit            = 1e-6*[vec(Linit);vec(Rinit)];
sigma            = 1e-4*norm(b(:),2);
tau              = norm(xinit,1);
xout             = spgLR(@NLfunforwardCS,b(:),tau,sigma,xinit,opts,params);
e                = params.numr*params.nr;
L                = xout(1:e);
R                = xout(e+1:end);
L                = reshape(L,params.numr,params.nr);
R                = reshape(R,params.numc,params.nr);
rec              = permute(reshape(L*R',nrecx,nsrcx,nrecy,nsrcy),[1 3 2 4]);
-20*log10(norm(true(:)-rec(:))/norm(true(:)))
