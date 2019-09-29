function [output] = spglr_exp( Dtrue,MH,R,iter,nr,nc,rank, alpha, dataclass)
params.numr          = nr;
params.numc          = 2*nc-1;
params.funForward    = @NLfunForward;
params.afunT         = @(x)reshape(x,nr,2*nc-1);
params.mode          = 1;
params.ls            = 1;
params.logical       = 0;
b                    = MH*(R'*(R*vec(Dtrue)));
params.Ind           = find(b==0);
params.afun          = @(x)afun(x,params.Ind,params);
params.nr            = rank;
[U, E, V]            = svds(reshape(b, [params.numr, params.numc]),params.nr);
Linit                = vec(U*sqrt(E));
Rinit                = vec(V*sqrt(E)');
xinit                = 1e-6*[Linit;Rinit];
tau                  = norm(xinit,1);
switch(dataclass)%alpha = r*numel(find(noisemat~=0))*alpha^n
    case{'l2'}
        sigma = alpha/norm(b,'fro')^2;
    case{'l1', 'l0'}
        sigma = alpha/norm(b,'fro')^2;
    case{'linf'}
        sigma = alpha/norm(b,'fro')^2;   
end 
% sigma                = 1e-3*norm(vec(b),'fro');
params.funPenalty    = @funLS;

opts    = spgSetParms('project', @TraceNorm_project_hassan, ...
                      'primal_norm', @TraceNorm_primal, ...
                      'dual_norm', @TraceNorm_dual, ...
                      'proxy', 1, ...
                      'ignorePErr', 1, ...
                      'iterations', iter,...
                      'verbosity', 1); 

xLS                  = spgl1(@NLfunForward,b(:),tau,sigma,xinit,opts,params);
e                    = params.numr*params.nr;
L1                   = xLS(1:e);
R1                   = xLS(e+1:end);
L1                   = reshape(L1,params.numr,params.nr);
R1                   = reshape(R1,params.numc,params.nr);
output               = reshape(MH'*vec(L1*R1'),nr,nc);
end

