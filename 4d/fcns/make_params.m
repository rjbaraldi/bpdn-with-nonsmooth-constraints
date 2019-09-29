function [params, sigval] = make_params(params, out, dataclass, r, psitype)

b = out;
ind = find(out==0);
k = 75; %rank(true) = 66
params.k = k; 
alpha = params. alpha; 
params.eta = .01; 
params.count = 0; 
params.printevery = 20;
stop_crit = 50;
iter_crit = 10; 
params.stop_crit = stop_crit; 
params.iter_crit = iter_crit;
params.method = 'nondist';
params.obs = find(b(:)~=0); 
params.no_obs = find(b(:)==0);
params.converged = 1e-10;
params.numr       = params.nsrcx*params.nrecx;
params.numc       = params.nrecy*params.nsrcy;
params.nr         = k; %k=100 set earlier
params.ind        = ind;
params.mode       = 1;
params.ls         = 1;
params.funForward = @NLfunforwardCS;
params.L = 1*randn(params.nsrcx*params.nrecx,k)+1*1i*randn(params.nsrcx*params.nrecx,k);
params.R = 1*randn(params.nrecy*params.nsrcy,k)+1*1i*randn(params.nrecy*params.nsrcy,k);
params.Ao = opRestriction(params.nsrcx*params.nrecx*params.nsrcy*params.nrecy, params.obs);
params.An = opRestriction(params.nsrcx*params.nrecx*params.nsrcy*params.nrecy, params.no_obs);

%% determine which method
switch(dataclass)
    
    case{'l2'}
        coeff = 2*r*numel(find(noisemat~=0)); 
        switch(psitype)
            case{'spglr'}
                sigval = coeff*alpha^2/norm(b,'fro'); 
            case{'l2'}
                sigval = coeff*alpha^2/norm(b,2);
                
            case{'l1'}
                sigval = coeff*1/norm(b,1);
            case{'linf'}
                sigval = coeff*alpha^3/norm(b, 'inf');
            case{'l0'}
                sigval = round(coeff/r); 
    
        end
    case{'l1', 'l0'}
        coeff = 2*r*10; 
        switch(psitype)
            case{'spglr'}
                sigval = coeff*1/norm(b,'fro'); 
            case{'l2'}
                sigval = coeff*alpha/normest(b,2);
                
            case{'l1'}
                sigval = coeff*1;
            case{'linf'}
                sigval = coeff*alpha^3/norm(b, 'inf');
            case{'l0'}
                sigval = round(coeff*.5); 
        end
    case{'linf'}
        coeff = 2*r*numel(find(noisemat~=0));  
        switch(psitype)
            case{'spglr'}
                sigval = coeff*alpha/norm(b,'fro'); 
            case{'l2'}
                sigval = coeff*alpha/norm(b,2);
                
            case{'l1'}
                sigval = coeff*alpha/norm(b,1);
            case{'linf'}
                sigval = coeff*alpha^2/norm(b, 'inf');
            case{'l0'}
                sigval = round(coeff/r); 
        end
        
end




end