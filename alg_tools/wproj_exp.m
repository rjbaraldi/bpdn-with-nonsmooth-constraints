function [outlr, outw, obj_value] = wproj_exp( Dtrue,MH,R,itermax,nr,nc,rank, projtype, alpha, dataclass)
params.numr = nr;
params.numc = 2*nc-1;
b = MH*(R'*(R*vec(Dtrue)));
scale = 1;
[U, E, V]= svds(reshape(b, [params.numr, params.numc]),nr);
params.L = vec(U*sqrt(E));
params.R= vec(V*sqrt(E)');
b = b/scale; 
params.k = rank;
params.eta = .001; 
params.printevery = 20; 
params.iter_crit = 10;
params.stop_crit = round(itermax/params.iter_crit); 
  
params.obs = find(vec(b)~=0); 
params.no_obs = find(vec(b)==0);
% params.eta_fact = 1.03; 
eta_fact = 3;%4
% params.eta_fact = 2.5; 
switch(dataclass)
    case{'l2'}
        sig = alpha*[ 1/norm(b,2), 1/norm(b,1),1/norm(b, 'inf'), 1/nr];
    case{'l1', 'l0'}
        sig = 1*[ 1/norm(b,2), norm(b,1)/alpha, 1/norm(b, 'inf')^3, 1/numel(find(b==0))/2];
    case{'linf'}
        sig = alpha*[ 1/norm(b,2), 1/norm(b,1)^2,1/norm(b, 'inf'), 1/nr];   
end 
switch(projtype)
    case('l2')
        sig = sig(1); 
    case('l1')
        sig = sig(2); 
    case('linf')
        sig = sig(3); 
    case('l0')
        sig = round(sig(4)); 
end 
params.method = '2d';
params.Ao = opRestriction(params.numr*params.numc, params.obs);
params.An = opRestriction(params.numr*params.numc, params.no_obs);
params.converged = 1e-5;
[L, R, w, obj_value] = wproj(b, params, projtype, sig, eta_fact);
L = scale*reshape(L,params.numr,params.k);
R = scale*reshape(R,params.numc,params.k);
w = scale*w; 
outlr = reshape(MH'*vec(L*R'),nr,nc);
outw  = reshape(MH'*vec(w),nr,nc); 
end

