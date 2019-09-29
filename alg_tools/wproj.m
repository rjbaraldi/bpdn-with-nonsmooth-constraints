function [L, R, w, obj_value] = wproj(b, params, projtype, sig, eta_factor)

err = 100; 
eta = params.eta;
% eta_factor = params.eta_fact;
iter_crit = params.iter_crit;
stop_crit = params.stop_crit;
Ao = params.Ao; 
An = params.An; 
L = params.L; 
R = params.R; 
obs = params.obs; 
no_obs = params.no_obs;
converged = params.converged;
k = params.k;
printevery = params.printevery;  
switch(params.method)
    case{'nondist'}
    nsrcx = params.nsrcx; 
    nsrcy = params.nsrcy; 
    nrecx = params.nrecx; 
    nrecy = params.nrecy; 
    case{'2d'}
        numr = params.numr; 
        numc = params.numc; 
    [U, E, V] = svds(reshape(b, [numr, numc]),k);
    L = U*sqrt(E);
    R = V*sqrt(E)';  
end
w = L*R';
count = 0;  
% werr = 0; 
obj_value = norm(L*R', 'fro');
bo = Ao*b(:); 
for i = 1:iter_crit
   while err>converged && count<stop_crit
            Lold      = L;
            Rold      = R;
            %need to speed this up, so replace with chol decomp
            temp = chol(eye(k) + eta*(R'*R)); 
            L = (temp\(temp'\(eta*R'*w')))'; 
            temp = chol((eye(k) + eta*(L'*L))); 
            R = (temp\(temp'\(eta*L'*w)))';
%note: can do svd too, but commented out
%             [U,S,V] = rsvd(w,k);
%             L = U*sqrt(S);
%             R = V*sqrt(S);
            d = vec(L*R');
            wold= w;
            w = w(:); 
            w(no_obs) = An*d;
%             w(no_obs) = d(no_obs); 
            v = Ao*d;
            switch(projtype)
                %note here we switch w to vector for speed (use opres)
                case{'l2'}
                    constraint = norm(v - bo); 
                    if(constraint <= sig)
                        w(obs) = v;
                    else
                        w(obs) = sig*(v-bo)/constraint + bo;
                    end
                    feas = norm(Ao*w - bo) - sig;
                case{'l1'} 
                    if(norm(v-bo,1) <= sig)
                        w(obs) = v;
                    else
                        w(obs) = oneProjector(v-bo,1,sig) + bo; %note link files to OneProjector
                    end
                    feas = norm(Ao*w - bo,1) - sig;
                case{'Hassan'} %testing other functions
                    w(obs) = TraceNorm_project_hassan(v - bo, [], sig, []) + bo; 
                    feas = norm(Ao*w - bo) - sig;
                case{'linf'} 
                    ind = abs(v - bo)<=sig; 
                    w(obs(ind)) = v(ind); 
                    w(obs(~ind)) = sig*sign(v(~ind) - bo(~ind))+ bo(~ind);
                    feas = norm(Ao*w - bo, Inf) - sig;
                case{'l0'}
                    if sig == 0
                        w(obs) = bo;
                        feas = 0; 
                    else
                        z = v-bo;
                        [~, nind] = maxk2(z, sig, 'b'); 
%                         [~, nind] = sort(z, 'descend'); 
%                         nind = setdiff(1:numel(z), max_vals_ind); 
%                         z(nind(sig+1:end)) = 0;
                        z(nind) = 0; 
                        w(obs) = z+bo;
                        feas = abs(numel(find(z~=0))-sig);
                    end 
            end
            werr = norm(Ao*w - v);
            obj_value = [obj_value, norm(d,'fro')]; 
            if(mod(count, printevery)==0)
                fprintf('%s batch: % d iter: %d, w-lr: %7.3f, feas: %7.3e, err: %7.3e\n',projtype, i, count, werr, feas, err);
            end
            %err = [err,norm(w,2)];
            err = normest(Rold(:) - R(:)) + normest(Lold(:)-L(:)) + normest(wold(:) - w(:));
            count = count+1;
            switch(params.method)
                case{'nondist'}
                    w = reshape(w, nrecx*nsrcx, nrecy*nsrcy); 
%                     eta = eta*eta_factor; 
                case{'2d'}
%                     eta = eta*eta_factor;
                    w = reshape(w, numr, numc);
            end
   end
   if strcmp(params.method, '2d')
       eta = eta*eta_factor;
   else
       eta = eta*eta_factor; %reset eta
   end
        count = 0; 
end
        



end