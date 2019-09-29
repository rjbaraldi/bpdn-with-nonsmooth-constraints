function [x, w1, w2, obj_value] = bpdn(b, params, projtype, sig)
%eta_factor should have two values, as should eta
err = 100; 
eta1 = params.eta(1);
eta2 = params.eta(2);
eta1_factor = params.eta_factor(1);
eta2_factor = params.eta_factor(2); 
iter_crit = params.iter_crit;
stop_crit = params.stop_crit;
R = params.R;
C = params.C; 
x = params.x(:); 
w1 = params.w1; 
w2 = params.w2; 
converged = params.converged;
printevery = params.printevery; 
epsilon = params.epsilon; 
if isfield(params, 'maxiter')
    maxiter = params.maxiter;
else
    maxiter = numel(b); 
end

%initialize the rest of the algorithm
RtR = R'*R;
CtC = C'*C; 
EyeTestC = sqrt(sum(sum(CtC - eye(size(CtC,1)) < epsilon)));
Rtb = R'*b; 
count = 0;   
obj_value = norm(x, 1);
for i = 1:iter_crit
   if  EyeTestC==size(R,2) && size(R,2)<200
%        disp('Computing Eigenvals')
        if isa(RtR, 'opFoG') && isa(R, 'opRestriction')
            [V, D] = eigs(eta1/eta2*RtR, size(RtR,2)); 
        else
            [V, D] = eig(eta1/eta2*RtR); 
        end
   else
%        disp('CG')
        Abig = eta1/eta2*RtR + CtC;
    
   end
   while err>converged && count<stop_crit
            %x update
            xold = x; 
            bbig = eta1/eta2*R'*w2 + C'*w1 + eta1/eta2*Rtb;
            if EyeTestC==size(R,2) && size(R,2)<200
                x = (V*diag(1./(diag(D)+1))*V')*bbig; 
            else
                x = fastcg(Abig, xold,  bbig, epsilon, maxiter); 
            end
            %w1 update
            w1old = w1; 
            w1 = sign(C*x).*max(0, abs(C*x)-eta1); 
            
            %w2 update
            w2old = w2; 
            w2 = R*x - b; %+ or - b? 
            
            switch(projtype)
                %note here we switch w to vector for speed (use opres)
                case{'l2'}
                    t = norm(w2,2);
                    w2 = w2/max(1, t/sig);
                    feasw = norm(w2) - sig;
                    feas = norm(R*x - b, 2)-sig;
                case{'l1'} 
                    w2 = double(w2); 
                    w2 = oneProjector(w2,1,sig); 
                    feasw = norm(w2,1) - sig;
                    feas = norm(R*x - b, 1)-sig;
                case{'Hassan'} %testing other functions
                    w2 = TraceNorm_project_hassan(w2, [], sig, []); 
                    feasw = norm(w2) - sig;
                    feas = norm(R*x - b, 2)-sig;
                case{'linf'} 
                    ind = abs(w2)>sig; 
                    w2(ind) = sig*sign(w2(ind));
                    feasw = norm(w2, Inf) - sig;
                    feas = norm(R*x - b, Inf)-sig;
                case{'l0'}
                    if sig == 0
                        w2 = b;
                        feasw = 0; 
                    else
                        [~, nind] = maxk2(abs(w2), sig, 'b');
                        w2(nind) = 0;
                        feasw = abs(numel(find(w2~=0))-sig);
                    end 
                    feas = abs(numel(find((R*x-b)~=0))-sig); 
            end
            werr = norm(w1 - x);
            w2err = norm(R*x - b - w2); 
            obj_value = [obj_value, norm(x,1)];
            
            if(mod(count, printevery)==0)
               fprintf('%s batch: % d iter: %d, w1-x: %7.3f, w2 - (Ax - b): %7.3f, w-feas: %7.3e, x-feas: %7.3e, err: %7.3e\n', projtype, i, count, werr,w2err, feasw,feas, err);
            end
%             if feas<0
%                 break
%             end
            %err = [err,norm(w,2)];
            err = norm(w2old - w2) + norm(w1old-w1) + norm(xold - x);
            count = count+1;
   end
   eta1 = eta1*eta1_factor;
   eta2 = eta2*eta2_factor;
   count = 0; 
   err = 100; 
end
 



end


function x = fastcg(A, x, b, epsilon, maxiter)
x = double(x); 
b = double(b); 
r = double(b - A*x);
p = double(r);  
rsold = r'*r; 
for i = 1:maxiter
    Ap = A*p; 
    alpha = rsold/(p'*Ap); 
    x = x+alpha*p; 
    r = r-alpha*Ap; 
    rsnew = r'*r; 
    if sqrt(rsnew)<epsilon
        break; 
    end
    p = r+(rsnew/rsold)*p; 
    rsold = rsnew; 
    
    
end



end