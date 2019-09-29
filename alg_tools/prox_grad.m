function [x, w1, w2, out] = prox_grad(b, params,phitype, psitype, sig)
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
x = params.x; 
w1 = params.w1; 
w2 = params.w2;
k = params.k; 
converged = params.converged;
printevery = params.printevery; 
x_switch = params.x_switch; 
w2err = norm(R*x - b - w2); 
werr = norm(w1 - C*x);

if isfield(params, 'maxiter')
    maxiter = params.maxiter;
else
    maxiter = numel(b); 
end

%initialize the rest of the algorithm
RtR = R'*R;
CtC = C'*C;
Rtb = R'*b; 
epsilon = params.epsilon; 
count = 0;   
out.obj_value = [norm(w1,1); norm(x, 1)];
out.witer_err =zeros(3,1); 
out.werr = [werr; w2err];
out.constants = [params.alpha; eta1; eta2];
out.feas = zeros(2,1);
out.obj_total = norm(w1, 1)+1/(2*eta1)*werr^2 + 1/(2*eta2)*w2err^2; 
for i = 1:iter_crit
        nCtC = (1/eta1)*CtC; 
        nRtR = (1/eta2)*RtR;
        if isfield(params,'alpha') && strcmp(params.x_switch, 'grad_des')
            alpha = params.alpha*eta1;
%         elseif isfield(params, 'alpha') && strcmp(phitype, 'l0')
%             alpha = params.alpha; 
        else 
            alpha = min(eta1, eta2); 
        end
   while err>converged && count<stop_crit
            %x update
            xold = x; 
            bbig = 1/eta2*R'*w2 + 1/eta1*C'*w1 + 1/eta2*Rtb;
            switch x_switch
                case{'grad_des'}
                    x = xold - alpha*(nCtC*x - 1/eta1*C'*w1+nRtR*x - 1/eta2*R'*w2 - 1/eta2*Rtb);
                case{'cg'}
                    x = fastcg(nCtC + nRtR, xold,  bbig, epsilon, maxiter);
%                     w = warning('off', 'all'); 
%                     warning(w); 
%                     [x,~] = pcg(nCtC + nRtR, bbig, epsilon, maxiter, [], [], xold); 
            end
            %w1 update gb
            w1old = w1; 
            w1p = w1 - alpha/eta1*(w1 - C*x); 
            
            switch(phitype)
                case{'l1'}
                    w1 = sign(w1p).*max(0, abs(w1p)-alpha); 
                case{'l2'}
                    w1 = (1 - (1/alpha)/max(norm(w1p), (1/alpha)))*w1p; 
                case{'l0'}
                    %hard thresholding? 
                    w1 = l0_prox(w1p, alpha, k, 'unconstrained');
%                     w1 = l0_prox(w1p, alpha, k, 'constrained'); 
%                     if mod(count, 4)==0
%                         figure;hold on
%                         plot(w1);
%                         plot(w1old)
%                         hold off
%                         figure;hold on
%                         plot(x);
%                         plot(xold);
%                         hold off
%                         figure;
%                         plot(abs(w1p));
%                     end
                case{'linf'}
                    w1p = double(w1p); 
                    w1 = w1p -  oneProjector(w1p,1,1/alpha); 
                    
            
            end
            
            %w2 update
            w2old = w2;
            w2 = w2old - alpha/eta2*(w2old-(R*x - b)); 
              %+ or - b? 
            
            switch(psitype)
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
            werr = norm(w1 - C*x);
            w2err = norm(R*x - b - w2); 
            out.obj_value = [out.obj_value, [norm(w1,1); norm(x,1)]]; %note obj value is wrt w1
            out.obj_total = [out.obj_total, norm(w1,1)+1/(2*eta1)*werr^2 + 1/(2*eta2)*w2err^2]; 
            out.witer_err = [out.witer_err, [norm(w1-w1old); norm(w2-w2old); norm(xold-x)]];
            out.constants = [out.constants, [alpha; eta1;eta2]]; 
            out.werr = [out.werr, [werr; w2err]]; 
            out.feas = [out.feas, [feas; feasw]]; 
            err = norm(w2old - w2) + norm(w1old-w1) + norm(xold - x);
            if(mod(count, printevery)==0)
               fprintf('%s batch: % d iter: %d, w1-Cx: %7.3f, w2 - (Ax - b): %7.3f, w-feas: %7.3e, x-feas: %7.3e, err: %7.3e, eta: %7.3e \n', psitype, i, count, werr,w2err, feasw,feas, err, eta1);
            end
%             if feas<0
%                 break
%             end
            %err = [err,norm(w,2)];
            
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
count = 0; 
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
    count=count+1;    
    
end

end