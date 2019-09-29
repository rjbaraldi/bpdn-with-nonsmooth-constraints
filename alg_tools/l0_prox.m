function [y,trace] = l0_prox(x, alpha, b, meth)
n = size(x,1); 
%unconstrained
if strcmp(meth, 'unconstrained')
% ind = (.5*w1p.^2)>(1/alpha); %used to be eta 1
% ind = abs(w1p)>sqrt(2*alpha);
    if nargout > 1
        y = (abs(x) > sqrt(2*alpha)); % indices nonzero entries for y. allocate y.
        numnonzeros = nnz(y);
        y = x.*y;
        trace.y = y;
        trace.c = alpha*numnonzeros + 0.5*norm(y-x,2)^2;
        trace.m = numel(x) - numnonzeros;
    else
        y = x .* (abs(x) > sqrt(2*alpha));
    end
    
%nonnegative
elseif strcmp(meth, 'nonneg')
    if nargout > 1    
        y = (x > sqrt(2*alpha)); % indices nonzero entries for y. allocate y.
        numnonzeros = nnz(y);
        y = x.*y;
        trace.y = y;
        trace.c = alpha*numnonzeros + 0.5*norm(y-x,2)^2;
        trace.m = numel(x) - numnonzeros;
    else
        y = x.*(x > sqrt(2*alpha));
    end
    

elseif strcmp(meth, 'constrained')
    if (b > 0)
    %
    elseif (b < 0)
        error('prox_l0_simplex: b should be nonnegative');
    else % b == 0
        y = zeros(size(x));
        if nargout > 1
            trace.y = y;
            trace.c = 0.5*norm(x,2)^2;
            trace.m = numel(x);
        end
        return
    end
    if isscalar(x)
        y = b;
        if nargout > 1
            trace.y = y;
            trace.c = alpha + 0.5*(b-x)^2;
            trace.m = 0;
        end
        return
    else
        x = x(:);
        nx = length(x);
    end

    % sorting
    [xs,idx] = sort(x,'ascend'); % xs = x(idx)
    iidx(idx) = 1:1:nx; % x = xs(iidx)
    % number of zeros
    mv = (0:1:(nx-1))'; % 0 <= m <= n-1
    % multiplier
    lv = (b - cumsum(xs,'reverse'))./(nx-mv);
    % feasibility
    feas = (xs+lv > 0.0);
    mv = mv(feas); % m in M
    lv = lv(feas);
    % cumulative sum of squares
    sv = cumsum([0;xs(1:end-1)].^2);
    sv = sv(feas);
    % cost
    cv = alpha*(nx-mv) + 0.5*sv + 0.5*(nx-mv).*lv.^2;
    % minimization
    [c,i] = min(cv); % least sparse solution
    % [c,i] = min(flipud(cv)); % most sparse solution
    % solution
    m = mv(i);
    l = lv(i);
    ys = [zeros(m,1); xs(m+1:end) + l];
    y = ys(iidx); % re-ordering
    if nargout > 1
        trace.m = m;
        trace.c = c;
        trace.l = l;
        trace.y = y;
        trace.M = mv;
        trace.C = cv;
        trace.L = lv;
        nm = length(mv);
        trace.Y = zeros(nx,nm);
        for i=1:nm
            ytmp = [zeros(mv(i),1); xs(mv(i)+1:end) + lv(i)];
            trace.Y(:,i) = ytmp(iidx);
        end
    end

    
    
end






end