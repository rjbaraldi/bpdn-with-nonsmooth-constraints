function [f1,f2] = NLfunforwardCS(x,g,params)
e = params.numr*params.nr;
L = x(1:e);
R = x(e+1:end);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);
if isempty(g)
    f1 = L*R';
    f1(params.ind) = 0;
    f1 = vec(f1);
    f2 = 0;
else
    fp = reshape(g,params.numr,params.numc);
    f1 = [vec(fp*R); vec(fp'*L)];
    f2 = vec(fp);
end
end
