clear; close all; 
addpath(genpath(pwd));
addpath(genpath('../alg_tools')); 
compound = 1;
rng(2); 
m = compound*120; n = compound*512; k = compound*20; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
R  = randn(m,n); [Q,~] = qr(R',0);  R = Q';
amp = 1;
perc = .1; 
noise = amp*randn(m,1); 
[~, d] = maxk2(noise, perc*m, 'b'); 
noise(d) = .001*randn(sum(d),1); 
b  = R*x0 + noise;

psitype = {'true', 'spgl1', 'l2', 'l1', 'linf', 'l0'};
phitype = {'l1'}; 
splitnum = numel(psitype); 
sigvals = [0, amp*norm(noise,2), amp*norm(noise,2), amp*norm(noise,1), amp*norm(noise, inf), perc*m]; 
params.converged = 1e-8; 
params.R = R; 
params.k = k; 
params.epsilon = 1e-6; 
params.x_switch = 'grad_des'; 
switch params.x_switch
    case{'grad_des'}
        params.eta = .2*ones(2,1); 
        params.eta_factor = .5*ones(2,1); 
        params.alpha = 10*norm(R, 'fro')^(-1);
        svfile = './figs/';
    case{'cg'}
        params.eta = [1000 1000]*.9^80;
        params.eta_factor = [.9, .9]; 
        params.alpha = norm(R, 'fro')^(-2);
        svfile = './figs/';
end
 
% params.iter_crit = 250 ;
% params.stop_crit = 500;
params.iter_crit = 100; 
params.stop_crit = 100;
params.printevery = 500;%floor(params.stop_crit/2); 
params.x = randn([size(params.R,2),1]); 
params.w1 = randn([size(params.R,2),1]); 
params.w2 = randn([size(params.R,1),1]);

params.C = eye(size(params.R, 2)); 
datasets = cell(1, splitnum);  
outparams = cell(1, splitnum); 
dw = cell(1, splitnum); 
times = zeros(1, splitnum); 
snr = zeros(1, splitnum);
snrw1 = zeros(1, splitnum); 
snrw2 = zeros(1, splitnum); 
for i = 1:splitnum
    switch(psitype{i})
        case 'true'
            x = x0; 
            datasets{i}=x;
            dw{i} = params.R*x; 
            
        case 'spgl1'
            opts = spgSetParms('optTol',params.converged,'iterations',params.iter_crit*params.stop_crit, 'verbosity',1);
            [x,r,g,info]=spgl1(params.R*params.C, b, 0, sigvals(i), [], opts);
            outparams{i} = info;
            dw{i} = params.R*x; 
        % recovered data
        case{'l2', 'l1', 'linf', 'l0'}
            [x, w1, w2, out] = prox_grad(b, params,phitype{1}, psitype{i}, sigvals(i));
            dw{i} = w2; 
            snrw1(i) = -20*log10(norm(x0-w1)/norm(vec(x0)));
            snrw2(i) = -20*log10(norm(params.R*(params.C*x0)-w2)/norm(vec(params.R*(params.C*x0))));
            outparams{i} = out;
    end
    if numel(x)~=numel(x0)
        x = params.x; 
    end
    snr(i) = -20*log10(norm(x0-x)/norm(vec(x0)));
    datasets{i} =x; 
    
    
end
% colors = {'k','b', 'r', 'g', 'c', 'm'};
colors = {'k-', 'k-*', 'k-.', 'k--', 'k->', 'k:', 'k-x'};
names = {'True', 'SPGL1', 'l_2', 'l_1', 'l_infty', 'l_0'}; 
l1deno_fig_create(b, colors{1},'obs',svfile, 'w2');
for i = 1:splitnum
    l1deno_fig_create(datasets{i}, colors{i},names{i},svfile, 'data');
    l1deno_fig_create(dw{i}, colors{i}, names{i}, svfile, 'w2'); 
    if i>2
        l1deno_fig_create(outparams{i}.werr', colors{i},names{i}, svfile, 'werr');
        l1deno_fig_create(outparams{i}.witer_err', colors{i},names{i}, svfile, 'witer');
        l1deno_fig_create(outparams{i}.obj_value', colors{i},names{i}, svfile, 'obj');
        l1deno_fig_create(outparams{i}.constants', colors{i},names{i}, svfile, 'const');
        l1deno_fig_create(outparams{i}.feas', colors{i},names{i}, svfile, 'feas');
    end
end
disp(snr)
disp(snrw1)
disp(snrw2)

