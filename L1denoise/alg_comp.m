clear; close all; 
addpath(genpath(pwd));
addpath(genpath('../alg_tools')); 
compound = 1;
rng(2); 
  m = compound*120; n = compound*512; k = compound*20; % m rows, n cols, k nonzeros.
  p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
  R  = randn(m,n); %[Q,~] = qr(R',0);  R = Q';
  amp = 1;
  perc = .1; 
  noise = amp*randn(m,1); 
  [~, d] = maxk2(noise, perc*m, 'b'); 
  noise(d) = .001*randn(sum(d),1); 
  b  = R*x0 + noise;

psitype = {'l1', 'l1', 'l1', 'l1', 'l1'};
splitnum = numel(psitype); 

% switchparms = {'alg1', 'inter1', 'inter5', 'inter20', 'max'};
% alpha = norm(R, 'fro').^(linspace(-1,-2, 5));
alpha = norm(R, 'fro').^(-2)*ones(1, splitnum); 
maxiter = [0, 1, 2, 20, numel(b)]; 
svfile = '/Users/bobby/Dropbox (uwamath)/TextOnly/IEEEversion_final/L1deno/comparison/';

sigvals = amp*norm(noise,1)*ones(size(psitype)); 
params.converged = 1e-16; 
params.R = R; 
params.eta = 1e-4*ones(2,1); 
params.eta_factor = ones(2,1);
params.k = k; 
params.printevery = 500;%floor(params.stop_crit/2); 
params.iter_crit = 1 ;
params.stop_crit = 2000;


params.x = randn([size(params.R,2),1]); 
params.w1 = randn([size(params.R,2),1]); 
params.w2 = randn([size(params.R,1),1]);
params.C = eye(size(params.R, 2)); 
params.epsilon = 1e-16; 
params.eta_factor = ones(2,1); 
datasets = cell(1, splitnum); 
outparams = cell(1, splitnum); 
dw = cell(1, splitnum); 
times = zeros(1, splitnum); 
snr = zeros(1, splitnum);
snrw1 = zeros(1, splitnum); 
snrw2 = zeros(1, splitnum); 
for i = 1:splitnum
    if i==1
        params.x_switch = 'grad_des'; 
    else
        params.x_switch = 'cg'; 
    end
    disp('--------------------------------')
    params.alpha = alpha(i); 
    params.maxiter = maxiter(i); 
    [x, w1, w2, out] = prox_grad(b, params,'l1', psitype{i}, sigvals(i));
    dw{i} = w2; 
    snrw1(i) = -20*log10(norm(x0-w1)/norm(vec(x0)));
    snrw2(i) = -20*log10(norm(params.R*(params.C*x0)-w2)/norm(vec(params.R*(params.C*x0))));
    outparams{i} = out;
    if numel(x)~=numel(x0)
        x = params.x; 
    end
    snr(i) = -20*log10(norm(x0-x)/norm(vec(x0)));
    datasets{i} =x; 
    
    
end
colors = {'k', 'k', 'k', 'k', 'k'};
dot = {'-', '--', '-.', ':', '-*'}; 
names = {'Alg (1)', '1', '5', '20', 'Alg (3)'}; 
h1 = figure;
hold on
for i = 1:splitnum
    loglog(outparams{i}.obj_total,[colors{i},dot{i}], 'LineWidth', 1.5)
end
    pbaspect([5,1,1])
    xlabel('k^{th} Iteration', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
    ylabel('Objective Value: Eq (4)', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
    set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold','XScale', 'log', 'YScale', 'log');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'Position', [500 500 1000 500]);
    set(gcf, 'PaperSize', [1000 200]); 
    set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    set(gcf,'color','w');
    legend(names, 'Location', 'northeast')
    img = frame2im(getframe(h1)); 
    imwrite(img, [svfile, 'algcomp_results.jpg']); 
hold off
disp(snr)
disp(snrw1)
disp(snrw2)
