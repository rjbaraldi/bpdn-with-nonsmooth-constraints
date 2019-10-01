clear;clc;
addpath(genpath(pwd));
addpath(genpath('../alg_tools'));
svfile = './figs/'; 
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
% load data
d = load('data.mat');
data = d.data(1:600,:);
[nt,ns] = size(data);
curv_figs(data, nt, ns, fullfile(svfile, 'true_data'), 'data')

% remove missing data
percM = .5;
nomissingindex = jittersamp_exact(ns, percM);
% create restriction operator
Res = opKron(opRestriction(ns,nomissingindex),opDirac(nt));
% add noisy traces
amp = .1;
percN = .01; 
% % create missing and noisy data
tgain = (0:0.004:(nt-1)*0.004);
tgain = repmat(tgain',1,ns);
data  = data.*tgain;
true = data; 
btrue = Res*true(:); 

noise = randn(numel(btrue),1)/amp;
[~, d] = maxk2(noise, floor(numel(noise)*percN), 'b'); 
noise(d) = 0;

b = Res*data(:); 
b = b + noise;

curv_figs(-Res'*noise,nt,ns, fullfile(svfile, 'noise_data'), 'other');
curv_figs(Res'*b(:),nt,ns, fullfile(svfile, 'noise_missing_data'), 'other');

% create curvelet operator
nbs = max(1,ceil(log2(min(nt,ns)) - 3));
nang = 16;
opC = opCurvelet(nt,ns,nbs,nang,1,'ME',0);
p1  = -8/1500;
p2  = -3/1500;
p3  = 3/1500;
p4  = 8/1500;
dipF = opdip(nt,ns,0.004,10,p1,p2,p3,p4);
params.R = Res*dipF;
params.C = opC; 
 %initiate pool
psitype = {'spgl1', 'l2', 'l1', 'linf', 'l0'};
splitnum = numel(psitype); 
phitype = cell(1,splitnum);
phitype(:) = {'l1'}; 


poolobj = parpool(splitnum);

nnnn = norm(b-btrue,2);
nnn1 = norm(b-btrue,1); 
ninf = norm(b-btrue,'inf'); 
sigvals = [nnnn, nnnn, nnn1, ninf, floor(numel(noise)*percN)]; 

params.iter_crit = 30; 
params.stop_crit = 10; 
params.converged = 1e-10;
params.epsilon = 1e-10;
params.x_switch = 'cg'; 

%compute alpha
params.alpha = norm(params.R*ones(size(params.R,2),1),'fro')^(-2); %10*
params.k = size(params.C,1) - numel(find(abs(params.C*true(:))<.1));


params.maxiter = 5;
params.printevery = 10;  
%set eta stuff
params.eta = .1*[1,1]; 
params.eta_factor = [.9,.9]; 



params.x = zeros(size(params.C,2),1); 
params.w1 = zeros(size(params.C,1),1); 
params.w2 = zeros(size(params.R,1),1);


datasets = cell(1, splitnum); 
times = zeros(1, splitnum); 
snr = zeros(1, splitnum);
snrw1 = zeros(1, splitnum); 



parfor i = 1:splitnum 
    tic; 
    switch(psitype{i})
        case{'spgl1'}
            opts = spgSetParms('optTol',params.converged,'iterations',params.iter_crit*params.stop_crit, 'verbosity',1);
            [x,r,g,info]=spgl1(params.R*params.C', b, 0, sigvals(i), [], opts); % Find BP sol'n.1e-3
            x = opC'*x; 
        % recovered data
        case{'l2', 'l1', 'linf', 'l0'}
            [x, w1, w2, out] = prox_grad(b, params, phitype{i}, psitype{i}, sigvals(i));
            w1 = real(opC'*w1);
            Drecw1= reshape(w1,nt,ns);
            snrw1(i) = -20*log10(norm(vec(true)-Drecw1(:))/norm(vec(true)));
    end
    times(i) = toc;
    Drec= reshape(real(x),nt,ns);
    snr(i) = -20*log10(norm(vec(true)-Drec(:))/norm(vec(true)));
    datasets{i}.Drec = Drec; 
    datasets{i}.diff = data - Drec; 
    curv_figs(datasets{i}.Drec,nt,ns, fullfile(svfile, [psitype{i} '_res']), 'data');
    curv_figs(datasets{i}.diff,nt,ns, fullfile(svfile, [psitype{i} '_diff']), 'data');
    
end

disp(snr)
disp(snrw1)
disp(times)
save([svfile '/' 'runinfo' '_l1'], 'snr','snrw1', 'times', '-v7.3');
delete(gcp('nocreate')); 