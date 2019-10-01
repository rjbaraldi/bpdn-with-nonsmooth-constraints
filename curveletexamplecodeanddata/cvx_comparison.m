clear;clc;
addpath(genpath('../alg_tools'));
svfile = './figs_cvx/';
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
curv_figs(Res'*b(:),nt,ns, fullfile(svfile, 'noise_missing_data'), 'data');

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
params.x = zeros(size(params.C,2),1); 
params.w1 = zeros(size(params.C,1),1); 
params.w2 = zeros(size(params.R,1),1);
 %initiate pool
psitype = {'spgl1', 'l2', 'l2', 'l0', 'l0', 'cvx'};
phitype={' ', 'l1', 'l0', 'l1', 'l0', 'l1'}; 
splitnum = numel(psitype);


poolobj = parpool(splitnum);


nnnn = norm(b-btrue,2);
nnn1 = norm(b-btrue,1); 
ninf = norm(b-btrue,'inf'); 
sigvals = [nnnn, nnnn, nnnn,floor(numel(noise)*percN), floor(numel(noise)*percN), nnn1]; 

params.iter_crit = 30; 
params.stop_crit = 10; 
params.converged = 1e-10;
params.epsilon = 1e-10;
params.x_switch = 'cg'; 
params.k = size(params.C,1) - numel(find(abs(params.C*true(:))<.1)); 

%compute alpha
params.alpha = norm(params.R*ones(size(params.R,2),1),'fro')^(-2); %10*
params.maxiter = 5;
params.printevery = 10;  
%set eta stuff
params.eta = .1*[1,1]; 
params.eta_factor = [.9,.9];  

datasets = cell(1, splitnum); 
times = zeros(1, splitnum); 
snr = zeros(1, splitnum);
snrw1 = zeros(1, splitnum); 
parfor i = 1:splitnum 
% for i = splitnum
    tic; 
    switch(psitype{i})
        case{'spgl1'}
            opts = spgSetParms('optTol',params.converged,'iterations',params.iter_crit*params.stop_crit, 'verbosity',1);
            [xx,r,g,info]=spgl1(params.R*params.C', b, 0, sigvals(i), [], opts); % Find BP sol'n.1e-3
            xx = opC'*xx; 
        % recovered data
        case{'l2', 'l1', 'linf', 'l0'}
            [xx, w1, w2, out] = prox_grad(b, params,phitype{i}, psitype{i}, sigvals(i));
            w1 = real(opC'*w1);
            Drecw1= reshape(w1,nt,ns);
            snrw1(i) = -20*log10(norm(vec(true)-Drecw1(:))/norm(vec(true)));
        case 'cvx'
            cvx_begin
                variable xx(size(params.C,1))
                minimize(norm(xx, 1))
                subject to
                    norm(params.R*params.C'*xx - b, 1) <= sigvals(i);
                
            cvx_end
            xx = opC'*xx;
            w1 = real(xx);
            Drecw1= reshape(w1,nt,ns);
            snrw1(i) = -20*log10(norm(vec(true)-Drecw1(:))/norm(vec(true))); 
    end
    times(i) = toc;
    Drec= reshape(real(xx),nt,ns);
    snr(i) = -20*log10(norm(vec(true)-Drec(:))/norm(vec(true)));
    datasets{i}.Drec = Drec; 
    datasets{i}.diff = data - Drec;
    curv_figs(datasets{i}.Drec,nt,ns, fullfile(svfile, [phitype{i} '_' psitype{i} '_res']), 'data');
    curv_figs(datasets{i}.diff,nt,ns, fullfile(svfile, [phitype{i} '_' psitype{i} '_diff']), 'data');
end

disp(snr)
disp(snrw1)
disp(times)
save([svfile '/' 'runinfo'], 'snr','snrw1', 'times', '-v7.3');
% delete(gcp('nocreate')); 