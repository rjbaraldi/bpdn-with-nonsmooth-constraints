%% load data
close all;  clear, clc; 
svfile = './figs/'; %file for saving figures

addpath(genpath('./fcns'));
%declare the number of sources and recievers
nsrcx= 102;
nsrcy= 102;
% nsrcx = 51; 
% nsrcy = 51; 
nrecx= 101;
nrecy = 101;
params.nsrcx = nsrcx; 
params.nsrcy = nsrcy;
params.nrecx = nrecx;
params.nrecy = nrecy;
%load data - ask for this, very large file
% load('./freq_40.mat')
%add things you need - you will have to download these
% addpath(genpath('mbin'));
% addpath(genpath('minConf'));
% addpath(genpath('spot'));
% addpath(genpath('pSPOT'));
% addpath(genpath('spgl1SISC'));

addpath(genpath('../alg_tools')); 
addpath(genpath('./fncs')); 

% clear out; 
out = out(1:2:end,1:2:end,1:nsrcy*nsrcx);
%pick which method you use
true = out;

% define restriction

% alpha = [1e-2, 1e-1, 1, 10, 100]; 
alpha = 100*ones(1,5); 
amp = 1; 

%note - change these around for interpolation/denoise combinations
mode.interp =0; 
mode.noise = 1;


perc = linspace(.5, .9, 5);
num_noise = linspace(5, 125, 5);
% num_noise = 10*ones(1,5); 
dataclass = 'l0'; 
splitnum = numel(perc);
psitype = cell(1,splitnum);
psitype(:) = {'l0'}; 



% picmat = zeros(nums*nrecx, nums*nrecy, splitnum+3);
nums = 5;  

poolobj = parpool(splitnum);
eta_fact = 5*ones(1,splitnum);
snr = zeros(1, splitnum); 
snrw = zeros(1, splitnum); 
noise = zeros(1, splitnum); 
snrw(1) = 0;
objvalue = cell(1,splitnum); 
dat_mat_cell = cell(1, splitnum);
true_mat_cell = cell(1, splitnum);  
noise_mat_cell = cell(1, splitnum); 
rec_mat_cell = cell(1, splitnum); 
diff_mat_cell = cell(1, splitnum); 
% file_cell = cell(1, splitnum); 


times = zeros(1,splitnum);

parfor i = 1:splitnum
%% 
params_i = params; 
params_i.alpha = alpha(i);  
[b, params_i.indm, svfile_i] = make_data(params, mode, perc(i), out, svfile);
% add noise

[b, nstot, svfile_i, r] = make_noise(params_i, mode, perc(i),num_noise(i), svfile_i, b, dataclass);
nm =  num2str(perc(i)); 

%% 
disp(svfile_i)
b       = seismic_figs_4d(b, params_i, amp, [1,1], [svfile_i '/' 'obs_4d_' dataclass '_' nm(end)]); 
true_t  = seismic_figs_4d(true, params, amp, [1,1], [svfile_i '/' 'true_4d_' dataclass '_'  nm(end)]);
nstot   = seismic_figs_4d(nstot, params_i, amp, [1,1], [svfile_i '/' 'noise_4d_' dataclass '_'  nm(end)]);
noise(i) = sum(sum(abs(nstot))); 
nstot(nstot~=0)=-1*amp; 
% file_cell{i} = svfile_i; 





%% 
[params_i, sigval] = make_params(params_i, b, dataclass, r, psitype{i}); 


xinit = 1e-6*[vec(params_i.L);vec(params_i.R)];
%numsources in pic

%save data
dat_mat_cell{i} = b(1:nums*nrecx, 1:nums*nrecy);
noise_mat_cell{i} = nstot(1:nums*nrecx, 1:nums*nrecy);
true_mat_cell{i} = true_t(1:nums*nrecx, 1:nums*nrecy);

%save figures
dat_mat_cell{i} = seismic_figs_4d(dat_mat_cell{i}, params_i, amp, [0,1], [svfile_i '/' 'obs_4d_s_' dataclass '_'  nm(end)]); 
true_mat_cell{i} = seismic_figs_4d(true_mat_cell{i}, params_i, amp, [0,1], [svfile_i '/' 'true_4d_s_' dataclass '_'  nm(end)]); 
noise_mat_cell{i} = seismic_figs_4d(noise_mat_cell{i}, params_i, amp, [0,1], [svfile_i '/' 'noise_4d_s_' dataclass '_' nm(end)]); 



    switch(psitype{i})
        case{'spglr'}
            tic
            tau = norm(xinit,1);
            opts              = spgSetParms('project', @TraceNorm_project_hassan, ...
                                'primal_norm', @TraceNorm_primal, 'dual_norm', @TraceNorm_dual, ...
                                'iterations',params_i.stop_crit*params_i.iter_crit,...
                                'proxy', 1, ...
                                'ignorePErr', 1, ...
                                'verbosity',1,...
                                'weights', []);
            opts.funPenalty  = @funLS;
            [xout, r,g, infor]= spgLR(@NLfunforwardCS,vec(b),tau,sigval,xinit,opts,params_i);
            e                = params_i.numr*params_i.nr;
            L                = xout(1:e);
            R                = xout(e+1:end);
            L                = reshape(L,params_i.numr,params_i.nr);
            R                = reshape(R,params_i.numc,params_i.nr);
            objvalue{i} = infor.rNorm; 
            times(i) = toc;
            w = 0; 
        case{'l2', 'l1', 'linf', 'l0'}
            tic
            Ao = opRestriction(nsrcx*nrecx*nsrcy*nrecy, params_i.obs);
            An = opRestriction(nsrcx*nrecx*nsrcy*nrecy, params_i.no_obs);
            [L, R,w, objvalue{i}] = wproj(b, params_i, psitype{i}, sigval,eta_fact(i)); 
            times(i) = toc;
    end
    rec = L*R';
    recw = w; 
    snr(i) = -20*log10(norm(vec(true_t)-rec(:))/norm(vec(true_t)));
    snrw(i) = -20*log10(norm(vec(true_t)-recw(:))/norm(vec(true_t)));
    rec_mat_cell{i} = rec(1:nums*nrecx, 1:nums*nrecy) ; 
    
    
    


    rec_mat_cell{i} = seismic_figs_4d(rec_mat_cell{i}, params, amp, [0,1], [svfile '/' psitype{i} '_4d_z_' dataclass '_'  nm(end)]);
    
    diff_mat_cell{i} = abs(rec_mat_cell{i}-true_mat_cell{i}); 
    diff_mat_cell{i} = seismic_figs_4d(diff_mat_cell{i}, params, amp, [0,1], [svfile '/' 'diff_' psitype{i} '_4d_z_' dataclass '_'  nm(end)]); 
   


end 
disp(snr)
disp(snrw)
disp(noise)
delete(gcp('nocreate')); 
save([svfile '/' 'runinfo' dataclass],'dat_mat_cell', 'true_mat_cell', 'noise_mat_cell', 'rec_mat_cell','rec_mat_cell', 'snr','snrw', 'times', '-v7.3');

