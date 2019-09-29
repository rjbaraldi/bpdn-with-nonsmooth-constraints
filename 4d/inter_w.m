%% load data
close all;  clear, clc; 

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
%load
load('./freq_40.mat')
%add things you need
addpath(genpath('../../Inverpolation/Inverp_trial/mbin'));
addpath(genpath('../../Inverpolation/Inverp_trial/minConf'));
addpath(genpath('../../Inverpolation/Inverp_trial/spot'));
addpath(genpath('../../Inverpolation/Inverp_trial/pSPOT'));
addpath(genpath('../../Inverpolation/Inverp_trial/spgl1SISC'));
addpath(genpath('../alg_tools')); 

% clear out; 
out = out(1:2:end,1:2:end,1:nsrcy*nsrcx);
%pick which method you use
true = out;
nstot = zeros(size(out));
% define restriction
amp = 1; 
mode.interp =0; 
if mode.interp==1
    rng(1)
    perc  = 0.8;
    index = randperm(nsrcx*nsrcy);%generate random index for 
    indm  = index(1:floor(nsrcx*nsrcy*perc));%take the first
% loind = indm(1); 
% left_out = squeeze(out(:,:,loind)); 
    out(:,:,indm) = 0;
    svfile = 'inter'; 
else
    indm = [] ;
    svfile = []; 
end
dataclass = 'l0'; 
algtype = {'spglr', 'l2', 'l1', 'linf', 'l0'};
% algtype = {'l0'}; 
splitnum = numel(algtype); 
% add noise
mode.noise = 1;
if mode.noise==1
    rng(2) 
    svfile = [svfile 'deno']; 
%     tot_noise = 0; 
    perc = .8; %used to be 0.01
    idx = setdiff(1:nsrcx*nsrcy,indm);%set of guy's you've observed
    indn = randperm(length(idx)); %random permutation of them
    r = floor(length(indn)*perc); %define length for ease
    %here, we corrupt with really strong points of noise  
    for i = 1:r
        switch(dataclass)
            case{'l2'}
                alpha = 1e-2;
                normb = norm(out(:,:,idx(indn(i))),2);
                noisemat = alpha*normb*(randn(nrecx, nrecy)+1i*randn(nrecx,nrecy));
            case{'l1', 'l0'}
                %here we want to add very large noise values to a few
                %points
                alpha = 1e-1; 
                normb = norm(out(:,:,idx(indn(i))),1); 
                noisemat = vec(alpha*normb*(randn(nrecx, nrecy)+1i*randn(nrecx,nrecy))); %formerly 1e3
                [~, nind] = sort(noisemat, 'descend');
                noisemat(nind(10+1:end))=0; %lose everything but the strongest 10
                noisemat=reshape(noisemat, nrecx, nrecy);
            case{'linf'}
                alpha = 5e-3; 
                normb = norm(out(:,:,idx(indn(i))),inf); 
                noisemat = alpha*normb*exp(randn([nrecx, nrecy])*1i);
        end
        out(:,:,idx(indn(i))) = out(:,:,idx(indn(i))) + noisemat; %add back to values? 
        nstot(:,:,idx(indn(i))) = noisemat;
    end
%     out(:,:,idx(indn(1:floor(length(indn)*perc))))= out(:,:,idx(indn(1:floor(length(indn)*perc))))...
%         + 1e3*(randn(nrecx,nrecy,length(indn(1:floor(length(indn)*perc))))...
%         +1i*randn(nrecx,nrecy,length(indn(1:floor(length(indn)*perc)))));
else 
    perc = .8; %used to be 0.01
    idx = setdiff(1:nsrcx*nsrcy,indm);%set of guy's you've observed
    indn = randperm(length(idx)); %random permutation of them
    r = floor(length(indn)*perc); %define length for ease 
    switch(dataclass)
        case{'l2'}
            alpha = 1e-2;
        case{'l1', 'l0'}
            alpha = 1e-1; 
        case{'linf'}
            alpha = 5e-3; 
    end
end
% load ./data_interp
% save('data_interp', 'out'); 
disp(svfile)
out = reshape(permute(reshape(out,nrecx,nrecy,nsrcx,nsrcy),[1 3 2 4]),nrecx*nsrcx,nrecy*nsrcy);
h1 = figure;set(h1, 'Visible', 'off');imagesc(real(out(1:nsrcx*nrecx,1:nsrcy*nrecy)));colormap gray;caxis([-1 1]*amp);
xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
saveas(h1, ['./' svfile '/' 'obs_4d_' dataclass],'jpg');
true = reshape(permute(reshape(true,nrecx,nrecy,nsrcx,nsrcy),[1 3 2 4]),nrecx*nsrcx,nrecy*nsrcy);
h1 = figure;set(h1, 'Visible', 'off');imagesc(real(true(1:nsrcx*nrecx,1:nsrcy*nrecy)));colormap gray;caxis([-1 1]*amp);
xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
saveas(h1, ['./' svfile '/' 'true_4d_' dataclass],'jpg');

nstot = reshape(permute(reshape(nstot,nrecx,nrecy,nsrcx,nsrcy),[1 3 2 4]),nrecx*nsrcx,nrecy*nsrcy);
h1 = figure;set(h1, 'Visible', 'off');imagesc(real(nstot(1:nsrcx*nrecx,1:nsrcy*nrecy)));colormap gray;caxis([-1 1]*amp);
xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
saveas(h1, ['./' svfile '/' 'noise_4d_' dataclass],'jpg');
nstot(nstot~=0)=-1*amp; 


b = out;
ind = find(out==0);
k = 75; %rank(true) = 66
params.k = k; 
params.eta = .01; 
count = 0; 
params.printevery = 20;
stop_crit = 50;
iter_crit = 3; 
params.stop_crit = stop_crit; 
params.iter_crit = iter_crit;
params.method = 'nondist';
params.obs = find(b(:)~=0); 
params.no_obs = find(b(:)==0);
params.converged = 1e-10;
params.numr       = nsrcx*nrecx;
params.numc       = nrecy*nsrcy;
params.nr         = k; %k=100 set earlier
params.ind        = ind;
params.mode       = 1;
params.ls         = 1;
params.funForward = @NLfunforwardCS;
% picked = vec(true(1:2:end, 1:2:end, :));
%% determine which method
switch(dataclass)
    case{'l2'}
        sigval = 2*r*numel(find(noisemat~=0))*[alpha^2/norm(b,'fro'), alpha^2/norm(b,2), 1/norm(b,1),alpha^3/norm(b, 'inf'), 1/r];
    case{'l1', 'l0'}
%         sigval = 2*r*numel(find(noisemat~=0))*[1/norm(b,'fro'), alpha/norm(b,2), 1, alpha^3/norm(b, 'inf'), .5];
          sigval = 2*r*10*[1/norm(b,'fro'), alpha/normest(b,2), 1, alpha^3/norm(b, 'inf'), .5]; %you picked the strongest 10
    case{'linf'}
        sigval = 2*r*numel(find(noisemat~=0))*[alpha/norm(b,'fro'), alpha/norm(b,2), alpha/norm(b,1),alpha^2/norm(b, 'inf'), 1/r];
        
end
sigval(end) = round(sigval(end)); 
eta_fact = 5*ones(1,splitnum);
snr = zeros(1, splitnum); 
snrw = zeros(1, splitnum); 
snrw(1) = 0;
objvalue = cell(1,splitnum); 
poolobj = parpool(splitnum);

times = zeros(1,splitnum);
params.L = 1*randn(nsrcx*nrecx,k)+1*1i*randn(nsrcx*nrecx,k);
params.R = 1*randn(nrecy*nsrcy,k)+1*1i*randn(nrecy*nsrcy,k);
params.Ao = opRestriction(nsrcx*nrecx*nsrcy*nrecy, params.obs);
params.An = opRestriction(nsrcx*nrecx*nsrcy*nrecy, params.no_obs);
xinit = 1e-6*[vec(params.L);vec(params.R)];
%numsources in pic
nums = 5; 
picmat = zeros(nums*nrecx, nums*nrecy, splitnum+3); 
picmat(:,:,1) = out(1:nums*nrecx, 1:nums*nrecy); 
picmat(:,:,2) = true(1:nums*nrecx, 1:nums*nrecy);
picmat(:, :, 3) = nstot(1:nums*nrecx, 1:nums*nrecy); 
h1 = figure;set(h1, 'Visible', 'off');imagesc(real(squeeze(picmat(:,:,1))));colormap gray;caxis([-1 1]*amp);
xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
saveas(h1, ['./' svfile '/' 'obs_4d_s_' dataclass],'jpg');
h1 = figure;set(h1, 'Visible', 'off');imagesc(real(squeeze(picmat(:,:,2))));colormap gray;caxis([-1 1]*amp);
xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
saveas(h1, ['./' svfile '/' 'true_4d_s_' dataclass],'jpg');

h1 = figure;set(h1, 'Visible', 'off');imagesc(real(squeeze(picmat(:,:,3))));colormap gray;caxis([-1, 0]*amp);
xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
saveas(h1, ['./' svfile '/' 'noise_4d_s_' dataclass],'jpg');



parfor i = 1:splitnum
    switch(algtype{i})
        case{'spglr'}
            tic
            tau = norm(xinit,1);
            opts              = spgSetParms('project', @TraceNorm_project_hassan, ...
                                'primal_norm', @TraceNorm_primal, 'dual_norm', @TraceNorm_dual, ...
                                'iterations',stop_crit*iter_crit,...
                                'proxy', 1, ...
                                'ignorePErr', 1, ...
                                'verbosity',1,...
                                'weights', []);
            opts.funPenalty  = @funLS;
            [xout, r,g, infor]= spgLR(@NLfunforwardCS,vec(b),tau,sigval(i),xinit,opts,params);
            e                = params.numr*params.nr;
            L                = xout(1:e);
            R                = xout(e+1:end);
            L                = reshape(L,params.numr,params.nr);
            R                = reshape(R,params.numc,params.nr);
            objvalue{i} = infor.rNorm; 
            times(i) = toc;
            w = 0; 
        case{'l2', 'l1', 'linf', 'l0'}
            tic
            Ao = opRestriction(nsrcx*nrecx*nsrcy*nrecy, params.obs);
            An = opRestriction(nsrcx*nrecx*nsrcy*nrecy, params.no_obs);
            [L, R,w, objvalue{i}] = wproj(b, params, algtype{i}, sigval(i),eta_fact(i)); 
            times(i) = toc;
    end
    rec = L*R';
    recw = w; 
    snr(i) = -20*log10(norm(vec(true)-rec(:))/norm(vec(true)));
    snrw(i) = -20*log10(norm(vec(true)-recw(:))/norm(vec(true)));
    picmat(:,:,i+3) = rec(1:nums*nrecx, 1:nums*nrecy) 
end 
disp(snr)
disp(snrw)
delete(gcp('nocreate')); 
save(['./' svfile '/' 'runinfo' dataclass],'picmat', 'snr','snrw', 'times', '-v7.3');
for i = 1:splitnum
    h1 = figure;set(h1, 'Visible', 'off');imagesc(real(squeeze(picmat(:,:,i+3))));colormap gray;caxis([-1 1]*amp);
    xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
    ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
    set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
    saveas(h1, ['./' svfile '/' algtype{i} '_4d_z_' dataclass],'jpg');
    h1 = figure;set(h1, 'Visible', 'off');imagesc(real(squeeze(abs(picmat(:,:,2)-picmat(:,:,i+3)))));colormap gray;caxis([-1 1]*amp);
    xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
    ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
    set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
    saveas(h1, ['./' svfile '/' 'diff_' algtype{i} '_4d_z_' dataclass],'jpg');
    
end
%%
