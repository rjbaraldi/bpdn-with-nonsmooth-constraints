clear;clc;close all;
% addpath(genpath('/home/slim/rkumar/Self_Function'));
addpath(genpath('../../Inverpolation/Inverp_trial/mbin'));
addpath(genpath('../../Inverpolation/Inverp_trial/minConf'));
addpath(genpath('../../Inverpolation/Inverp_trial/spot'));
addpath(genpath('../../Inverpolation/Inverp_trial/pSPOT'));
addpath(genpath('../../Inverpolation/Inverp_trial/spgl1SISC'));
addpath(genpath('./oneProjector')); 
%% load data
load ./Green_Function_true1_mat.mat
nt             = 101;%1001
d              = d(1:nt,:,:);
dtrue          = d;
nr             = size(d,2);
nc             = nr;
t = linspace(0,4, nt); 
dataclass = 'l0'; 
% clear d; 
perc           = .4; %used to be .4
data_problem_perc = .1; %used to be .1
pos            = jittersamp_exact(nc,perc);%sample perc%
indperm        = randperm(length(pos)); %random get indices
indperm        = pos(indperm(1:floor(length(pos)*data_problem_perc)));
pos            = setdiff(pos,indperm);
% noiseamp       = .01;
r = length(indperm); 
% d(:,:,indperm) = noiseamp*randn(nt,nr,length(indperm));
mode.noise = 1;
if mode.noise==1
    r = length(indperm); %define length for ease
    for i = 1:nt%only change x-coordinate? 
        switch(dataclass)
            case{'l2'}
                alpha = 1e-2;
                normb = norm(squeeze(d(i,:,:)),2);
                noisemat = alpha*normb*randn(nr, r);
            case{'l1', 'l0'}
                %here we want to add very large noise values to a few
                %points
                alpha = 1e-1; 
                normb = norm(squeeze(d(i,:,:)),1); 
                noisemat = vec(alpha*normb*(randn(nr,r))); %formerly 1e3
                [~, nind] = sort(noisemat, 'descend');
                noisemat(nind(70+1:end))=0; %keep everything but the strongest
                noisemat=reshape(noisemat, nr, r); 
            case{'linf'}
                alpha = 5e-3; 
                normb = norm(squeeze(d(i,:,:)),inf); 
                noisemat = alpha*normb*exp(randn([nr, r])*1i);   
        end
       d(i,:,indperm) = squeeze(d(i,:,indperm)) + noisemat; %add back to values? 
    end
end


%%
R1             = opRestriction(nc,pos);
d              = fft(d);
nf             = floor(nt/2)+1;
d              = d(1:nf,:,:);
d              = permute(reshape(d,nf,nr,nc),[2 3 1]);
%% Restriction operator
R2             = opKron(R1,opDirac(nr));
MH             = opMH(nr,nc);
normstr = {'spgl1', 'l2', 'l1', 'linf', 'l0'};
algtype = cell(12, numel(normstr)); 
for kk = 1:12
    for jj = 1:numel(normstr)
        algtype{kk,jj} = normstr{jj};
    end
end
nk = size(algtype,2);
times = zeros(1,nk);
SNRlr = zeros(1,nk);
SNRw = zeros(1,nk);
alpha = r*numel(find(noisemat==0))*alpha^4*nt; 
% for j = 1:numel(normstr)
for j=3
    tic
    projtype       = distributed(algtype(:,j));
    rank           = distributed(floor(linspace(80,120,nf)));
    freq           = distributed(1:nf);
    dtemp          = distributed(d);
spmd
    codistr      = codistributor1d(3,codistributor1d.unsetPartition,[nr,nc,nf]);
    freqloc      = getLocalPart(freq);
    rankloc      = getLocalPart(rank);
    Dtrueloc     = getLocalPart(dtemp);
    meth         = getLocalPart(projtype); 
    outputloc    = zeros(nr,nc,length(freqloc));
    outputlocw    = zeros(nr,nc,length(freqloc));
    for i   = 1:length(freqloc)
        method_chosen = meth{1}; 
        switch(method_chosen)
            case{'spgl1'}
                iter           = 300; %used to be 300
                [output]= spglr_exp(Dtrueloc(:,:,i),MH,R2,iter,nr,nc,rankloc(i), alpha, dataclass);
                outputw = output; 
            case{'l2', 'l1', 'l0', 'linf'}
                iter           = 300;
                [outlr, outw, obj_value]= wproj_exp(Dtrueloc(:,:,i),MH,R2,iter,nr,nc,rankloc(i), method_chosen, alpha, dataclass);
                output = outlr; 
                outputw = outw; 
        end
        outputloc(:,:,i)  = output;
        outputlocw(:,:,i)  = outputw;
    end
    outputloc= codistributed.build(outputloc,codistr,'noCommunication');
    outputlocw= codistributed.build(outputlocw,codistr,'noCommunication');
end
clear rank freq dtemp codistr freqloc rankloc Dtrueloc
final    = gather(outputloc);
finalw    = gather(outputlocw);
clear outputloc projtype meth method_chosen i iter output outputw outputlocw;
final    = permute(final,[3 1 2]);
finalw   = permute(finalw, [3 1 2]); 
final    = cat(1,final,conj(final(end:-1:2,:,:)));
finalw   = cat(1,finalw, conj(finalw(end:-1:2,:,:))); 
final    = real(ifft(final));
finalw   = real(ifft(finalw)); 
% fs = final(end,:,:); 
% if strcmp(normstr{j}, 'spgl1') || strcmp(normstr{j}, 'l1')
%     save(normstr{j}, 'final','fs', '-v7.3')
% end
% save(normstr{j}, 'fs', '-v7.3'); 
SNRlr(j)    = -20*log10(norm(dtrue(:)-final(:))/norm(dtrue(:)));
SNRw(j)     = -20*log10(norm(dtrue(:)-finalw(:))/norm(dtrue(:)));
times(j)  = toc;
clear final finalw;
delete(gcp('nocreate'))
end
save('final5d', 'SNRlr','SNRw', 'times'); 
