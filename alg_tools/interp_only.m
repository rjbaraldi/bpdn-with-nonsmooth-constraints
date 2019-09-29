clear;clc;close all;
% addpath(genpath('/home/slim/rkumar/Self_Function'));
addpath(genpath('../../Inverpolation/Inverp_trial/spgl1SISC'));
addpath(genpath('../../Inverpolation/Inverp_trial/minConf'));
addpath(genpath('../../Inverpolation/Inverp_trial/mbin'));
addpath(genpath('../../Inverpolation/Inverp_trial/spot'));
addpath(genpath('../../Inverpolation/Inverp_trial/pSPOT'));
addpath(genpath('./oneProjector')); 
%% load data
load ./Green_Function_true1_mat.mat
d              = d(1:1001,:,:);
dtrue          = d;
nt             = 1001;
nr             = size(d,2);
nc             = nr;
t              = 0:0.004:4;
clear d; 
load ./5d_data_interponly.mat

%%
R1             = opRestriction(nc,pos);
d              = fft(d);
nf             = floor(nt/2)+1;
d              = d(1:nf,:,:);
d              = permute(reshape(d,nf,nr,nc),[2 3 1]);
%% Restriction operator
R2             = opKron(R1,opDirac(nr));
MH             = opMH(nr,nc);
normstr = {'spgl1', 'l2', 'l1', 'l0', 'infty'};

algtype = cell(12, numel(normstr)); 
for kk = 1:12
    for jj = 1:numel(normstr)
        algtype{kk,jj} = normstr{jj};
    end
end
nk = size(algtype,2);
times = zeros(1,nk);
SNR = zeros(1,nk);

for j = 3:nk
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
    for i   = 1:length(freqloc)
        method_chosen = meth{1}; 
        switch(method_chosen)
            case{'spgl1'}
                iter           = 300; %used to be 300
                [output]= spglr_exp(Dtrueloc(:,:,i),MH,R2,iter,nr,nc,rankloc(i));
            case{'l2', 'l1', 'l0', 'infty'}
                iter           = 300;
                [outlr, outw, obj_value]= wproj_exp(Dtrueloc(:,:,i),MH,R2,iter,nr,nc,rankloc(i), method_chosen,j);
                if strcmp(method_chosen, 'l2')
                    output = outw; 
                else
                    output = outlr; 
                end
        end
        outputloc(:,:,i)  = output;
    end
    outputloc= codistributed.build(outputloc,codistr,'noCommunication');
end
clear rank freq dtemp codistr freqloc rankloc Dtrueloc
final    = gather(outputloc);
clear outputloc projtype meth method_chosen i iter output;
final    = permute(final,[3 1 2]);
final    = cat(1,final,conj(final(end:-1:2,:,:)));
final    = real(ifft(final));
fs = final(end,:,:); 
if strcmp(normstr{j}, 'spgl1') || strcmp(normstr{j}, 'l1')
    save([normstr{j} '_interponly'], 'final','fs', '-v7.3')
end
% save(normstr{j}, 'fs', '-v7.3'); 
SNR(j)    = -20*log10(norm(dtrue(:)-final(:))/norm(dtrue(:)));
times(j)  = toc;
clear final;
delete(gcp('nocreate'))
end
save('final5d', 'SNR', 'times'); 
