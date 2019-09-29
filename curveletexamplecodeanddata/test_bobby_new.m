clear;clc;
addpath(genpath(pwd));
addpath(genpath('../alg_tools'));
svfile = '/home/rbaraldi/Dropbox (uwamath)/TextOnly/IEEEversion/curv_temp/';
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
% load /d/home/rajivk/Downloads/shots_base.mat 
% data                                    = shots_all(1:400,20:140,80);
load data.mat
data                                    = data(1:600,:);
[nt,ns]                                 = size(data);
figure(1);imagesc(data,[-0.5,0.5]);colormap gray
% remove missing data
percM                                   = 0.5;
nomissingindex                          = jittersamp_exact(ns,percM);
% create restriction operator
Res                                     = opKron(opRestriction(ns,nomissingindex),opDirac(nt));
% create missing and noisy data
tgain                                   = (0:0.004:(nt-1)*0.004);
tgain                                   = repmat(tgain',1,ns);
data                                    = data.*tgain;
b                                       = Res*data(:);
figure(3);imagesc(reshape(Res'*b(:),nt,ns),[-0.5,0.5]);colormap gray

% create curvelet operator
nbs                                     = max(1,ceil(log2(min(nt,ns)) - 3));
nang                                    = 16;
opC= opCurvelet(nt,ns,nbs,nang,1,'ME',0);
p1                                      = -8/1500;
p2                                      = -3/1500;
p3                                      = 3/1500;
p4                                      = 8/1500;
dipF                                    = opdip(nt,ns,0.004,10,p1,p2,p3,p4);
A                                       = Res*dipF*opC';
iter                                    = 200;
opts                                    = spgSetParms('iterations',iter);
[x,r,g,info]                            = spgl1(A, b, 0, 0, [], opts); % Find BP sol'n.
% recovered data
Drec                                    = reshape(real(opC'*x),nt,ns);
SNR                                     = -20*log10(norm(data(:)-Drec(:))/norm(data(:))) %NOTE: this compares it to corrupted data
figure(4);imagesc(Drec,[-0.5,0.5]);colormap gray
figure(5);imagesc(data-Drec,[-0.5,0.5]);colormap gray
% figure(6);imagesc(reshape(Res'*(Res*(data(:)-Drec(:))),nt,ns),[-0.1,0.1]);colormap gray
