clear all;clc;close all;

% add SLIM toolbox
SLIM_COMP.CurveLab
SLIM_APPS.tools.utilities.SPOT_SLIM
SLIM_APPS.tools.utilities.segyMAT
SLIM_APPS.tools.Miscellaneous

% load any data 
D             = ReadSegyFast('/localData/shared/BDR04-1120F1-219_cable5_gun1.sgy');

% time samples
nt            = 3451; 

% sources
ns            = 882;

% channels
nc            = 478; 

% reshape data
D             = reshape(D,nt,nc,ns);

% extract a shot gather 
sourceindex   = 100;
D1            = D(:,:,sourceindex);

% Define number of scales and angles
nbscales      = max(1,ceil(log2(min(nt,nc)) - 3));
nbangles      = 16;

% Forward curvelet transform
Cs            = mefcv2(D1,nt,nc,nbscales,nbangles);

% display coefficients
img           = fdct_me_dispcoef(Cs,nt,nc,nbscales);

% since the display function creates lots of zeros to fit all the scales
% and angles, its better to make its value very small for visualization
% purpose
ind           = find(abs(img)==0.5);
img(ind)      = 0.01;

% display the coefficients 
h = figure;imagesc(abs(img));caxis([0 0.07]);pbaspect([1 1 1])
set(gca,'plotboxaspectratio',[1 1.5 1]);
xlabel('# of Curvelet coefficients','FontSize',16, 'FontName','helvetica','FontWeight','demi')
ylabel('# of Curvelet coefficients','FontSize',16, 'FontName','helvetica','FontWeight','demi')
set(gca,'fontsize',16,'FontName','helvetica','FontWeight','demi')
