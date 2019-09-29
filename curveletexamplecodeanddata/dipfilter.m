function datafilt = dipfilter(data,dt,dx,p1,p2,p3,p4,n1,n2)

%DIPFILTER - apply a dip-filter in the FK domain to the input data
%
%   Syntax datfilt = dipfilter(data,dt,dx,p1,p2,p3,p4)
%      data      : input data in x-t domain
%      dt        : time step
%      dx        : offset step
%      p1        : all slopes p<p1 are set to zero
%      p2        : slopes p1<p<p2 are tapered from 0 to 1
%      p3        : slopes p2<p<p3 are passed
%      p4        : slopes p3<p<p4 are tapered from 1 to 0
%
%   The four p-values define a pass-band filter like below:
%
%                      ----------
%                     /          \
%                    /            \
%                   /              \
%               -------------------------
%                  p1 p2        p3 p4
%
%   The edges of the band-pass filter are cos^2 functions
%   The output array datfilt is the dip-filtered data in x,t domain
% 
%   e.g; one set of values of p1,p2,p3 and p4 is 
%                             (-1/1500,-0.1/1500,0.1/1500,1/1500)
%
%   Notes : internally, the data is zero-padded: nkx=2*nx;
%           in this way spatial wrap-around is avoided
%
%   Author: Eric Verschuur
%   Date  : April 2011
%

%-------------------------------------------------------------
% check number of input arguments
%-------------------------------------------------------------

if ( nargin < 7 ) 
   error('Syntax datafilt = dipfilter(data,dt,dx,p1,p2,p3,p4)')
end

%-------------------------------------------------------------
% make sure p1<p2<p3<p4
%-------------------------------------------------------------

if (p2 < p1 || p3 < p2 || p4 < p3)
   error('parameters p1, p2, p3, p4 must be monotonic')
end

%-------------------------------------------------------------
% get size of input data
%-------------------------------------------------------------
data = reshape(data,n1,n2);
nt = size(data,1);
nx = size(data,2);

%-------------------------------------------------------------
% get number of frequencies and kx-values
% pad double number of traces to avoid wrap around
%-------------------------------------------------------------

nf = nt;
df = 1/(dt*nt);
nkx = 2*nx;
dkx = 1/(dx*nkx);

%-------------------------------------------------------------
% define frequency axes as a column vector
% make f=0 into f=0.1*dp to avoid division by 0
% create a 2D grid version of this axis
%-------------------------------------------------------------

% original by Eric
% f=(0:nf-1)'*df;
% f(1)=0.1*df;
% f(nf/2+2:nf)=-f(nf/2:-1:2);

% added by AJ
fquist=0.5/dt;
nfq=nquist(nt);
f(1:nfq) = 0:df:fquist;
if mod(nfq,2) == 0
    f(nfq+1:nf)=-f(nfq:-1:2);
else
    f(nfq+1:nf)=-f(nfq-1:-1:2);
end    
f = f';
f(1)=0.1*df;
% make fgrid
fgrid=f*ones(1,nkx);

%-------------------------------------------------------------
% define kx axes as a row vector
% flip pos and neg kx values, because fft2 has flipped them
% create a 2D grid version of this axis
%-------------------------------------------------------------

kx=(0:-1:-nkx+1)*dkx;
kx(nkx/2+2:nkx)=-kx(nkx/2:-1:2);
kxgrid=ones(nf,1)*kx;

pgrid=kxgrid./fgrid;

%-------------------------------------------------------------
% create dipfilter based on pgrid for pos. freqs only
%-------------------------------------------------------------

% pgridhalf=pgrid(1:nf/2+1,:);
pgridhalf=pgrid(1:nfq,:);

pfilthalf=pgridhalf;
aa= pgridhalf<=p1;
pfilthalf(aa)=0.0;
aa=find(pgridhalf>p1 & pgridhalf<=p2);
pfilthalf(aa)=sin((pgridhalf(aa)-p1)*pi/(2*(p2-p1))).^2;
aa= pgridhalf>p2 & pgridhalf<=p3;
pfilthalf(aa)=1.0;
aa=find(pgridhalf>p3 & pgridhalf<=p4);
pfilthalf(aa)=cos((pgridhalf(aa)-p3)*pi/(2*(p4-p3))).^2;
aa= pgridhalf>=p4;
pfilthalf(aa)=0.0;

pfilt=zeros(nf,nkx);
% pfilt(1:nf/2+1,:)=pfilthalf;
pfilt(1:nfq,:)=pfilthalf;

%-------------------------------------------------------------
% copy to full grid in pfilt
%-------------------------------------------------------------

for ifreq=2:nf/2
   pfilt(nf+2-ifreq:nf+2-ifreq,2:nkx)=pfilt(ifreq:ifreq,nkx:-1:2);
end

%-------------------------------------------------------------
% apply filter to data in FK domain
%-------------------------------------------------------------

tmp=fft2(ifft2(data,nt,nkx).*pfilt); % Don't know why Eric first applied ifft2 then fft2, may be it doesn't matter. 

% tmp=ifft2(fft2(data,nt,nkx).*pfilt);

datafilt=real(tmp(1:nt,1:nx));
datafilt=datafilt(:);
end

