function [Out] = CLrealini(M,N,AC,NBSCALES,NBANGLES_COARSE)

% CLrealini - Initialization in order to use CLreal
%
% [Out] = CLrealini(M,N,AC,NBSCALES,NBANGLES_COARSE)
%    Set global variable M,N,AC,NBSCALES,NBANGLES_COARSE and
%    compute the length of the curvelet vector CN (global) 
%    and the curvelet structure HDR (global) according to
%    the given parameters.
% 
% [Out] = CLrealini(M,N)
%    Set global variable M,N. Use CurveLab-2.0 default values
%    and set global variable AC,NBSCALES,NBANGLES_COARSE
%    compute the length of the curvelet vector CN (global) 
%    and the curvelet structure HDR (global) according to
%    the given parameters.
%
% See also fdct_wrapping.

clear global HDR CN M N ISREAL AC NBSCALES NBANGLES_COARSE;

usedefault=0;
  
if(nargchk(5,5,nargin))
  if(nargchk(2,2,nargin))
    disp('Wrong number of input arguments, see correct format below:')
    help CLrealini
    return
  end
    disp('Using default parameter values')
    usedefault=1;
end

if(M<1||N<1 || rem(M,1)~=0 || rem(N,1)~=0 )
  disp('Error: Image dimensions must be positive integers.')
  return
end

if(usedefault==0)
  if( (AC~=0) & (AC~= 1) )
    disp('Error: Allcurvelets must be 0 (wavelets at finest level) or 1 (curvelets).')
    return
  end  
  if(NBSCALES<1 || rem(NBSCALES,1)~=0 )
    disp('Error: Number of scales must be positive integer.')
    return
  end  
  if( NBANGLES_COARSE < 8 || rem(NBANGLES_COARSE,4) ~= 0 )
    disp('Error: Number of angles at 2nd coarsest scale must be 8 or more and be a multiple of 4.')
    return
  end
end%end case where usedefault==0

warning off

global HDR
global CN
global M
global N
global ISREAL
global AC
global NBSCALES
global NBANGLES_COARSE

warning on

ISREAL=1;

if(usedefault==1)
  AC=0;
  NBSCALES=floor(log2(min(M,N)))-3;
  NBANGLES_COARSE=16;
end %end case where usedefault==1

In = randn(M,N);
Out = fdct_wrapping_mex(M,N,NBSCALES,NBANGLES_COARSE,AC,double(In));
Out = fdct_wrapping_c2r(Out);
[CN,HDR] = reshapeFDCTini(Out);
Out = FDCTCoef2Vec_mex(Out,CN);