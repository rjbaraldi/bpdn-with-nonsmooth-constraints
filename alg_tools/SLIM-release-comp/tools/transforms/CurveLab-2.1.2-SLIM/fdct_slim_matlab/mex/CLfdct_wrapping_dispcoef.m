function img = CLfdct_wrapping_dispcoef(C)

% CLfdct_wrapping_dispcoef - Returns image of curvelet coefficients
%                            given by CLreal
%
% PRE-REQUISITE
% [Out] = CLrealini(M,N,AC,NBSCALES,NBANGLES_COARSE) or
% [Out] = CLrealini(M,N)
%
% [Out] = CLreal(In)
%
% USE
% [Out] = CLfdct_wrapping_dispcoef(In)
%    Returns image of curvelet coefficients.
%
% See also CLrealini, CLreal, fdct_wrapping_dispcoef.

  global ISREAL

  if(size(ISREAL)==[0,0])
    disp('Please run CLrealini first') %prevents crash
    return
  end

  global M
  global N

  global AC
  global NBSCALES
  global HDR

  C  = FDCTVec2Coef_mex(C, HDR);
  C  = fdct_wrapping_r2c(C);

  m = M;
  n = N;
  nbscales = NBSCALES;
  
  img = C{1}{1};  img = img/max(max(abs(img))); %normalize
  if(~AC)
    nbscales = nbscales-1;
  end
  for sc=2:nbscales
    nd = length(C{sc})/4;
    wcnt = 0;
    
    ONE = [];
    for w=1:nd
      ONE = [ONE, C{sc}{wcnt+w}];
    end
    wcnt = wcnt+nd;
    
    TWO = [];
    for w=1:nd
      TWO = [TWO; C{sc}{wcnt+w}];
    end
    wcnt = wcnt+nd;
    
    THREE = [];
    for w=1:nd
      THREE = [C{sc}{wcnt+w}, THREE];
    end
    wcnt = wcnt+nd;
    
    FOUR = [];
    for w=1:nd
      FOUR = [C{sc}{wcnt+w}; FOUR];
    end
    wcnt = wcnt+nd;
    
    [p,q] = size(img);
    [a,b] = size(ONE);
    [g,h] = size(TWO);
    m = 2*a+g;    n = 2*h+b; %size of new image
    scale = max(max( max(max(abs(ONE))),max(max(abs(TWO))) ), max(max(max(abs(THREE))), max(max(abs(FOUR))) )); %scaling factor
    
    new = 0.5 * ones(m,n); %background value
    new(a+1:a+g,1:h) = FOUR /scale;
    new(a+g+1:2*a+g,h+1:h+b) = THREE /scale;
    new(a+1:a+g,h+b+1:2*h+b) = TWO /scale;
    new(1:a,h+1:h+b) = ONE /scale; %normalize
    
    dx = floor((g-p)/2);    dy = floor((b-q)/2);
    
    new(a+1+dx:a+p+dx,h+1+dy:h+q+dy) = img;
    
    img = new;
  end
