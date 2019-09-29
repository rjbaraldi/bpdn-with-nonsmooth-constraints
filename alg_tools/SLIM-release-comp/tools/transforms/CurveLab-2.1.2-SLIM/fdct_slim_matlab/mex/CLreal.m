function [Out] = CLreal(In,flag)

% CLreal - Forward or inverse fdct_wrapping
%
% PRE-REQUISITE
% [Out] = CLrealini(M,N,AC,NBSCALES,NBANGLES_COARSE) or
% [Out] = CLrealini(M,N)
%
% USE
% [Out] = CLreal(In)
%    Forward fdct_wrapping curvelet transform.
% 
% [Out] = CLreal(In,'transp')
%    Inverse fdct_wrapping curvelet transform.
%
% OPTION
%    If global variable VEC_MODE == 1, In and Out
%    always vectors (could be usefull e.g. lsqr).
%
% See CLrealini.

%TO DO:
%   Check size of input according to flag (Gilles)

  global ISREAL

  if(size(ISREAL)==[0,0])
    disp('Please run CLrealini first') %prevents crash
    return
  end

  global HDR
  global CN
  global M
  global N

  global AC
  global NBSCALES
  global NBANGLES_COARSE

  global VEC_MODE
  
  if (nargin > 1)
    if strcmp(flag,'transp')
      % Curvelet coefficient vector
      In  = FDCTVec2Coef_mex(In, HDR);
      In  = fdct_wrapping_r2c(In);
      Out = real(ifdct_wrapping_mex(M,N,NBSCALES,NBANGLES_COARSE,AC,In));
      if (VEC_MODE == 1)
        Out = Out(:);
      end
    else
      disp('Wrong flag for inverse')
      return
    end
  else
    if (VEC_MODE == 1)
      In = reshape(In,M,N);
    end
    % Curvelet transform
    Out = fdct_wrapping_mex(M,N,NBSCALES,NBANGLES_COARSE,AC,double(In));
    Out = fdct_wrapping_c2r(Out);
    Out = FDCTCoef2Vec_mex(Out,CN);
  end
  