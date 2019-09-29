function [Out] = CL(In,flag)

% CL - Forward or inverse fdct_wrapping (COMPLEX)
%
% PRE-REQUISITE
% [Out] = CLini(M,N,AC,NBSCALES,NBANGLES_COARSE) or
% [Out] = CLini(M,N)
%
% USE
% [Out] = CL(In)
%    Forward fdct_wrapping curvelet transform.
% 
% [Out] = CL(In,'transp')
%    Inverse fdct_wrapping curvelet transform.
%
% OPTION
%    If global variable VEC_MODE == 1, In and Out
%    always vectors (could be usefull e.g. lsqr).
%
% See CLini.

%TO DO:
%   Check size of input according to flag (Gilles)

  global ISREAL

  if(size(ISREAL)==[0,0])
    disp('Please run CLini first') %prevents crash
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
      if(ISREAL==1)
	In  = fdct_wrapping_r2c(In);
      end
      Out = ifdct_wrapping_mex(M,N,NBSCALES,NBANGLES_COARSE,AC,In);
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
    Out = FDCTCoef2VecCpx_mex(Out,CN);
  end
