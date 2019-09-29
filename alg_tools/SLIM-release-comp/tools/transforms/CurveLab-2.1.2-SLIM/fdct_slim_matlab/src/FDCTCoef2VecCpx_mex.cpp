/*
 * =============================================================
 * FDCTCoef2VecCpx_mex.cpp
 *
 * AUTHOR
 * Gilles Hennenfent
 * University of British Columbia
 * ghennenfent@eos.ubc.ca
 * 
 * DEVELOPMENT
 * May 5th, 2005: First version
 * Aug 30th, 2005: Updated to accommodate complex coefficients (Colin Russell)
 * =============================================================  
 */

#include "mex.h"
#include "math.h"

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int n = mxGetN(prhs[0]);
  double*jp = mxGetPr(prhs[1]);// CN = # of coeffs
  
  plhs[0] = mxCreateDoubleMatrix((int)*jp,1, mxCOMPLEX);
  double*dr = mxGetPr(plhs[0]);
  double*di = mxGetPi(plhs[0]);
  
  int cnt = 0;
  // Loop over scales
  for (int i = 0; i < n; i++) {
    const mxArray*y = mxGetCell(prhs[0],i);
    int nn = mxGetN(y);
    for (int j = 0; j < nn; j++){
      const mxArray*yy = mxGetCell(y,j);
      /* Check that input is complex. */
      if ( !mxIsComplex(yy) )
	mexErrMsgTxt("Input must be complex.\n");
      double*xr = mxGetPr(yy);
      double*xi = mxGetPi(yy);
      int mm = mxGetM(yy)*mxGetN(yy);
      for (int k = 0; k < mm; k++){
	*(dr+cnt+k) = *(xr+k);
	*(di+cnt+k) = *(xi+k);
      }
      cnt += mm;
    }
  }
  return;
}
