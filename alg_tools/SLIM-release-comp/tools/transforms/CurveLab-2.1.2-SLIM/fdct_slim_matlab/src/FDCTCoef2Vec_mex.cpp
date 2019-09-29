/*
 * =============================================================
 * FDCTCoef2Vec_mex.cpp
 *
 * AUTHOR
 * Gilles Hennenfent
 * University of British Columbia
 * ghennenfent@eos.ubc.ca
 * 
 * DEVELOPMENT
 * May 5th, 2005: First version
 * =============================================================  
 */

#include "mex.h"
#include "math.h"

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
 int n = mxGetN(prhs[0]);
 double*jp = mxGetPr(prhs[1]);

 plhs[0] = mxCreateDoubleMatrix((int)*jp,1, mxREAL);
 double*dr = mxGetPr(plhs[0]);

 int cnt = 0;
 // Loop over scales
 for (int i = 0; i < n; i++) {
        const mxArray*y = mxGetCell(prhs[0],i);
	int nn = mxGetN(y);
	  for (int j = 0; j < nn; j++){
	    const mxArray*yy = mxGetCell(y,j);
	    double*xr = mxGetPr(yy);
	    int mm = mxGetM(yy)*mxGetN(yy);
	    for (int k = 0; k < mm; k++){
	      *(dr+cnt+k) = *(xr+k);
	    }
	    cnt += mm;
	  }
 }
}
