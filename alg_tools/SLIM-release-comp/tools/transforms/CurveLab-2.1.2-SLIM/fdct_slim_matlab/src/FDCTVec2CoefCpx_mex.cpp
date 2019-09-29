/*
 * =============================================================
 * FDCTVec2Coef.cpp
 *
 * AUTHOR
 * Gilles Hennenfent
 * University of British Columbia
 * ghennenfent@eos.ubc.ca
 * 
 * DEVELOPMENT
 * Feb 3rd, 2005: First version
 * Aug 30th, 2005: Updated to accommodate complex coefficients (Colin Russell)
 * =============================================================
 */

#include "mex.h"
#include "math.h"

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* 
 * prhs[0] is a vect full of C coef
 * prhs[1] is the header i.e. the cell array structure
 */

  int tmp;

 // Determine how many scale in FCT
 int n = mxGetN(prhs[1]);
 mxArray*res = mxCreateCellMatrix(1,n);

 double*vr = mxGetPr(prhs[0]);
 double*vi = mxGetPi(prhs[0]);

 int cnt = 0;
 for (int i = 0; i < n; i++) {
   mxArray*y = mxGetCell(prhs[1],i);
   // Determine # of angle for the ith scale
   int nn = mxGetN(y);
   mxSetCell(res, i, mxCreateCellMatrix(1,nn));
   mxArray*z = mxGetCell(res,i);
   for (int j = 0; j < nn; j++){
     mxArray*yy = mxGetCell(y,j);
     double*yr = mxGetPr(yy);
     mxSetCell(z, j, mxCreateDoubleMatrix((int)*yr,(int)*(yr+1),mxCOMPLEX));
     tmp = (int)*yr * (int)*(yr+1);
     mxArray*zz = mxGetCell(z,j);
     double*zr = mxGetPr(zz);
     double*zi = mxGetPi(zz);
     for (int k = 0; k < tmp; k++){
       *(zr+k) = *(vr+cnt+k);
       *(zi+k) = *(vi+cnt+k);
     }
     cnt += tmp;
   }
 }
 plhs[0] = res;
}
