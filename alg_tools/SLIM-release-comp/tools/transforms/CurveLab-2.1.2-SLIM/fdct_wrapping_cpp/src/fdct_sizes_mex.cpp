/*
   Copyright (C) 2010 SLIM
	Written by Nameet Kumar
*/

#include "mex.h"
#include "matrix.h"

#include "fdct_wrapping.hpp"
#include "fdct_sizes.hpp"
using namespace fdct_wrapping_ns;

#include "mexaux.hpp"

extern void _main();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if(nrhs!=5)
	 mexErrMsgTxt("5 inputs required");
  //if(nlhs!=2 && nlhs != 1)
	// mexErrMsgTxt("1 or 2 outputs required");
  
  int m; mex2cpp(prhs[0], m);
  int n; mex2cpp(prhs[1], n);
  int nbscales; mex2cpp(prhs[2], nbscales);
  int nbdstz_coarse; mex2cpp(prhs[3], nbdstz_coarse);
  int ac; mex2cpp(prhs[4], ac);
  
  
  vector<int> map;  //vector<int> extra;
  size_t totalcoeffs;
  
  fdct_sizes( map, nbscales, nbdstz_coarse, ac, m, n,  totalcoeffs);
  
  cpp2mex(map, plhs[0]);
  if(nlhs != 1)
  plhs[1] = mxCreateDoubleScalar(long(totalcoeffs));
  
  return;
}

