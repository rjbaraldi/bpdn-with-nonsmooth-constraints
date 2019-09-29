/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#include "mex.h"
#include "matrix.h"

#include "fdct_wrapping.hpp"

#include "mexaux.hpp"

using namespace std;
using namespace fdct_wrapping_ns;

//digital curvelet transform
extern void _main();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  typedef pair<int,int> intpair;
  typedef pair< intpair,vector<intpair> > CtrwInfo;
  typedef vector< CtrwInfo > OnceInfo;
  
  if(nrhs!=5)
	 mexErrMsgTxt("5 inputs required");
  if(nlhs!=1)
	 mexErrMsgTxt("1 outputs required");
  
  int m; mex2cpp(prhs[0], m);
  int n; mex2cpp(prhs[1], n);
  int nbscales; mex2cpp(prhs[2], nbscales);
  int nbangles_coarse; mex2cpp(prhs[3], nbangles_coarse);
  vector< vector< CpxOffMat> > c; mex2cpp(prhs[4], c);
  
  vector< OnceInfo > grps;
  fdct_grp(m, n, nbscales, nbangles_coarse, c, grps);
  
  cpp2mex(grps, plhs[0]);
  
  return;
}
