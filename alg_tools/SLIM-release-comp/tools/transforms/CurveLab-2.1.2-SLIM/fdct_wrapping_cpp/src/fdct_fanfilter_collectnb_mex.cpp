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
  
  if(nrhs!=4)
	 mexErrMsgTxt("4 inputs required");
  if(nlhs!=1)
	 mexErrMsgTxt("1 outputs required");
  
  int s; mex2cpp(prhs[0], s);
  int w; mex2cpp(prhs[1], w);
  int nbscales; mex2cpp(prhs[2], nbscales);
  vector< vector< CpxOffMat> > c; mex2cpp(prhs[3], c);
  
  vector<intpair> sws;
  fdct_collectnb(s, w, nbscales, c, sws);
  
  cpp2mex(sws, plhs[0]);
  
  return;
}
