/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

//THE HEADER FILE
#ifndef _FDCT_WRAPPING_HPP_
#define _FDCT_WRAPPING_HPP_

#include "fdct_wrapping_inc.hpp"
#include "numvec.hpp"
#include "nummat.hpp"
#include "offvec.hpp"
#include "offmat.hpp"
using std::vector;
using std::map;
using std::sort;
using std::pair;
using std::max;
using std::min;
using std::abs;

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

//New SLIM Functions (5 Total)

int fdct_wrapping_fft(CpxNumMat& x, CpxOffMat& O, bool inv);
int fdct_wrapping_skipfft(CpxNumMat& x, CpxOffMat& O, bool inv);
//this default function performs a fft or ifft of a N1 by N2 matrix depending on the flag inv.
//INPUTS:
// N1,N2 -- the size of the input image
// inv -- a bool which tells which output to produce and input to expect.
//OUTPUTS:
// O -- N1 by N2 matrix offset stored in a CpxOffMat class
// x -- N1 by N2 matrix stored in CpxNumMat class

int fdct_wrapping_fft2(CpxNumMat& x, CpxOffMat& O, bool inv);
//this default function performs a fft or ifft of a N1 by N2 matrix depending on the flag inv.  Usally it defaults to the fdct_wrapping_fft function. Unless otherwise specified.
//INPUTS:
// N1,N2 -- the size of the input image
// inv -- a bool which tells which output to produce and input to expect.
//OUTPUTS:
// O -- N1 by N2 matrix offset stored in a CpxOffMat class
// x -- N1 by N2 matrix stored in CpxNumMat class

int fftw_cleanmap();
//this function will clean up the MapPlan optimization.  The MapPlan optimization is not nessasary, and this function can just pass if no optimization is used, or anything global vars has to be cleaned.
//INPUTS:
// none
//OUTPUTS:
// should clean any global variables or destory any stored map plans.

int fdct_wrapping_clt(int N1, int N2, int nbscales, int nbangles_coarse, int allcurvelets, CpxOffMat& O, vector< vector<CpxNumMat> >& c);
//this function performs the forward curvelet transform. Allowing this to be access seperatly gives the user to perform there own FFT transfroamtion outside C++ and just run the curvlet information.
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  O -- N1 by N2 matrix stored in CpxOffMat class
//OUTPUTS:
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j) is accessed by c[s][w](i,j).

int fdct_wrapping_invclt(int N1, int N2, int nbscales, int nbangles_coarse, int allcurvelets, vector< vector<CpxNumMat> >& c, CpxOffMat& O);
//this function performs the inverse curvelet transform.  Allowing this to be access seperatly gives the user to perform there own FFT transfroamtion outside C++ and just run the curvlet information.
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j) is accessed by c[s][w](i,j). 
//OUTPUTS:
//  O -- N1 by N2 matrix stored in CpxOffMat class



//Original Functions

int fdct_wrapping(int N1, int N2, int nbscales, int nbangles_coarse, int ac, CpxNumMat& x, vector< vector<CpxNumMat> >& c);
int fdct_wrapping_nofft(int N1, int N2, int nbscales, int nbangles_coarse, int ac, CpxNumMat& x, vector< vector<CpxNumMat> >& c);
//this function performs the forward curvelet transform
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  x -- N1 by N2 matrix stored in CpxNumMat class
//OUTPUTS:
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j) is accessed by c[s][w](i,j). 

int ifdct_wrapping(int N1, int N2, int nbscales, int nbangles_coarse, int ac, vector< vector<CpxNumMat> >& c, CpxNumMat& x);
int ifdct_wrapping_nofft(int N1, int N2, int nbscales, int nbangles_coarse, int ac, vector< vector<CpxNumMat> >& c, CpxNumMat& x);
//this function performs the inverse curvelet transform
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j) is accessed by c[s][w](i,j). 
//OUTPUTS:
//  x -- N1 by N2 matrix stored in CpxNumMat class

int fdct_wrapping_param(int N1, int N2, int nbscales, int nbangles_coarse, int ac,
								vector< vector<double> >& sx, vector< vector<double> >& sy,
								vector< vector<double> >& fx, vector< vector<double> >& fy,
								vector< vector<int> >& nx, vector< vector<int> >& ny);
//this function obtains auxiliary information about curvelet transform 
//INPUTS:
//  N1,N2 -- the size of the input image
//  nbscales -- the total number of scales for subband decomposition
//  nbangles_coarse -- the number of angles in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//OUTPUTS:
//  sx -- for scale s and wege w, sx[s][w] is the first-cooridnate spacing of curvelet grid of scale s and wedge w
//  sy -- for scale s and wege w, sy[s][w] is the second-cooridnate spacing of curvelet grid of scale s and wedge w
//  fx -- for scale s and wege w, fx[s][w] is the first coordinate of the frequency center of curvelet of scale s and wedge w
//  fy -- for scale s and wege w, fy[s][w] is the second coordinate of the frequency center of curvelet of scale s and wedge w 
//  nx -- for scale s and wege w, nx[s][w] is the first dimension of the curvelet matrix of scale s and wedge w
//  ny -- for scale s and wege w, ny[s][w] is the second dimension of the curvelet matrix of scale s and wedge w 
//COMMENTS:
//  In computing sx and sy, we assume that the input grid are Nyquist sampling grid of domain [-1/2,1/2) by [-1/2,1/2) with
//  the center sampling point at origin (0,0). Therefore, for any scale s, and wedge w, nx[s][w] * sx[s][w] = 1 and
//  ny[s][w] * sy[s][w] = 1.

typedef pair<int,int> intpair;
typedef pair< intpair,vector<intpair> > CtrwInfo;
typedef vector< CtrwInfo > OnceInfo;

int fdct_collectnb(int s, int w, int nbscales, vector< vector<CpxOffMat> >& c, vector<intpair>& sws);
int fdct_grp(int N1, int N2, int nbscales, int nbangles_coarse, vector< vector<CpxOffMat> >& c, vector< OnceInfo >& grps);


FDCT_WRAPPING_NS_END_NAMESPACE

#endif
