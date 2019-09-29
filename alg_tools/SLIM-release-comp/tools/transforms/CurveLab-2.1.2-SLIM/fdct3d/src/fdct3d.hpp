/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/


#ifndef _FDCT3D_HPP_
#define _FDCT3D_HPP_

#include "commoninc.hpp"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"
#include "offvec.hpp"
#include "offmat.hpp"
#include "offtns.hpp"

//New SLIM Functions (5 Total)

int fdct3d_fft(CpxNumTns& x, CpxOffTns& O, bool inv);
//this default function performs a fft or ifft of a N1 by N2 matrix depending on the flag inv.
//INPUTS:
// N1,N2 -- the size of the input image
// inv -- a bool which tells which output to produce and input to expect.
//OUTPUTS:
// O -- N1 by N2 matrix offset stored in a CpxOffMat class
// x -- N1 by N2 matrix stored in CpxNumMat class

int fdct3d_fft2(CpxNumTns& x, CpxOffTns& O, bool inv);
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

int fdct3d_forward_clt(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac, CpxOffTns& O, vector< vector<CpxNumTns> >& c);
int fdct3d_inverse_clt(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac, vector< vector<CpxNumTns> >& c, CpxOffTns& O);
//these fuctions are just partial wrappers of the fdct3d_forward and inverse functions.  They allow the user to access the curvelet transformation without the initial fourier transformation.



//Original Functions

int fdct3d_forward(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac, CpxNumTns& x, vector< vector<CpxNumTns> >& c);
//this function performs the forward curvelet transform
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  x -- N1 by N2 by N3 tensor stored in CpxNumTns class
//OUTPUTS:
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j,k) is accessed by c[s][w](i,j,k).


int fdct3d_inverse(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac, vector< vector<CpxNumTns> >& c, CpxNumTns& x);
//this function performs the inverse curvelet transform
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//  c -- the curvelet coefficients data structure coeffient at scale s, wedge w and indices (i,j,k) is accessed by c[s][w](i,j,k). 
//OUTPUTS:
//  x -- N1 by N2 by N3 tensor stored in CpxNumTns class


int fdct3d_param(  int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac,
						 vector< vector<double> >& fxs, vector< vector<double> >& fys, vector< vector<double> >& fzs,
						 vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs);
//this function obtains auxiliary information about curvelet transform 
//INPUTS:
//  N1,N2,N3 -- the size of the input data
//  nbscales -- the total number of scales for subband decomposition
//  nbdstz_coarse -- the number of discretizations in each direction on each face in the 2nd coarest scale
//  ac -- ac==1 use curvelet at the finest level, ac==0, use wavelet at the finest level 
//OUTPUTS:
//  fx,fy,fz -- for scale s and wege w, fx[s][w], fy[s][w] and fz[s][w] give the coordinate of the frequency center of curvelet of scale s and wedge w
//  nx,ny,nz -- for scale s and wege w, nx[s][w], ny[s][w] and nz[s][w] give the dimensions of the curvelet matrix of scale s and wedge w


#endif
