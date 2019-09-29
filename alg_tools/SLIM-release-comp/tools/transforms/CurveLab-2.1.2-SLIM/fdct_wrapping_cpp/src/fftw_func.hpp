/*
  New FFTW function written for Madagascar and using FFTW3.
*/

#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"
#include "fftw3.h"
#include <cstdio>
#include <cmath>

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

#define MAX(x,y) ((x) > (y) ? (x) : (y))

int fdct_wrapping_fft(CpxNumMat& x, CpxOffMat& O, bool inv)
{
  if(inv)
  {
    int n1 = O.m();		int n2 = O.n();
    x.resize(n1,n2);
    fftw_complex *in = (fftw_complex*)O.data();
    fftw_complex *out = (fftw_complex*)x.data();

    fftw_complex *workin, *workout;
     fftw_plan pl;
     int *p1, *p2, i, j, f1, f2, maxn;
     double sqrtprod;

     p1 = (int *)malloc(n1*sizeof(int));
     p2 = (int *)malloc(n2*sizeof(int));
     f1 = n1/2; f2 = n2/2;

     for(i=0; i<n1; i++)  p1[i] = (i-f1+n1)%n1;
     for(i=0; i<n2; i++)  p2[i] = (i-f2+n2)%n2;

     maxn = MAX(n1,n2);
     workin = (fftw_complex *)fftw_malloc(maxn*sizeof(fftw_complex));
     workout = (fftw_complex *)fftw_malloc(maxn*sizeof(fftw_complex));

     pl = fftw_plan_dft_1d(n1, workin, workout, FFTW_BACKWARD, FFTW_ESTIMATE);
         for(j=0; j<n2; j++) {
             for(i=0; i<n1; i++) {
                 workin[p1[i]][0] = in[j*n1+i][0];
                 workin[p1[i]][1] = in[j*n1+i][1];
             }
             fftw_execute(pl);
             for(i=0; i<n1; i++) {
                 out[j*n1+i][0] = workout[i][0];
                 out[j*n1+i][1] = workout[i][1];
         }
     }
     fftw_destroy_plan(pl);  
     sqrtprod = 1.0/sqrt((double)(n1*n2));

     pl = fftw_plan_dft_1d(n2, workin, workout, FFTW_BACKWARD, FFTW_ESTIMATE);
         for(i=0; i<n1; i++) {
             for(j=0; j<n2; j++) {
                 workin[p2[j]][0] = out[j*n1+i][0];
                 workin[p2[j]][1] = out[j*n1+i][1];
             }
             fftw_execute(pl);
             for(j=0; j<n2; j++) {
                 out[j*n1+i][0] = workout[j][0]*sqrtprod;
                 out[j*n1+i][1] = workout[j][1]*sqrtprod;
             }
         }
     fftw_destroy_plan(pl);  

    //fftw_cleanup_threads();
    free(p1);
    free(p2);
    fftw_free(workin);
    fftw_free(workout);
  } else
  {	
    int n1 = x.m();		int n2 = x.n();
    O.resize(n1, n2);
    fftw_complex *in = (fftw_complex*)x.data();
    fftw_complex *out = (fftw_complex*)O.data();

     fftw_complex *workin, *workout;
     fftw_plan pl;
     int *p1, *p2, i, j, f1, f2, maxn;
     double sqrtprod;

     p1 = (int *)malloc(n1*sizeof(int));
     p2 = (int *)malloc(n2*sizeof(int));
     f1 = n1/2; f2 = n2/2;

     for(i=0; i<n1; i++)  p1[i] = (i-f1+n1)%n1;
     for(i=0; i<n2; i++)  p2[i] = (i-f2+n2)%n2;

     maxn = MAX(n1,n2);
     workin = (fftw_complex *)fftw_malloc(maxn*sizeof(fftw_complex));
     workout = (fftw_complex *)fftw_malloc(maxn*sizeof(fftw_complex));

     pl = fftw_plan_dft_1d(n1, workin, workout, FFTW_FORWARD, FFTW_ESTIMATE);
         for(j=0; j<n2; j++) {
             for(i=0; i<n1; i++) {
                 workin[i][0] = in[j*n1+i][0];
                 workin[i][1] = in[j*n1+i][1];
             }
             fftw_execute(pl);
             for(i=0; i<n1; i++) {
                 out[j*n1+i][0] = workout[p1[i]][0];
                 out[j*n1+i][1] = workout[p1[i]][1];
             }
     }
     fftw_destroy_plan(pl);  

     sqrtprod = 1.0/sqrt((double)(n1*n2));
     pl = fftw_plan_dft_1d(n2, workin, workout, FFTW_FORWARD, FFTW_ESTIMATE);
         for(i=0; i<n1; i++) {
             for(j=0; j<n2; j++) {
                 workin[j][0] = out[j*n1+i][0];
                 workin[j][1] = out[j*n1+i][1];
             }
             fftw_execute(pl);
             for(j=0; j<n2; j++) {
                 out[j*n1+i][0] = workout[p2[j]][0]*sqrtprod;
                 out[j*n1+i][1] = workout[p2[j]][1]*sqrtprod;
             }
         }
    fftw_destroy_plan(pl);  

    free(p1);
    free(p2);
    fftw_free(workin);
    fftw_free(workout);
  }  
  return 0;
}

int fdct_wrapping_skipfft(CpxNumMat& x, CpxOffMat& O, bool inv)
{
  // std::cout << "fft copy start " << inv << "\n"; // debug
  if(inv)
  {
    int n1 = O.m();		int n2 = O.n();
    x.resize(n1,n2);
    fftw_complex *in = (fftw_complex*)O.data();
    fftw_complex *out = (fftw_complex*)x.data();

     int i, j;

     for(j=0; j<n2; j++) {
         for(i=0; i<n1; i++) {
             // if (i>13 && i<19 && j>13 && j<19) // debug
                // printf("inv in: %d %d %.4f %.4f \n",j,i,in[j*n1+i][0],in[j*n1+i][1]);
             out[j*n1+i][0] = in[j*n1+i][0];
             out[j*n1+i][1] = in[j*n1+i][1];
         }
     }

    //fftw_cleanup_threads();
  } else
  {	
    int n1 = x.m();		int n2 = x.n();
    O.resize(n1, n2);
    fftw_complex *in = (fftw_complex*)x.data();
    fftw_complex *out = (fftw_complex*)O.data();

     int i, j;

     for(i=0; i<n1; i++) {
         for(j=0; j<n2; j++) {
             out[j*n1+i][0] = in[j*n1+i][0];
             out[j*n1+i][1] = in[j*n1+i][1];
             // if (i>13 && i<19 && j>13 && j<19) // debug
                // printf("for out: %d %d %.4f %.4f \n",j,i,out[j*n1+i][0],out[j*n1+i][1]);
         }
     }

    // std::cout << "fft copy end\n"; // debug
  }  
  return 0;
}

int fdct_wrapping_fft2(CpxNumMat& x, CpxOffMat& O, bool inv)
{
  //Default to first function.
  fdct_wrapping_fft(x, O, inv);
  return 0;
}

int fftw_cleanmap()
{
  //There is no PlanMap oiptimization used in this code. Just pass.
  return 0;
}

FDCT_WRAPPING_NS_END_NAMESPACE
