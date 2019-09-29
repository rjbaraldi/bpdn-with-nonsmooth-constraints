/*
  New FFTW function written for Madagascar and using FFTW3.
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"
#include "fftw3.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))

int fdct3d_fft_base(fftw_complex *in, fftw_complex *out, int n1, int n2, int n3, bool inv);

int fdct3d_fft(CpxNumTns& x, CpxOffTns& O, bool inv)
{
  if(inv)
  {
    x.resize(O.m(),O.n(),O.p());
    fdct3d_fft_base((fftw_complex*)O.data(),(fftw_complex*)x.data(),O.m(),O.n(),O.p(),inv);
  } else
  {	
    O.resize(x.m(),x.n(),x.p());
    fdct3d_fft_base((fftw_complex*)x.data(),(fftw_complex*)O.data(),x.m(),x.n(),x.p(),inv);
  }  
  return 0;
}

int fdct3d_fft2(CpxNumTns& x, CpxOffTns& O, bool inv)
{
  //Default to first function.
  fdct3d_fft(x, O, inv);
  return 0;
}

int fdct3d_fft_base(fftw_complex *in, fftw_complex *out, int n1, int n2, int n3, bool inv)
{
  if(inv)
  {
     fftw_complex *workin, *workout;
     fftw_plan pl;
     int *p1, *p2, *p3, i, j, k, f1, f2, f3, maxn;
     double sqrtprod;

     p1 = (int *)malloc(n1*sizeof(int));
     p2 = (int *)malloc(n2*sizeof(int));
     p3 = (int *)malloc(n3*sizeof(int));
     f1 = n1/2; f2 = n2/2; f3 = n3/2;

     for(i=0; i<n1; i++)  p1[i] = (i-f1+n1)%n1;
     for(i=0; i<n2; i++)  p2[i] = (i-f2+n2)%n2;
     for(i=0; i<n3; i++)  p3[i] = (i-f3+n3)%n3;

     maxn = MAX(n1,MAX(n2,n3));
     workin = (fftw_complex *)fftw_malloc(maxn*sizeof(fftw_complex));
     workout = (fftw_complex *)fftw_malloc(maxn*sizeof(fftw_complex));

     //fftw_plan_with_nthreads(2);
     pl = fftw_plan_dft_1d(n1, workin, workout, FFTW_BACKWARD, FFTW_ESTIMATE);
     //pl = fftw_create_plan(n1,  FFTW_BACKWARD, FFTW_ESTIMATE);
     for(k=0; k<n3; k++) {
         for(j=0; j<n2; j++) {
             for(i=0; i<n1; i++) {
                 workin[p1[i]][0] = in[k*n2*n1+j*n1+i][0];
                 workin[p1[i]][1] = in[k*n2*n1+j*n1+i][1];
             }
             fftw_execute(pl);
             for(i=0; i<n1; i++) {
                 out[k*n2*n1+j*n1+i][0] = workout[p1[i]][0];
                 out[k*n2*n1+j*n1+i][1] = workout[p1[i]][1];
             }
         }
     }
     fftw_destroy_plan(pl);  

     pl = fftw_plan_dft_1d(n2, workin, workout, FFTW_BACKWARD, FFTW_ESTIMATE);
     //pl = fftw_create_plan(n2,  FFTW_BACKWARD, FFTW_ESTIMATE);
     for(k=0; k<n3; k++) {
         for(i=0; i<n1; i++) {
             for(j=0; j<n2; j++) {
                 workin[p2[j]][0] = out[k*n2*n1+j*n1+i][0];
                 workin[p2[j]][1] = out[k*n2*n1+j*n1+i][1];
             }
             fftw_execute(pl);
             for(j=0; j<n2; j++) {
                 out[k*n2*n1+j*n1+i][0] = workout[p2[j]][0];
                 out[k*n2*n1+j*n1+i][1] = workout[p2[j]][1];
             }
         }
     }
     fftw_destroy_plan(pl);  

     sqrtprod = sqrt((double)(n1*n2*n3));
     sqrtprod = 1.0/sqrt((double)(n1*n2*n3));

     pl = fftw_plan_dft_1d(n3, workin, workout, FFTW_BACKWARD, FFTW_ESTIMATE);
     //pl = fftw_create_plan(n3,  FFTW_BACKWARD, FFTW_ESTIMATE);
      for(j=0; j<n2; j++) {
         for(i=0; i<n1; i++) {
             for(k=0; k<n3; k++) {
                 workin[p3[k]][0] = out[k*n2*n1+j*n1+i][0];
                 workin[p3[k]][1] = out[k*n2*n1+j*n1+i][1];
             }
             fftw_execute(pl);
             for(k=0; k<n3; k++) {
                 out[k*n2*n1+j*n1+i][0] = workout[p3[k]][0]*sqrtprod;
                 out[k*n2*n1+j*n1+i][1] = workout[p3[k]][1]*sqrtprod;
             }
         }
     }
     fftw_destroy_plan(pl);  

    //fftw_cleanup_threads();
    free(p1);
    free(p2);
    free(p3);
    fftw_free(workin);
    fftw_free(workout);
  } else
  {	
     fftw_complex *workin, *workout;
     fftw_plan pl;
     int *p1, *p2, *p3, i, j, k, f1, f2, f3, maxn;
     double sqrtprod;

     p1 = (int *)malloc(n1*sizeof(int));
     p2 = (int *)malloc(n2*sizeof(int));
     p3 = (int *)malloc(n3*sizeof(int));
     f1 = n1/2; f2 = n2/2; f3 = n3/2;

     for(i=0; i<n1; i++)  p1[i] = (i-f1+n1)%n1;
     for(i=0; i<n2; i++)  p2[i] = (i-f2+n2)%n2;
     for(i=0; i<n3; i++)  p3[i] = (i-f3+n3)%n3;

     maxn = MAX(n1,MAX(n2,n3));
     workin = (fftw_complex *)fftw_malloc(maxn*sizeof(fftw_complex));
     workout = (fftw_complex *)fftw_malloc(maxn*sizeof(fftw_complex));

     pl = fftw_plan_dft_1d(n1, workin, workout, FFTW_FORWARD, FFTW_ESTIMATE);
     //pl = fftw_create_plan(n1, FFTW_FORWARD, FFTW_ESTIMATE);
     for(k=0; k<n3; k++) {
         for(j=0; j<n2; j++) {
             for(i=0; i<n1; i++) {
                 workin[p1[i]][0] = in[k*n2*n1+j*n1+i][0];
                 workin[p1[i]][1] = in[k*n2*n1+j*n1+i][1];
             }
             fftw_execute(pl);
             for(i=0; i<n1; i++) {
                 out[k*n2*n1+j*n1+i][0] = workout[p1[i]][0];
                 out[k*n2*n1+j*n1+i][1] = workout[p1[i]][1];
             }
         }
     }
     fftw_destroy_plan(pl);  

     pl = fftw_plan_dft_1d(n2, workin, workout, FFTW_FORWARD, FFTW_ESTIMATE);
     //pl = fftw_create_plan(n2, FFTW_FORWARD, FFTW_ESTIMATE);
     for(k=0; k<n3; k++) {
         for(i=0; i<n1; i++) {
             for(j=0; j<n2; j++) {
                 workin[p2[j]][0] = out[k*n2*n1+j*n1+i][0];
                 workin[p2[j]][1] = out[k*n2*n1+j*n1+i][1];
             }
             fftw_execute(pl);
             for(j=0; j<n2; j++) {
                 out[k*n2*n1+j*n1+i][0] = workout[p2[j]][0];
                 out[k*n2*n1+j*n1+i][1] = workout[p2[j]][1];
             }
         }
     }
     fftw_destroy_plan(pl);  

     sqrtprod = 1.0/sqrt((double)(n1*n2*n3));

     pl = fftw_plan_dft_1d(n3, workin, workout, FFTW_FORWARD, FFTW_ESTIMATE);
     //pl = fftw_create_plan(n3, FFTW_FORWARD, FFTW_ESTIMATE);
      for(j=0; j<n2; j++) {
         for(i=0; i<n1; i++) {
             for(k=0; k<n3; k++) {
                 workin[p3[k]][0] = out[k*n2*n1+j*n1+i][0];
                 workin[p3[k]][1] = out[k*n2*n1+j*n1+i][1];
             }
             fftw_execute(pl);
             for(k=0; k<n3; k++) {
                 out[k*n2*n1+j*n1+i][0] = workout[p3[k]][0]*sqrtprod;
                 out[k*n2*n1+j*n1+i][1] = workout[p3[k]][1]*sqrtprod;
             }
         }
     }
    fftw_destroy_plan(pl);  

    free(p1);
    free(p2);
    free(p3);
    fftw_free(workin);
    fftw_free(workout);
  }  
  return 0;
}

int fftw_cleanmap()
{
  //There is no PlanMap oiptimization used in this code. Just pass.
  return 0;
}

