#include <cstdio>
#include <cmath>

#include "fdct_wrapping.hpp"
using namespace std;
using namespace fdct_wrapping_ns;

static bool dbg = false;

extern "C" {
#include <complex.h>
void jl_ifdct_wrapping_cpx(int, int, int, int, int, int, int, size_t, cpx*, cpx*);
}

void jl_ifdct_wrapping_cpx(int m, int n, int nbscales, int nbangles_coarse, int all_crvlts,
        int real_crvlts, int zero_finest,
        size_t cl, cpx *C, cpx *X) {
    vector< vector<CpxNumMat> > c;
    CpxNumMat x;
    if (dbg) printf("I am in inverse\n");
    if (dbg) printf("%d %d %d %d %d %d %d %ld\n",
        m,n,nbscales,nbangles_coarse,all_crvlts,real_crvlts,zero_finest,cl);
    vector< vector<int> > nx, ny;
    {
        vector< vector<double> > sx, sy;
        vector< vector<double> > fx, fy;
        fdct_wrapping_param(m, n, nbscales, nbangles_coarse, all_crvlts, sx, sy, fx, fy, nx, ny);
    }
    { // resize c
        int m_, n_;
        c.resize(nbscales);
        c[0].resize(1);
        m_=nx[0][0]; n_=ny[0][0];
        c[0][0].resize(m_,n_);
        for( int i=1; i<nbscales-1+all_crvlts; i++ ) {
            int nba = (int)(nbangles_coarse * pow(2.,i/2))/4;
            c[i].resize(nba*4);
            for( int q=0; q<4; q++)
                for( int j=0; j<nba; j++){
                    m_=nx[i][j+nba*q]; n_=ny[i][j+nba*q];
                    c[i][j+nba*q].resize(m_,n_);
                }
        }
        if (!all_crvlts) {
            c[nbscales-1].resize(1);
            c[nbscales-1][0].resize(m,n);
        }
    }
    { // copy C(C) to c(C++)
        int cnt = 0;
        int nbs = c.size();
        if (dbg) cout << nbs << endl;
        for(int ks=0; ks<nbs; ks++) {
            int nba = c[ks].size();
            if (dbg) cout << nba << " : ";
            for(int ka=0; ka<nba; ka++) {
                int bm = c[ks][ka].m();
                int bn = c[ks][ka].n();
                if (dbg) cout << bm << "," << bn << " ";
                for(int j=0; j<bn; j++)
                    for(int i=0; i<bm; i++) {
                        c[ks][ka](i,j) = C[cnt];
                        cnt++;
                    }
            }
            if (dbg) cout << endl;
        }
    }
    if (zero_finest) { // zero finest scales
        int nbs = c.size()-1;
        int nba = c[nbs].size();
        if (dbg) cout << nbs << " : zero finest scales" << endl;
        for(int ka=0; ka<nba; ka++) clear(c[nbs][ka]);
    }
    if (real_crvlts) { // convert from real
        int nbs = c.size();
        if (dbg) cout << nbs << endl;
        { // first scale
            if (dbg) cout << 1 << " : ";
            int bm = c[0][0].m();
            int bn = c[0][0].n();
            if (dbg) cout << bm << "," << bn << " ";
            for(int j=0; j<bn; j++)
                for(int i=0; i<bm; i++)
                    c[0][0](i,j) = real(c[0][0](i,j));
        }
        if (dbg) cout << endl;
        for(int kc=1; kc<nbs-1+all_crvlts; kc++) { // curvelet scales
            int nba = c[kc].size();
            int nba2 = nba/2;
            cpx II(0.,1.);
            if (dbg) cout << nba << " : ";
            for(int ka=0; ka<nba2; ka++) {
                double isqrt2 = 1./sqrt(2.);
                cpx dummy1, dummy2;
                int ka2 = ka+nba2;
                int bm = c[kc][ka].m();
                int bn = c[kc][ka].n();
                if (dbg) cout << isqrt2 << " ";
                if (dbg) cout << bm << "," << bn << " ";
                if (dbg) cout << c[kc][ka](0,0) << "," << c[kc][ka2](0,0) << "/";
                for(int j=0; j<bn; j++)
                    for(int i=0; i<bm; i++) {
                        dummy1 = c[kc][ka](i,j);
                        dummy2 = c[kc][ka2](i,j);
                        c[kc][ka](i,j) = isqrt2*(dummy1+dummy2*II);
                        c[kc][ka2](i,j) = isqrt2*(dummy1-dummy2*II);
                    }
                if (dbg) cout << c[kc][ka](0,0) << "," << c[kc][ka2](0,0) << endl;
            }
            if (dbg) cout << endl;
        }
        if (!all_crvlts) { // wavelet scale
            if (dbg) cout << 1 << " : ";
            int kw = nbs-1;
            int bm = c[kw][0].m();
            int bn = c[kw][0].n();
            if (dbg) cout << bm << "," << bn << " ";
            for(int j=0; j<bn; j++)
                for(int i=0; i<bm; i++)
                    c[kw][0](i,j) = real(c[kw][0](i,j));
        }
    }

    ifdct_wrapping(m, n, nbscales, nbangles_coarse, all_crvlts, c, x);

    { // copy x(C++) to X(C)
        int cnt = 0;
        for(int j=0; j<n; j++)
            for(int i=0; i<m; i++) {
                X[cnt] = x(i,j);
                cnt++;
            }
    }
}
