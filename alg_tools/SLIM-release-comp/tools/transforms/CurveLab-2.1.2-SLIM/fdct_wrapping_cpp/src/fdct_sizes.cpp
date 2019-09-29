#include "fdct_wrapping.hpp"
#include "fdct_sizes.hpp"

using namespace fdct_wrapping_ns;

int fdct_sizes( vector<int>& map, int& nbscales, int& nbdstz_coarse, int& ac,
    int& n1, int& n2, size_t &totalcoeffs ){
    vector< vector<int> > nx, ny;
    {
      vector< vector<double> > sx, sy;
      vector< vector<double> > fx, fy;
      fdct_wrapping_param(n1, n2, nbscales, nbdstz_coarse, ac, sx, sy, fx, fy, nx, ny);
    }
    unsigned int map_o=0;
    totalcoeffs=0L;
    map.resize(2+(nbscales-2+ac)*5+(1-ac)*2);

    // sclae 0
    map[map_o+0]=nx[0][0];
    map[map_o+1]=ny[0][0];
    map_o=+2;
    totalcoeffs+=nx[0][0]*ny[0][0];

    for( int i=1; i<nbscales-1+ac; i++ ){//curvelet scales
        int nbangles_ = (int)(nbdstz_coarse * pow(2.,i/2))/4;//angles per quadrant
        map[map_o+0] = nbangles_;
        for( int q=0; q<4; q++) {
            if( q==0 || q==2 ){
                if (q==0) {
                    map[map_o+1] = nx[i][nbangles_*q];
                    map[map_o+2] = ny[i][nbangles_*q];
                } else {
                   assert(map[map_o+1] == nx[i][nbangles_*q]);
                   assert(map[map_o+2] == ny[i][nbangles_*q]);
                }
                totalcoeffs+=map[map_o+1]*map[map_o+2]*nbangles_;
            }
            else {
                if (q==1) {
                    map[map_o+3] = nx[i][nbangles_*q];
                    map[map_o+4] = ny[i][nbangles_*q];
                } else {
                   assert(map[map_o+3] == nx[i][nbangles_*q]);
                   assert(map[map_o+4] == ny[i][nbangles_*q]);
                }
                totalcoeffs+=map[map_o+3]*map[map_o+4]*nbangles_;
            }
        }
        map_o+=5;
    }

    if (!ac) { //wavelet scale
        map[map_o+0]=n1;
        map[map_o+1]=n2;
        map_o+=2;
        totalcoeffs+=n1*n2;
    }

    assert(map_o==map.size());
    return 0;
}
