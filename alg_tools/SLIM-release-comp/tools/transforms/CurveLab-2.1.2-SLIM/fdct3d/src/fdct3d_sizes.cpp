#include "commoninc.hpp"
#include "fdct3dinline.hpp"
#include "fdct3d_sizes.hpp"

int fdct3d_sizes( vector<int>& map, int& nbscales, int& nbdstz_coarse, int& ac,
      int& n1, int& n2, int& n3, size_t &totalcoeffs ){

      int L = nbscales;
      totalcoeffs=0L;
      unsigned int map_o=0;
      map.resize(3+(nbscales-2+ac)*10+(1-ac)*3);

      if(ac==1) {
	{
	  int s = 0;
	  double L1 = 4.0*n1/3.0 / pow2(L-1-s);	 double L2 = 4.0*n2/3.0 / pow2(L-1-s);	 double L3 = 4.0*n3/3.0 / pow2(L-1-s);
	  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;
	  fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
          map[map_o+0]=S1;
          map[map_o+1]=S2;
          map[map_o+2]=S3;
          map_o+=3;
	  totalcoeffs += S1*S2*S3;
	}
	for(int s=1; s<L; s++) {
	  double L1 = 4.0*n1/3.0 / pow2(L-1-s);	 double L2 = 4.0*n2/3.0 / pow2(L-1-s);	 double L3 = 4.0*n3/3.0 / pow2(L-1-s);
	  int nd = nbdstz_coarse * pow2(s/2);
          map[map_o+0]=nd*nd;
	  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;
	  fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
	  double W1 = L1/nd;  double W2 = L2/nd;  double W3 = L3/nd;
	  //face 0: x,y,z
	  for(int h=0; h<1; h++) { //(y first z second)
	    for(int g=0; g<1; g++) {
	      double xs = R1/4-(W1/2)/4;		double xe = R1;
	      double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
	      double zs = -R3 + (2*h-1)*W3/2;		double ze = -R3 + (2*h+3)*W3/2;
	      map[map_o+1] = int(ceil(xe-xs));
	      map[map_o+2] = int(ceil(ye-ys));
	      map[map_o+3] = int(ceil(ze-zs));
	    }
	  }
	  //face 1. y z x
	  for(int f=0; f<1; f++) {
	    for(int h=0; h<1; h++) {
	      double ys = R2/4-(W2/2)/4;		  double ye = R2;
	      double zs = -R3 + (2*h-1)*W3/2;		  double ze = -R3 + (2*h+3)*W3/2;
	      double xs = -R1 + (2*f-1)*W1/2;		  double xe = -R1 + (2*f+3)*W1/2;
	      map[map_o+4] = int(ceil(xe-xs));
	      map[map_o+5] = int(ceil(ye-ys));
	      map[map_o+6] = int(ceil(ze-zs));
	    }
	  }
	  //face 2. z x y
	  for(int g=0; g<1; g++) {
	    for(int f=0; f<1; f++) {
	      double zs = R3/4-(W3/2)/4;		double ze = R3;
	      double xs = -R1 + (2*f-1)*W1/2;		double xe = -R1 + (2*f+3)*W1/2;
	      double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
	      map[map_o+7] = int(ceil(xe-xs));
	      map[map_o+8] = int(ceil(ye-ys));
	      map[map_o+9] = int(ceil(ze-zs));
	    }
	  }
	  /* totalcoeffs += (cubedims[0]*cubedims[1]*cubedims[2]
                         +cubedims[3]*cubedims[4]*cubedims[5]
                         +cubedims[6]*cubedims[7]*cubedims[8])
                         *nd*nd*2;*/
          {
            size_t _cubedims_ = 0;
            _cubedims_ +=(size_t)(map[map_o+1]*map[map_o+2]*map[map_o+3]);
            _cubedims_ +=(size_t)(map[map_o+4]*map[map_o+5]*map[map_o+6]);
            _cubedims_ +=(size_t)(map[map_o+7]*map[map_o+8]*map[map_o+9]);
            size_t _multiply_ = (size_t)(nd*nd*2);
            totalcoeffs += _cubedims_*_multiply_;
          }
          map_o+=10;
	}
      } //end ac==1 case

      else { // ac==0 case
	
	{// wavelet scale [0]
	  int s = 0;
	  double L1 = 4.0*n1/3.0 / pow2(L-1-s);	 double L2 = 4.0*n2/3.0 / pow2(L-1-s);	 double L3 = 4.0*n3/3.0 / pow2(L-1-s);
	  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;
	  fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
          map[map_o+0]=S1;
          map[map_o+1]=S2;
          map[map_o+2]=S3;
          map_o+=3;
	  totalcoeffs += S1*S2*S3;
	}

	for(int s=1; s<L-1; s++) {// curvelet scales
	  double L1 = 4.0*n1/3.0 / pow2(L-1-s);	 double L2 = 4.0*n2/3.0 / pow2(L-1-s);	 double L3 = 4.0*n3/3.0 / pow2(L-1-s);
	  int nd = nbdstz_coarse * pow2(s/2);
          map[map_o+0]=nd*nd;
	  int S1, S2, S3;	 int F1, F2, F3;	 double R1, R2, R3;
	  fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
	  double W1 = L1/nd;  double W2 = L2/nd;  double W3 = L3/nd;
	  //face 0: x,y,z
	  for(int h=0; h<1; h++) { //(y first z second)
	    for(int g=0; g<1; g++) {
	      double xs = R1/4-(W1/2)/4;		double xe = R1;
	      double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
	      double zs = -R3 + (2*h-1)*W3/2;		double ze = -R3 + (2*h+3)*W3/2;
	      map[map_o+1] = int(ceil(xe-xs));
	      map[map_o+2] = int(ceil(ye-ys));
	      map[map_o+3] = int(ceil(ze-zs));
	    }
	  }
	  //face 1. y z x
	  for(int f=0; f<1; f++) {
	    for(int h=0; h<1; h++) {
	      double ys = R2/4-(W2/2)/4;		  double ye = R2;
	      double zs = -R3 + (2*h-1)*W3/2;		  double ze = -R3 + (2*h+3)*W3/2;
	      double xs = -R1 + (2*f-1)*W1/2;		  double xe = -R1 + (2*f+3)*W1/2;
	      map[map_o+4] = int(ceil(xe-xs));
	      map[map_o+5] = int(ceil(ye-ys));
	      map[map_o+6] = int(ceil(ze-zs));
	    }
	  }
	  //face 2. z x y
	  for(int g=0; g<1; g++) {
	    for(int f=0; f<1; f++) {
	      double zs = R3/4-(W3/2)/4;		double ze = R3;
	      double xs = -R1 + (2*f-1)*W1/2;		double xe = -R1 + (2*f+3)*W1/2;
	      double ys = -R2 + (2*g-1)*W2/2;		double ye = -R2 + (2*g+3)*W2/2;
	      map[map_o+7] = int(ceil(xe-xs));
	      map[map_o+8] = int(ceil(ye-ys));
	      map[map_o+9] = int(ceil(ze-zs));
	    }
	  }
	  /*totalcoeffs += (cubedims[0]*cubedims[1]*cubedims[2]
                         +cubedims[3]*cubedims[4]*cubedims[5]
                         +cubedims[6]*cubedims[7]*cubedims[8])
                         *nd*nd*2;*/
          {
            size_t _cubedims_ = 0;
            _cubedims_ +=(size_t)(map[map_o+1]*map[map_o+2]*map[map_o+3]);
            _cubedims_ +=(size_t)(map[map_o+4]*map[map_o+5]*map[map_o+6]);
            _cubedims_ +=(size_t)(map[map_o+7]*map[map_o+8]*map[map_o+9]);
            size_t _multiply_ = (size_t)(nd*nd*2);
            totalcoeffs += _cubedims_*_multiply_;
          }
          map_o+=10;
	}

        // wavelets scale [nbscales-1].
        map[map_o+0] = n1;
        map[map_o+1] = n2;
        map[map_o+2] = n3;
        map_o+=3;
        totalcoeffs += n1*n2*n3;
      }

      assert(map_o==map.size());
      return 0;
}
