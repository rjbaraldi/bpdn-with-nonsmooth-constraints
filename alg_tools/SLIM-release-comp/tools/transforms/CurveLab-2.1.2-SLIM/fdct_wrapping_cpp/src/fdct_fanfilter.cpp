#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

typedef pair<int,int> intpair;
typedef pair< intpair,vector<intpair> > CtrwInfo;
typedef vector< CtrwInfo > OnceInfo;

int fdct_grp(int N1, int N2, int nbscales, int nbangles_coarse, vector< vector<CpxOffMat> >& c, vector< OnceInfo >& grps)
{
  //1. simplest stuff
  /*
  for(int s=0; s<nbscales-1; s++)
	 for(int w=0; w<c[s].size(); w++) {
		vector<intpair> sws;
		fdct_collectnb(s, w, nbscales, c, sws);
		CtrwInfo ctrwinfo(intpair(s,w), sws);
		OnceInfo onceinfo;		onceinfo.push_back(ctrwinfo);
		grps.push_back( onceinfo );
		} */

  //2. more aggressive
  int s = 1;
  for(int w=0; w<c[s].size(); w++) {
	 vector<intpair> sws;
	 fdct_collectnb(s, w, nbscales, c, sws);
	 CtrwInfo ctrwinfo(intpair(s,w), sws);
	 OnceInfo onceinfo;		onceinfo.push_back(ctrwinfo);
	 grps.push_back( onceinfo );
  }

  for(int sg=0; sg<3; sg++)
	 for(int wg=0; wg<4; wg++) {
		OnceInfo onceinfo;
		for(int s=sg; s<nbscales-1; s=s+3)
		  for(int w=wg; w<c[s].size(); w=w+4) {
			 if(s!=1) {
				vector<intpair> sws;
				fdct_collectnb(s, w, nbscales, c, sws);
				CtrwInfo ctrwinfo(intpair(s,w), sws);
				onceinfo.push_back( ctrwinfo );
			 }
		  }
		grps.push_back( onceinfo);
	 }
  return 0;
}

int fdct_collectnb(int s, int w, int nbscales, vector< vector<CpxOffMat> >& c, vector<intpair>& sws)
{
  sws.clear();

  //search s-1
  if(s-1>=0) {
	 if(s-1==0) {
		sws.push_back( intpair(0,0) );
	 } else if(c[s].size()==c[s-1].size()) {
		int R = c[s-1].size();
		sws.push_back( intpair(s-1, (w-1+R)%R) );
		sws.push_back( intpair(s-1, w) );
		sws.push_back( intpair(s-1, (w+1+R)%R) );
	 } else {
		int R = c[s-1].size();
		if(w%2==0) {
		  sws.push_back( intpair(s-1, (w/2)%R) );
		  sws.push_back( intpair(s-1, (w/2-1+R)%R) );
		} else {
		  sws.push_back( intpair(s-1, (w/2)%R) );
		  sws.push_back( intpair(s-1, (w/2+1+R)%R) );
		}
	 }
  }
  //search s, same res
  if(s==0) {
	 sws.push_back( intpair(0,0) );
  } else {
	 int R = c[s].size();
	 sws.push_back( intpair(s, (w-1+R)%R) );
	 sws.push_back( intpair(s, w) );
	 sws.push_back( intpair(s, (w+1+R)%R) );
  }
  //search s+1, high res
  if(s+1<nbscales-1) {
	 if(s==0) {
		for(int k=0; k<c[1].size(); k++)
		  sws.push_back( intpair(1, k) );
	 } else if(c[s].size()==c[s+1].size()) {
		int R = c[s+1].size();
		sws.push_back( intpair(s+1, (w-1+R)%R) );
		sws.push_back( intpair(s+1, w) );
		sws.push_back( intpair(s+1, (w+1+R)%R) );				
	 } else {
		int R = c[s+1].size();
		sws.push_back( intpair(s+1, (2*w-1+R)%R) );
		sws.push_back( intpair(s+1, (2*w  +R)%R) );
		sws.push_back( intpair(s+1, (2*w+1+R)%R) );
		sws.push_back( intpair(s+1, (2*w+2+R)%R) );
	 }
  }

  return 0;
}


FDCT_WRAPPING_NS_END_NAMESPACE
