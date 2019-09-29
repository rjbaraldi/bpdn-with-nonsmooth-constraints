/*
 * Added for w/o wavelets (hm)
 */

//------------------------------------------------------------------------------------
int fdct3d_forward(int N1, int N2, int N3, int nbscales, int nbdstz_coarse,
					 CpxNumTnsBlkd& X,
					 CpxCrvletPrtd& C)
{
  int ac=1; // curvelets at fine scales
  //check the size of x, make sure it is okay
  //check N1, N2, N3 are powers of 2, and proc number is power of 2 as well
  time_t tm0, tm1;  tm0 = time(NULL);
  
  int b = X.b();
  int e = X.e();  int f = X.f();  int g = X.g();
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  
  iC( MPI_Barrier(MPI_COMM_WORLD) );  //iC( PetscPrintf(MPI_COMM_SELF, "%d forward 0\n", mpirank) );  iC( MPI_Barrier(MPI_COMM_WORLD) );
  //-------------------------------------------
  //1. fft on W
  CpxNumTnsBlkd W = X;
  
  BolNumTns newtnsexists(e,f,g);
  IntNumTns newtnsowners(e,f,g);
  
  //scatter w to contain z slices 
  fdct3d_partition_cpxnumtnsblkd_z(N1,N2,N3,b, newtnsexists,newtnsowners);
  //hm if (mpirank==0) for(int E=0;E<e;E++) for(int F=0;F<f;F++) for(int G=0;G<g;G++) cerr << E << ' ' << F << ' ' << G << " -> " << newtnsexists(E,F,G) << "/" << newtnsowners(E,F,G) << endl;
  iC( W.scatter(newtnsexists) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "fwd w scatter %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //shift w's owner to z slices
  iC( W.shift(newtnsowners) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "fwd w shift %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //discard w's nonowners
  iC( W.discard() );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "fwd w discard %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  iC( fdct3d_fft(W) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "fwd w fft %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //-------------------------------------------
  //2. compute wedges
  int L = nbscales;
  
  //setup c 1,2,3, 6*np/8 processors are computing. 0 processor contains also the center wedge
  vector< vector<bool> > newcrvexists;
  vector< vector<int > > newcrvowners;
  iC( fdct3d_partition_cpxcrvletprtd(N1,N2,N3, nbscales,nbdstz_coarse,ac, newcrvexists, newcrvowners) );
  vector< vector<double> > fxs, fys, fzs;
  vector< vector<int   > > nxs, nys, nzs;
  iC( fdct3d_param(N1,N2,N3, nbscales,nbdstz_coarse,ac, fxs,fys,fzs, nxs,nys,nzs) );
  iC( C.setup(nxs, nys, nzs, newcrvowners) );
  
  //find out the required blocks from w for each processor
  iC( fdct3d_dependency(N1,N2,N3,b, nbscales,nbdstz_coarse,ac, newcrvowners, newtnsexists) );
  //scatter w according to c's request
  iC( W.scatter(newtnsexists) );  iC( W.check() );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "fwd w scatter %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //compute c (standard computation, since data is local already)
  {
	 int s = 0;
	 double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
	 //hm if(mpirank==0) cerr << L1 << ' ' << L2 << ' ' << L3 << ' ' << nxs[s].size() << endl;
	 fdct3d_forward_center(N1,N2,N3,b, L1,L2,L3, s, W, C);
  }
  for(int s=1; s<nbscales-1; s++) {
	 double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
	 //hm if(mpirank==0) cerr << L1 << ' ' << L2 << ' ' << L3 << ' ' << nxs[s].size() << endl;
	 int nd = nbdstz_coarse * pow2(s/2);
	 fdct3d_forward_angles(N1,N2,N3,b, L1,L2,L3, s, nd, W, C);
  }
  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "fwd c compute %f\n", difftime(tm1,tm0)) );  tm0 = tm1;	 

  return 0;
}


