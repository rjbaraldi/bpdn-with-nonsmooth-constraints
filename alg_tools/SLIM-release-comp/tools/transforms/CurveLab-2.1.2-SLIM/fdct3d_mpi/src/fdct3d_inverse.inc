/*
 * Added for w/o wavelets (hm)
 */

//------------------------------------------------------------------------------------
int fdct3d_inverse(int N1, int N2, int N3, int b, int nbscales, int nbdstz_coarse,
					 CpxCrvletPrtd& C,
					 CpxNumTnsBlkd& X)
{
  int ac=1; // curvelets at fine scales
  //check the size of c,w, make sure it is okay
  time_t tm0, tm1;  tm0 = time(NULL);
  
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  
  iC( MPI_Barrier(MPI_COMM_WORLD) );  //iC( PetscPrintf(MPI_COMM_WORLD, "%d inverse 0\n", mpirank) );  iC( MPI_Barrier(MPI_COMM_WORLD) );
  //-------------------------------------------
  //1. fft on X
  // X = W;
  BolNumTns newtnsexists(N1/b,N2/b,N3/b);
  IntNumTns newtnsowners(N1/b,N2/b,N3/b);
  fdct3d_partition_cpxnumtnsblkd_z(N1,N2,N3,b, newtnsexists,newtnsowners);
  X.setup(N1,N2,N3,b, newtnsowners);
  int e = X.e();  int f = X.f();  int g = X.g();
  
  fdct3d_partition_cpxnumtnsblkd_z(N1,N2,N3,b, newtnsexists,newtnsowners);
  //scatter x to contain z slices 
  iC( X.scatter(newtnsexists) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "inv x scatter %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //shift x's owner to z slices
  iC( X.shift(newtnsowners) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "inv x shift %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //discard x's nonowners
  iC( X.discard() );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "inv x discard %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //ifft
  iC( fdct3d_fft(X) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "inv x fft %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //scale x with POU
  DblOffVec big1(N1);  fdct3d_lowpass(2.0*N1/3, big1);
  DblOffVec big2(N2);  fdct3d_lowpass(2.0*N2/3, big2);
  DblOffVec big3(N3);  fdct3d_lowpass(2.0*N3/3, big3);
  IntNumTns& Xowners = X.owners();
  for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
	 if(Xowners(i,j,k)==mpirank) {
		CpxNumTns& Xblk = X.block(i,j,k);
		int istt = i*b-N1/2;		int jstt = j*b-N2/2;		int kstt = k*b-N3/2;
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
		  double pou = big1(ioff+istt) * big2(joff+jstt)*big3(koff+kstt);
		  Xblk(ioff, joff, koff) *= sqrt(1-pou*pou);
		}
	 }
  }
  
  //-------------------------------------------
  //2. compute wedges
  int L = nbscales;
  //setup c, 1,2,3, 6*np/8 processors are computing. 0 processor contains also the center wedge
  vector< vector<bool> > newcrvexists;
  vector< vector<int > > newcrvowners;
  fdct3d_partition_cpxcrvletprtd(N1,N2,N3, nbscales,nbdstz_coarse,ac, newcrvexists, newcrvowners);
  //vector< vector<double> > fxs, fys, fzs;
  //vector< vector<int   > > nxs, nys, nzs;
  //fdct3d_param(N1,N2,N3, nbscales,nbdstz_coarse, fxs,fys,fzs, nxs,nys,nzs); //LEXING: might not be necessary
  
  //find out the required blocks from x for each processor
  fdct3d_dependency(N1,N2,N3,b, nbscales,nbdstz_coarse,ac, newcrvowners, newtnsexists);
  
  //expand x according to c's request
  iC( X.expand(newtnsexists) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "inv x expand %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //compute contribution to x from c
  {
	 int s = 0;
	 double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
	 fdct3d_inverse_center(N1,N2,N3,b, L1,L2,L3, s, C, X);
  }
  
  for(int s=1; s<nbscales-1; s++) {
	 double L1 = 4.0*N1/3.0 / pow2(L-1-s);	 double L2 = 4.0*N2/3.0 / pow2(L-1-s);	 double L3 = 4.0*N3/3.0 / pow2(L-1-s);
	 int nd = nbdstz_coarse * pow2(s/2);
	 fdct3d_inverse_angles(N1,N2,N3,b, L1,L2,L3, s, nd, C, X);
  }
  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "inv c compute %f\n", difftime(tm1,tm0)) );  tm0 = tm1;

  //combine x
  iC( X.combine() );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "inv x combine %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  
  //-------------------------------------------
  //3. ifft on X
  iC( fdct3d_ifft(X) );  iC( MPI_Barrier(MPI_COMM_WORLD) );  tm1 = time(NULL);  //iC( PetscPrintf(MPI_COMM_WORLD, "inv x ifft %f\n", difftime(tm1,tm0)) );  tm0 = tm1;
  //done
  { // Energy check
    //hm double xene = X.globalenergy(); double cene = C.globalenergy();
    //hm if (mpirank==0) cerr << "\tenergy here " << xene << ' ' << cene << endl;
  }
  
  return 0;
}

