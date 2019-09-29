/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
	Written by Lexing Ying
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

int main(int argc, char** argv)
{
  PetscInitialize(&argc,&argv,"options",NULL);  //PetscTruth flg = PETSC_FALSE;
  
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  iC( PetscPrintf(MPI_COMM_WORLD, "mpisize %d\n", mpisize) );  
  
  PetscTruth flg = PETSC_FALSE;
  int m;  iC( PetscOptionsGetInt("", "-m", &m, &flg) ); iA(flg==PETSC_TRUE);
  int n;  iC( PetscOptionsGetInt("", "-n", &n, &flg) ); iA(flg==PETSC_TRUE);
  int p;  iC( PetscOptionsGetInt("", "-p", &p, &flg) ); iA(flg==PETSC_TRUE);
  int b;  iC( PetscOptionsGetInt("", "-b", &b, &flg) ); iA(flg==PETSC_TRUE);
  int nbscales;  iC( PetscOptionsGetInt("", "-nbscales", &nbscales, &flg) ); iA(flg==PETSC_TRUE);
  int nbdstz_coarse;  iC( PetscOptionsGetInt("", "-nbdstz_coarse", &nbdstz_coarse, &flg) ); iA(flg==PETSC_TRUE);
  iC( PetscPrintf(MPI_COMM_WORLD, "-m %d\n", m) );
  iC( PetscPrintf(MPI_COMM_WORLD, "-n %d\n", n) );
  iC( PetscPrintf(MPI_COMM_WORLD, "-p %d\n", p) );
  iC( PetscPrintf(MPI_COMM_WORLD, "-b %d\n", b) );
  iC( PetscPrintf(MPI_COMM_WORLD, "-nbscales %d\n", nbscales) );
  iC( PetscPrintf(MPI_COMM_WORLD, "-nbdstz_coarse %d\n", nbdstz_coarse) );
  assert( m%b==0 && n%b==0 && p%b==0 );
  assert( (p/b)%mpisize==0 );
  
  double tm0, tm1, tma0, tma1;
  double cene; double wene; double yene; double eene;

  iC( PetscPrintf(MPI_COMM_WORLD, "\nW/ WAVELETS\n") );
  // W/ WAVELETS
  {
  CpxNumTnsBlkd X; CpxCrvletPrtd C;  CpxNumTnsBlkd W; CpxNumTnsBlkd Y;

  //1a. generate data
  srand48(0);
  BolNumTns newtnsexists(m/b,n/b,p/b);
  IntNumTns newtnsowners(m/b,n/b,p/b);
  iC( fdct3d_partition_cpxnumtnsblkd_z(m,n,p,b, newtnsexists, newtnsowners) );
  X.setup(m,n,p,b, newtnsowners);
  int e = X.e();  int f = X.f();  int g = X.g();
  for(int i=0; i<e; i++) for(int j=0; j<f; j++)	for(int k=0; k<g; k++) {
	 if(X.owners()(i,j,k)==mpirank) {
		CpxNumTns& Xblk = X.block(i,j,k);
		for(int ioff=0; ioff<b; ioff++) for(int joff=0; joff<b; joff++) for(int koff=0; koff<b; koff++) {
		  Xblk(ioff,joff,koff) = cpx(drand48(), 0); //cpx(drand48(), drand48());
		}
	 }
  }
  double xene = X.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "X energy %e\n", xene) );  //  if(mpirank==0)	 cerr<<"X energy  "<<xene<<endl;
  
  tma0 = MPI_Wtime();

  //2a. forward w/ wavelets
  tm0 = MPI_Wtime();
  iC( fdct3d_forward(m,n,p, nbscales,nbdstz_coarse, X, C,W) );
  tm1 = MPI_Wtime();  iC( PetscPrintf(MPI_COMM_WORLD, "FORWARD %d %d %d %e\n", p, b, mpisize, tm1-tm0) );

  cene = C.globalenergy();
  wene = W.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "CW energy %e %e %e\n", cene, wene, cene+wene) );
  
  //3a. inverse w/ wavelets
  tm0 = MPI_Wtime();
  iC( fdct3d_inverse(m,n,p, nbscales,nbdstz_coarse, C,W, Y) );
  tm1 = MPI_Wtime();  iC( PetscPrintf(MPI_COMM_WORLD, "INVERSE %d %d %d %e\n", p, b, mpisize, tm1-tm0) );

  yene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "Y energy %e\n", yene) );
  
  tma1 = MPI_Wtime();  iC( PetscPrintf(MPI_COMM_WORLD, "(FORWARD+INVERSE)/2 %d %d %d %e\n", p, b, mpisize, (tma1-tma0)/2.) );

  //4a. compute difference w/ wavelets
  for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
	 if(X.owners()(i,j,k)==mpirank) {
		iA(Y.owners()(i,j,k)==mpirank); 		//CHECK THEY HAVE THE SAME DISTRIBUTION
		CpxNumTns& Xblk = X.block(i,j,k);		CpxNumTns& Yblk = Y.block(i,j,k);
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
		  Yblk(ioff,joff,koff) -= Xblk(ioff,joff,koff);
		}
	 }
  }
  eene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "E energy %e\n", eene) );
  if(abs(eene)>0.00001*xene) iC( PetscPrintf(MPI_COMM_WORLD, "FATAL Energy Error\n") );
  }

  PetscFinalize(); exit(0);

  iC( PetscPrintf(MPI_COMM_WORLD, "\nW/O WAVELETS\n") );
  // W/O WAVELETS
  {
  CpxNumTnsBlkd X; CpxCrvletPrtd C;  CpxNumTnsBlkd Y;

  //1b. generate data
  srand48(0);
  BolNumTns newtnsexists(m/b,n/b,p/b);
  IntNumTns newtnsowners(m/b,n/b,p/b);
  iC( fdct3d_partition_cpxnumtnsblkd_z(m,n,p,b, newtnsexists, newtnsowners) );
  X.setup(m,n,p,b, newtnsowners);
  int e = X.e();  int f = X.f();  int g = X.g();
  for(int i=0; i<e; i++) for(int j=0; j<f; j++)	for(int k=0; k<g; k++) {
	 if(X.owners()(i,j,k)==mpirank) {
		CpxNumTns& Xblk = X.block(i,j,k);
		for(int ioff=0; ioff<b; ioff++) for(int joff=0; joff<b; joff++) for(int koff=0; koff<b; koff++) {
		  Xblk(ioff,joff,koff) = cpx(drand48(), 0); //cpx(drand48(), drand48());
		}
	 }
  }
  double xene = X.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "X energy %e\n", xene) );  //  if(mpirank==0)	 cerr<<"X energy  "<<xene<<endl;
  
  tma0 = MPI_Wtime();

  //2b. forward w/o wavelets
  tm0 = MPI_Wtime();
  iC( fdct3d_forward(m,n,p, nbscales,nbdstz_coarse, X, C) );
  tm1 = MPI_Wtime();  iC( PetscPrintf(MPI_COMM_WORLD, "FORWARD %d %d %d %e\n", p, b, mpisize, tm1-tm0) );

  cene = C.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "C energy %e\n", cene ) );
  
  //3b. inverse w/o wavelets
  tm0 = MPI_Wtime();
  iC( fdct3d_inverse(m,n,p,b, nbscales,nbdstz_coarse, C, Y) );
  tm1 = MPI_Wtime();  iC( PetscPrintf(MPI_COMM_WORLD, "INVERSE %d %d %d  %e\n", p, b, mpisize, tm1-tm0) );

  yene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "Y energy %e\n", yene) );
  
  tma1 = MPI_Wtime();  iC( PetscPrintf(MPI_COMM_WORLD, "(FORWARD+INVERSE)/2 %d %d %d %e\n", p, b, mpisize, (tma1-tma0)/2.) );

  //4b. compute difference w/o wavelets
  for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
	 if(X.owners()(i,j,k)==mpirank) {
		iA(Y.owners()(i,j,k)==mpirank); 		//CHECK THEY HAVE THE SAME DISTRIBUTION
		CpxNumTns& Xblk = X.block(i,j,k);		CpxNumTns& Yblk = Y.block(i,j,k);
		for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
		  Yblk(ioff,joff,koff) -= Xblk(ioff,joff,koff);
		}
	 }
  }
  eene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "E energy %e\n", eene) );
  if(abs(eene)>0.00001*xene) iC( PetscPrintf(MPI_COMM_WORLD, "FATAL Energy Error\n") );
  // iC( PetscSynchronizedPrintf(MPI_COMM_WORLD, "E energy %e\n", eene) );
  // iC( PetscSynchronizedFlush(MPI_COMM_WORLD) );
  }

  iC( PetscPrintf(MPI_COMM_WORLD, "\n") );
  PetscFinalize();
  return 0;
}
