lib:
	@cd fdct_wrapping_cpp/src; make libfdct_wrapping.a ../../libsos/libdfdct_wrapping.so
	@cd fdct_usfft_cpp/src; make libfdct_usfft.a
	@cd fdct3d/src; make libfdct3d.a ../../libsos/libdfdct3d.so
	@cd fdct3d_outcore/src; make libfdct3d.a

test:
	@cd fdct_wrapping_cpp/src; make test
	@cd fdct_usfft_cpp/src; make test
	@cd fdct3d/src; make test
	@cd fdct3d_outcore/src; make test

matlab:
	@cd fdct_wrapping_cpp/src; make matlab
	@cd fdct_usfft_cpp/src; make matlab
	@cd fdct3d/src; make matlab
	@cd fdct_slim_matlab/src; make matlab

clean:
	@cd fdct_wrapping_cpp/src; make clean
	@cd fdct_usfft_cpp/src; make clean
	@cd fdct3d/src; make clean
	@cd fdct3d_outcore/src; make clean
	rm -rf `find . -type f -name "*mex.mex*"`
	rm -rf `find . -type f -name "lib*"`

mpi:
	@cd fdct3d_mpi/src; make libfdct3d.a

testmpi:
	@cd fdct3d_mpi/src; make test

cleanmpi:
	@cd fdct3d_mpi/src; make clean
