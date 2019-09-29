include ../../makefile.slim

LIB_SRC = fdct_wrapping.cpp ifdct_wrapping.cpp \
			fdct_wrapping_param.cpp fdct_fanfilter.cpp \
			fdct_sizes.cpp

LIB_OBJ = $(LIB_SRC:.cpp=.o)

JL_SRC	=	fdct_sizes_jl.cpp \
			fdct_wrapping_cpx_jl.cpp fdct_wrapping_real_jl.cpp \
			ifdct_wrapping_cpx_jl.cpp ifdct_wrapping_real_jl.cpp

JL_OBJ	= $(JL_SRC:.cpp=.o)

TST_SRC = test.cpp

MEX_SRC = fdct_wrapping_mex.cpp	ifdct_wrapping_mex.cpp \
			fdct_wrapping_param_mex.cpp	fdct_fanfilter_grp_mex.cpp	fdct_fanfilter_collectnb_mex.cpp  \
			fdct_sizes_mex.cpp

DEP     = $(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libfdct_wrapping.a: $(LIB_OBJ)
	$(AR) $(ARFLAGS) libfdct_wrapping.a $(LIB_OBJ)
	$(RANLIB) libfdct_wrapping.a

../../libsos/libdfdct_wrapping.so: libfdct_wrapping.a $(JL_OBJ)
	cp ${FFTW3_DIR}/lib/libfftw3.a .
	$(CXX) -Wfatal-errors --shared -o $@ $(JL_OBJ) libfdct_wrapping.a libfftw3.a
	rm libfftw3.a

test: libfdct_wrapping.a test.o
	${CXX} ${CXXFLAGS} -o test test.o libfdct_wrapping.a ${LDFLAGS}
	./test options

matlab: libfdct_wrapping.a ${MEX_SRC}
	cp ${FFTW3_DIR}/lib/libfftw3.a .
	${RANLIB} libfftw3.a
	${MEX} fdct_wrapping_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw3.a
	${MEX} fdct_wrapping_nofft_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw3.a
	${MEX} ifdct_wrapping_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw3.a
	${MEX} ifdct_wrapping_nofft_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw3.a
	${MEX} fdct_wrapping_param_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw3.a
	${MEX} fdct_fanfilter_grp_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw3.a
	${MEX} fdct_sizes_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw3.a
	#${MEX} fdct_fanfilter_collectnb_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw3.a
	rm libfftw3.a
	mv *_mex.mex* ../mex

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.out libfdct_wrapping.a test

tags:
	etags *hpp *cpp
