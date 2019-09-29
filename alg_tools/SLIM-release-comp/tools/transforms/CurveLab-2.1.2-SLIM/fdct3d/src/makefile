include ../../makefile.slim

LIB_SRC = 	fdct3d_forward.cpp	fdct3d_inverse.cpp	fdct3d_param.cpp \
                fdct3d_sizes.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	test.cpp

MEX_SRC = 	fdct3d_forward_mex.cpp fdct3d_inverse_mex.cpp fdct3d_param_mex.cpp \
                fdct3d_sizes_mex.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libfdct3d.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

../../libsos/libdfdct3d.so:	libfdct3d.a
	$(CXX) -Wfatal-errors --shared -o $@ -Wl,--whole-archive $^ -Wl,--no-whole-archive

test: libfdct3d.a test.o
	${CXX} ${CXXFLAGS} -o test test.o libfdct3d.a ${LDFLAGS}
	./test options

matlab: libfdct3d.a	${MEX_SRC}
	cp ${FFTW3_DIR}/lib/libfftw3.a .
	${RANLIB} libfftw3.a
	${MEX} fdct3d_forward_mex.cpp ${INCLUDES} libfdct3d.a libfftw3.a
	${MEX} fdct3d_inverse_mex.cpp ${INCLUDES} libfdct3d.a libfftw3.a
	${MEX} fdct3d_param_mex.cpp ${INCLUDES} libfdct3d.a libfftw3.a
	${MEX} fdct3d_sizes_mex.cpp ${INCLUDES} libfdct3d.a libfftw3.a
	rm libfftw3.a
	mv *_mex.mex* ../mex

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.out libfdct3d.a test

tags:
	etags *.hpp *.cpp
