include ../../makefile.mpi

LIB_SRC = 	cpxnumtnsblkd.cpp		cpxcrvletprtd.cpp		fdct3d.cpp	fdct3d_param.cpp	fdct3d_forward.cpp	fdct3d_inverse.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	test.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libfdct3d.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@
	@echo PetSC settings:
	@echo .. PETSC_DIR=${PETSC_DIR}
	@echo .. PETSC_ARCH=${PETSC_ARCH}

test: libfdct3d.a test.o
	${CXX} ${CXXFLAGS} -o test test.o libfdct3d.a ${LDFLAGS}

test1: libfdct3d.a test1.o
	${CXX} ${CXXFLAGS} -o test1 test1.o libfdct3d.a ${LDFLAGS}

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.out libfdct3d.a	test test1

tags:
	etags *.hpp *.cpp

####
fdct3d_forward.cpp: fdct3d_forward.inc
fdct3d_inverse.cpp: fdct3d_inverse.inc
