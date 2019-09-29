include ../../makefile.opt

MEX_SRC = FDCTCoef2Vec_mex.cpp FDCTVec2Coef_mex.cpp FDCTCoef2VecCpx_mex.cpp FDCTVec2CoefCpx_mex.cpp

matlab:	${MEX_SRC}
	${MEX} FDCTCoef2Vec_mex.cpp
	${MEX} FDCTCoef2VecCpx_mex.cpp
	${MEX} FDCTVec2Coef_mex.cpp
	${MEX} FDCTVec2CoefCpx_mex.cpp
	mv *_mex.mex* ../mex

#------------------------------------------------------
clean:
	rm -rf *~ 

tags:
	etags *hpp *cpp
