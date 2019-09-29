#!/bin/csh -f

# test SLIM_COMP
if ( $?SLIM_COMP ) then
	echo SLIM_COMP = $SLIM_COMP
	test -e $SLIM_COMP/ibin/test_env4slim.csh || echo WARNING: cannot find myself in $SLIM_COMP
else
	echo FATAL ERROR: undefined environment SLIM_COMP || exit 1
endif

# test for matlab
which matlab >&/dev/null || echo ERROR: no matlab executable found
which mex >&/dev/null || echo ERROR: no mex executable found

# show MATLAB version
echo Checking MATLAB version
set matlab="matlab -nodesktop -nodisplay -nosplash"
cd $SLIM_COMP/tools/matlab_test/ || exit 1
$matlab -r slim_matlab_tests_version
