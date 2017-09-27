@echo off
set MATLAB=C:\MATLAB\R2017a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\MATLAB\R2017a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=nn_batch_mex
set MEX_NAME=nn_batch_mex
set MEX_EXT=.mexw64
call setEnv.bat
echo # Make settings for nn_batch > nn_batch_mex.mki
echo COMPILER=%COMPILER%>> nn_batch_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> nn_batch_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> nn_batch_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> nn_batch_mex.mki
echo LINKER=%LINKER%>> nn_batch_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> nn_batch_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> nn_batch_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> nn_batch_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> nn_batch_mex.mki
echo BORLAND=%BORLAND%>> nn_batch_mex.mki
echo OMPFLAGS= >> nn_batch_mex.mki
echo OMPLINKFLAGS= >> nn_batch_mex.mki
echo EMC_COMPILER=mingw64>> nn_batch_mex.mki
echo EMC_CONFIG=optim>> nn_batch_mex.mki
"C:\MATLAB\R2017a\bin\win64\gmake" -B -f nn_batch_mex.mk
