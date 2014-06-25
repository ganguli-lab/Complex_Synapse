@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2013b
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2013b\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=BWbetaNloop_mex
set MEX_NAME=BWbetaNloop_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for BWbetaNloop > BWbetaNloop_mex.mki
echo COMPILER=%COMPILER%>> BWbetaNloop_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> BWbetaNloop_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> BWbetaNloop_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> BWbetaNloop_mex.mki
echo LINKER=%LINKER%>> BWbetaNloop_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> BWbetaNloop_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> BWbetaNloop_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> BWbetaNloop_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> BWbetaNloop_mex.mki
echo BORLAND=%BORLAND%>> BWbetaNloop_mex.mki
echo OMPFLAGS= >> BWbetaNloop_mex.mki
echo OMPLINKFLAGS= >> BWbetaNloop_mex.mki
echo EMC_COMPILER=msvc100>> BWbetaNloop_mex.mki
echo EMC_CONFIG=optim>> BWbetaNloop_mex.mki
"C:\Program Files\MATLAB\R2013b\bin\win64\gmake" -B -f BWbetaNloop_mex.mk
