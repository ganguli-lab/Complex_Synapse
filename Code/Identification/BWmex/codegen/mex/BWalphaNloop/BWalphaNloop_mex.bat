@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2013b
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2013b\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=BWalphaNloop_mex
set MEX_NAME=BWalphaNloop_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for BWalphaNloop > BWalphaNloop_mex.mki
echo COMPILER=%COMPILER%>> BWalphaNloop_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> BWalphaNloop_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> BWalphaNloop_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> BWalphaNloop_mex.mki
echo LINKER=%LINKER%>> BWalphaNloop_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> BWalphaNloop_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> BWalphaNloop_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> BWalphaNloop_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> BWalphaNloop_mex.mki
echo BORLAND=%BORLAND%>> BWalphaNloop_mex.mki
echo OMPFLAGS= >> BWalphaNloop_mex.mki
echo OMPLINKFLAGS= >> BWalphaNloop_mex.mki
echo EMC_COMPILER=msvc100>> BWalphaNloop_mex.mki
echo EMC_CONFIG=optim>> BWalphaNloop_mex.mki
"C:\Program Files\MATLAB\R2013b\bin\win64\gmake" -B -f BWalphaNloop_mex.mk
