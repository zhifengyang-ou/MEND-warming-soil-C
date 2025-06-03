#!/bin/sh
FSRCS="src/MOD_Photosynthesis.F90 \
       src/MOD_MEND_TYPE.F90 \
       src/MOD_OPT_TYPE.F90 \
       src/MOD_USRFS.F90 \
       src/MOD_STRING.F90 \
       src/MOD_MEND.F90 \
       src/MOD_MCMC.F90 \
       src/MOD_OPT.F90 \
       src/MEND_IN.F90 \
       src/MEND_main.F90"

CPPFLAGS=''

echo $FSRCS

gfortran $FSRCS $CPPFLAGS -o biocon

./biocon

rm biocon
rm *.mod