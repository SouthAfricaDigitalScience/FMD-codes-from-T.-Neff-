#!/bin/bash -e
# this should be run after check-build finishes.
. /etc/profile.d/modules.sh
module load deploy
module add ncurses
module add gcc/5.2.0
module add openmpi/1.8.8-gcc-${GCC_VERSION}
module add lapack/3.6.0-gcc-${GCC_VERSION}

echo "All tests have passed, will now build into ${SOFT_DIR}"
echo $PWD
cd $WORKSPACE
make veryclean
make
export FMD=${SOFT_DIR}
make install
echo "Creating the modules file directory ${LIBRARIES_MODULES}"
mkdir -p ${PHYSICAL_MODULES}/${NAME}
(
cat <<MODULE_FILE
#%Module1.0
## $NAME modulefile
##
proc ModulesHelp { } {
    puts stderr "       This module does nothing but alert the user"
    puts stderr "       that the [module-info name] module is not available"
}

module add ncurses
module add gcc/5.2.0
module add openmpi/1.8.8-gcc-5.2.0
module add lapack/3.6.0-gcc-5.2.0

module-whatis   "$NAME $VERSION : See https://github.com/SouthAfricaDigitalScience/FMD-deploy"
setenv FMD_VERSION       $VERSION
setenv FMD_DIR           $::env(CVMFS_DIR)/$::env(SITE)/$::env(OS)/$::env(ARCH)/$NAME/$VERSION
prepend-path LD_LIBRARY_PATH   $::env(FMD_DIR)/lib
prepend-path GCC_INCLUDE_DIR   $::env(FMD_DIR)/include
prepend-path CFLAGS            "-I${FMD_DIR}/include"
prepend-path LDFLAGS           "-L${FMD_DIR}/lib"
MODULE_FILE
) > ${PHYSICAL_MODULES}/${NAME}/${VERSION}
