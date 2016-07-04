#!/bin/bash -e
. /etc/profile.d/modules.sh
module load ci
module add ncurses
module add gcc/${GCC_VERSION}
module add openmpi/1.8.8-gcc-${GCC_VERSION}
module add lapack/3.6.0-gcc-${GCC_VERSION}

# Top-level dir.
cd ${WORKSPACE}

#  find the executables.
NBINARIES=$(find -type f -executable -exec file -i '{}' \; | grep 'x-executable; charset=binary' | wc -l)
echo "there are $NBINARIES binaries"
echo $?

## run tests...

# TODO

## tests done

export FMD=${SOFT_DIR}
make install
mkdir -p ${REPO_DIR}
mkdir -p modules
(
cat <<MODULE_FILE
#%Module1.0
## $NAME modulefile
##
proc ModulesHelp { } {
    puts stderr "       This module does nothing but alert the user"
    puts stderr "       that the [module-info name] module is not available"
}

module-whatis   "$NAME $VERSION."

module add ncurses
module add gcc/${GCC_VERSION}
module add openmpi/1.8.8-gcc-${GCC_VERSION}
module add lapack/3.6.0-gcc-${GCC_VERSION}

setenv       FMD_VERSION       $VERSION
setenv       FMD_DIR           /apprepo/$::env(SITE)/$::env(OS)/$::env(ARCH)/$NAME/$VERSION
prepend-path LD_LIBRARY_PATH   $::env(FMD_DIR)/lib
prepend-path GCC_INCLUDE_DIR   $::env(FMD_DIR)/include
prepend-path CFLAGS            "-I${FMD_DIR}/include"
prepend-path LDFLAGS           "-L${FMD_DIR}/lib"
prepend-path PATH              "$::env(FMD_DIR)/bin"
MODULE_FILE
) > modules/$VERSION

mkdir -p ${PHYSICAL_MODULES}/${NAME}
cp modules/$VERSION ${PHYSICAL_MODULES}/${NAME}
module purge

module add ci
echo "checking the modulefile"

module add ${NAME}/${VERSION}

which minenergy

minenergy AV18-UsrgD2000f-v11ls-2.0.int He4.fmd
