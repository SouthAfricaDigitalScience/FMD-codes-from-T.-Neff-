#!/bin/bash -e
. /etc/profile.d/modules.sh
module add ci
module add ncurses
module add gcc/${GCC_VERSION}
module add  openmpi/1.8.8-gcc-${GCC_VERSION}
module add lapack/3.6.0-gcc-${GCC_VERSION}

echo "REPO_DIR is "
echo $REPO_DIR
echo "SRC_DIR is "
echo $SRC_DIR
echo "WORKSPACE is "
echo $WORKSPACE
echo "SOFT_DIR is"
echo $SOFT_DIR

mkdir -p ${WORKSPACE}
mkdir -p ${SRC_DIR}
mkdir -p ${SOFT_DIR}

cd ../
make
