#!/bin/sh

# Exit on error
set -ev
export CC=mpicc
export CXX=mpicxx
export FC=mpifort

# Print compiler information
echo $PATH
$CC --version
$CXX --version
$FC --version

# log the CMake version (need 3+)
cmake --version

git clone https://github.com/Reference-ScaLAPACK/scalapack.git
cd scalapack

export INSTALL_DIR=$1
export BUILD_DIR=$2
cmake -H. -B$BUILD_DIR -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC
make -C $BUILD_DIR -j2 install
