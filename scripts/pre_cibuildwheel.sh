#!/bin/bash
set -e -x

dnf install -y epel-release
dnf install -y https://repo.almalinux.org/almalinux/8/devel/x86_64/os/Packages/suitesparse-static-4.4.6-11.el8.x86_64.rpm

dnf -y install \
    unzip \
    wget \
    arpack-devel \
    arpack-static \
    f2c \
    eigen3-devel \
    boost-devel \
    openblas-serial64 \
    openblas-static \
    openblas-devel \
    lapack-devel \
    lapack-static \
    suitesparse-devel 

CWD=`pwd`

cd external
mkdir apbs_installed
mkdir gmx_installed

cd apbs/
if [ -d build ]; then
    rm -rf build
fi

mkdir build && cd build

cmake .. \
  -DCMAKE_INSTALL_INCLUDEDIR="include" \
  -DBUILD_DOC=OFF \
  -DAPBS_STATIC_BUILD=OFF  \
  -DBUILD_TOOLS=OFF \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=${CWD}/external/apbs_installed \
  -DENABLE_PYGBE=OFF \
  -DENABLE_BEM=OFF \
  -DENABLE_iAPBS=ON \
  -DENABLE_GEOFLOW=OFF \
  -DENABLE_OPENMP=ON \
  -DENABLE_PBAM=OFF \
  -DENABLE_PBSAM=OFF \
  -DENABLE_PYTHON=OFF \
  -DENABLE_TESTS=OFF \
  -DFETK_VERSION=57195e55351e04ce6ee0ef56a143c996a9aee7e2 \
  -DGET_NanoShaper=OFF \
  -DCMAKE_C_FLAGS="-fpermissive"

make 
make install

cd $CWD/external/gromacs
if [ -d build ]; then
    rm -rf build
fi
mkdir build && cd build

export GMX_PATH=${CWD}/external/gmx_installed
export GMX_SRC=${CWD}/external/gromacs

cmake .. -DGMX_SIMD=NONE -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on \
             -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_PATH}
            
make
make install

cd ${CWD}