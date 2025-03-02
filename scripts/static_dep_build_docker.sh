#!/bin/bash

CWD=$(pwd)

cd /io/external
mkdir apbs_installed
mkdir gromacs_installed

cd apbs/
if [ -d build ]; then
    rm -rf build
fi

mkdir build && cd build

cmake .. \
  -DCMAKE_INSTALL_INCLUDEDIR="include" \
  -DBUILD_DOC=OFF \
  -DAPBS_STATIC_BUILD=ON  \
  -DBUILD_TOOLS=OFF \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX=/io/external/apbs_installed \
  -DENABLE_PYGBE=OFF \
  -DENABLE_BEM=OFF \
  -DENABLE_iAPBS=OFF \
  -DENABLE_GEOFLOW=OFF \
  -DENABLE_OPENMP=ON \
  -DENABLE_PBAM=OFF \
  -DENABLE_PBSAM=OFF \
  -DENABLE_PYTHON=OFF \
  -DENABLE_TESTS=OFF \
  -DFETK_VERSION=57195e55351e04ce6ee0ef56a143c996a9aee7e2 \
  -DGET_NanoShaper=OFF \
  -DCMAKE_C_FLAGS="-fpermissive"

make -j12
make install

#cd /io/external/gromacs
#if [ -d build ]; then
#    rm -rf build
#fi
#mkdir build && cd build

#export GMX_PATH=/io/external/gmx_installed
#export GMX_SRC=/io/external/gromacs

#cmake .. -DGMX_SIMD=NONE -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on \
#             -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_PATH} \
#             -DGMX_PREFER_STATIC_LIBS=ON -DGMX_OPENMP=OFF -DBUILD_SHARED_LIBS=OFF
            
#make -j12
#make install
