#!/bin/bash
set -e -x

gcc --version


yum -y install epel-release
yum -y install cairo
yum -y install cairo-devel
yum -y install libxml2-devel
yum -y remove cmake28
yum -y install fftw3 fftw3-devel
yum -y install openssl-devel
yum -y install openblas openblas-devel blas blas-devel atlas atlas-devel lapack lapack-devel

# Go to mounted directory
cd io


mkdir build
cd build

mkdir external
cd external

# Install slightly newer version of cmake
curl -O https://cmake.org/files/v3.22/cmake-3.22.1.tar.gz
tar -zxvf cmake-3.22.1.tar.gz
cd cmake-3.22.1
./bootstrap
make -j2
make install
cd ..

# Install eigen3
curl -L -O https://gitlab.com/libeigen/eigen/-/archive/3.2.10/eigen-3.2.10.tar.gz
#mv 3.2.10.tar.gz eigen-3.2.10.tar.gz
#mkdir eigen-3.2.10
#tar --directory=eigen-3.2.10 --strip-components=1 -zxvf eigen-3.2.10.tar.gz
tar -zxvf eigen-3.2.10.tar.gz
cd eigen-3.2.10
mkdir build
cd build
cmake ..
make -j2 && make install
cd ..
cd ..

# Install GROMACS
curl -L -O http://ftp.gromacs.org/pub/gromacs/gromacs-2021.4.tar.gz
tar -zxvf gromacs-2021.4.tar.gz
cd gromacs-2021.4

# Patch CMakeLists.txt
patch="if(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS) \n"
patch+="    set(OLD_CMAKE_FIND_LIBRARY_SUFFIXES \${CMAKE_FIND_LIBRARY_SUFFIXES}) \n"
patch+="    set(CMAKE_FIND_LIBRARY_SUFFIXES \".so\") \n"
patch+="    find_package(ZLIB QUIET) \n"
patch+="    set(CMAKE_FIND_LIBRARY_SUFFIXES \${OLD_CMAKE_FIND_LIBRARY_SUFFIXES}) \n"
patch+="else(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS) \n"
patch+="    find_package(ZLIB QUIET) \n"
patch+="endif(GMX_PREFER_STATIC_LIBS AND BUILD_SHARED_LIBS)"
sed -i "s|find_package(ZLIB QUIET)|$patch|g" CMakeLists.txt

mkdir build
mkdir installed
export GMX_PATH=`pwd`/installed
export GMX_SRC=`pwd`
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$GMX_PATH -DGMX_GPU=off -DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=ON -DCMAKE_CXX_FLAGS="-static-libstdc++" -DGMX_SIMD=SSE2 -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on
make -j2 && make install
cd ..
cd ..

# Out of external
cd ..

# g_mmpbsa compilation
mkdir installed
cmake .. -DGMX_PATH=$GMX_PATH -DGMX_SRC=$GMX_SRC -DBUILD_STATIC=ON -DCMAKE_INSTALL_PREFIX=/io/build/installed
make
make install


