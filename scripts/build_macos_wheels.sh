#!/bin/bash
set -e -x

CWD=`pwd`

brew install gcc@14 libomp gsl fftw pyenv eigen arpack \
             boost openblas lapack suite-sparse
brew cleanup

eval "$(pyenv init -)"
pyenv install --list
PYVERS=("3.9" "3.10" "3.11" "3.12")
PYTHONS=()
for PYVER in ${PYVERS[@]}
do
    # Install the latest release of the specified Python version using pyenv.
    PYVER="$(pyenv install --list | grep -E "^\\s*$PYVER" | sort -n -t. -k3 | tail -n1)"
    pyenv install $PYVER
    pyenv global $PYVER
    PYTHONS+=($PYVER)
    echo $PYVER
    echo $(python --version)
    python -m pip install --upgrade pip
    python -m pip install -r ${CWD}/dev-requirements.txt
    python -m pip install delocate
done
pyenv global system

CWD=`pwd`

cd external
mkdir apbs_installed
mkdir gmx_installed

export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:/usr/local/opt/lapack/lib/pkgconfig"
export LDFLAGS="-L/usr/local/opt/libomp/lib"
export CPPFLAGS="-I/usr/local/opt/libomp/include"

export CC=gcc-14
export CXX=g++-14

export APBS_INSTALL=${CWD}/external/apbs_installed

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
  -DCMAKE_INSTALL_PREFIX=${APBS_INSTALL} \
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
  -DCMAKE_C_FLAGS="-fpermissive" \
  -DCMAKE_CC_COMIPLER=gcc-14 \
  -DCMAKE_CXX_COMIPLER=g++-14

make 
make install

cd $CWD/external/gromacs
if [ -d build ]; then
    rm -rf build
fi
mkdir build && cd build

export GMX_INSTALL=${CWD}/external/gmx_installed
export GMX_SRC=${CWD}/external/gromacs

cmake .. -DCMAKE_CC_COMIPLER=gcc-14 -DCMAKE_CXX_COMIPLER=g++-14 \
         -DGMX_SIMD=NONE -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on \
         -DGMX_FFT_LIBRARY=fftpack -DCMAKE_INSTALL_PREFIX=${GMX_INSTALL}
            
make
make install

cd ${CWD}

mkdir wheels
mkdir fixed_wheels
for PYTHON in ${PYTHONS[@]}
do
    pyenv global $PYTHON
    echo $(python --version)
    GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} python -m pip install -v --no-deps --no-cache-dir .
    # otool -L build/lib.*/gmx_clusterByFeatures/*.so
    install_name_tool -change @rpath/libgromacs.2.dylib ${GMX_INSTALL}/lib/libgromacs.dylib build/lib.*/g_mmpbsa/*.so
    GMX_INSTALL=${GMX_INSTALL} GMX_SRC=${GMX_SRC} APBS_INSTALL=${APBS_INSTALL} python -m pip wheel -v -w temp_wheels/ --no-deps --no-cache-dir .
    python -m pip uninstall -y g_mmpbsa
    python -m pip install -v --no-deps --no-cache-dir temp_wheels/*.whl
    cp temp_wheels/*.whl wheels/.
    rm -rf build
    rm -rf temp_wheels
done

delocate-listdeps wheels/*.whl
delocate-wheel -w fixed_wheels -v wheels/*.whl
delocate-listdeps fixed_wheels/*.whl
pyenv global system

ls -lrt fixed_wheels/