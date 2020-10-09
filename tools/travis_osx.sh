#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh;

export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
bash miniconda.sh -b -p $HOME/miniconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda_pkgs="cmake eigen nlopt ipopt boost boost-cpp tbb tbb-devel pybind11"

conda create -q -p $deps_dir -y
source activate $deps_dir
conda install $conda_pkgs -y

export CXX=clang++
export CC=clang

mkdir build

# Install audi
cd build
cmake -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=${AUDI_BUILD_TYPE} \
    -DAUDI_BUILD_AUDI=yes \
    -DAUDI_BUILD_TESTS=yes\
    ..
make -j2 VERBOSE=1
make install
ctest -j4 -V
cd ..

# Install pyaudi
mkdir build_pyaudi
cd build_pyaudi
cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=${AUDI_BUILD_TYPE} \
    -DAUDI_BUILD_AUDI=no \
    -DAUDI_BUILD_PYAUDI=yes \
    -Dpybind11_DIR=$PYAUDI_BUILD_DIR/share/cmake/pybind11/ \
    ..
make -j2 VERBOSE=1
make install
python -c "import pyaudi.test; pyaudi.test.run_test_suite()"

set +e
set +x

