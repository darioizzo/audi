#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
bash miniconda.sh -b -p $HOME/miniconda
conda config --add channels conda-forge --force

conda_pkgs="cmake>=3.3 clang clangdev eigen obake-devel mppp boost boost-cpp python=3.7"

conda create -q -p $deps_dir -y $conda_pkgs
source activate $deps_dir

export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
export PATH="$deps_dir/bin:$PATH"

export CXX=clang++
export CC=clang

mkdir build

# Install Pybind11 (making sure its the same used in our pipeline)
export PYAUDI_BUILD_DIR=`pwd`
git clone https://github.com/pybind/pybind11.git
cd pybind11
git checkout 4f72ef846fe8453596230ac285eeaa0ce3278bb4
mkdir build
cd build
pwd
cmake \
    -DPYBIND11_TEST=NO \
    -DCMAKE_INSTALL_PREFIX=$PYAUDI_BUILD_DIR \
    -DCMAKE_PREFIX_PATH=$PYAUDI_BUILD_DIR \
    ..
make install
cd ../..

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

