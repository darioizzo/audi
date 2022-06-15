#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install wget

# Install conda+deps.
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
bash miniconda.sh -b -p $HOME/miniconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -y -q -p $deps_dir c-compiler cxx-compiler cmake eigen obake-devel mppp boost boost-cpp pybind11 python=3.10
source activate $deps_dir

# Create the build dir and cd into it.
mkdir build

# Install audi
cd build
cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=Debug \
    -DAUDI_BUILD_AUDI=yes \
    -DAUDI_BUILD_TESTS=yes \
    ..
make VERBOSE=1 install
ctest -j4 -V
cd ..

# Install pyaudi 
mkdir build_pyaudi
cd build_pyaudi
cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=Debug \
    -DAUDI_BUILD_AUDI=no \
    -DAUDI_BUILD_PYAUDI=yes \
    ..
make VERBOSE=1 install
python -c "import pyaudi.test; pyaudi.test.run_test_suite()"

set +e
set +x
