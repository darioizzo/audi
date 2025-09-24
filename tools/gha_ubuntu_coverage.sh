#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install wget

# Install conda+deps.
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge3.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniforge3/bin:$PATH"
bash miniforge3.sh -b -p $HOME/miniforge3
conda env create --file=audi_devel.yml -q -p $deps_dir
source activate $deps_dir

conda install lcov -y
conda list

export CXXFLAGS="$CXXFLAGS --coverage"

# Create the build dir and cd into it.
mkdir build

# Install audi
cd build
cmake \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=Debug \
    -DAUDI_BUILD_AUDI=yes \
    -DAUDI_BUILD_TESTS=yes \
    ..
make VERBOSE=1 install
ctest -j4 -V

# Create lcov report
lcov --capture --directory . --output-file coverage.info

# Install pyaudi 
cd ..
mkdir build_pyaudi
cd build_pyaudi
cmake \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=Debug \
    -DAUDI_BUILD_AUDI=no \
    -DAUDI_BUILD_PYAUDI=yes \
    ..
make VERBOSE=1 install
python -c "import pyaudi.test; pyaudi.test.run_test_suite()"

# Build the documentation
cd ${GITHUB_WORKSPACE}/doc/doxygen
doxygen Doxyfile
cd ../sphinx
make html linkcheck

set +e
set +x
