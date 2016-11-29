#!/usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

export PATH="$deps_dir/bin:$PATH"

if [[ "${BUILD_TYPE}" == "Debug" ]]; then
    cmake -DBUILD_MAIN=no -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=yes ../; 
    make VERBOSE=1;
    ctest -VV
elif [[ "${BUILD_TYPE}" == "Release" ]]; then
    cmake -DBUILD_MAIN=no -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=yes ../;
    make install VERBOSE=1;
    ctest -VV
fi

if [[ "${BUILD_TYPE}" == "Python2" || "${BUILD_TYPE}" == "Python3" ]]; then
    cmake -DPIRANHA_WITH_MSGPACK=yes -DPIRANHA_WITH_BZIP2=yes -DPIRANHA_WITH_ZLIB=yes -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_BUILD_TYPE=Debug -DBUILD_PYRANHA=yes -DCMAKE_CXX_FLAGS_DEBUG=-g0 -DCMAKE_CXX_FLAGS=-Os -DCMAKE_INSTALL_PREFIX=$deps_dir  ../;
    make install VERBOSE=1;
    python -c "import pyranha.test; pyranha.test.run_test_suite()";
fi

set +e
set +x
