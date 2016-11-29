#!/usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

export PATH="$deps_dir/bin:$PATH"

# This variable will contain something if this build is tagged, otherwise it will be empty.
export AUDI_RELEASE_VERSION=`echo "${TRAVIS_TAG}"|grep -E 'v[0-9]+\.[0-9]+.*'|cut -c 2-`


if [[ "${BUILD_TYPE}" == "Debug" ]]; then
    cmake -DBUILD_MAIN=no -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=yes ../;
    make VERBOSE=1;
    ctest -VV
elif [[ "${BUILD_TYPE}" == "Release" ]]; then
    cmake -DBUILD_MAIN=no -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=yes ../;
    make install VERBOSE=1;
    ctest -VV
fi

if [[ "${BUILD_TYPE}" == "Python27" ]]; then
    export PYAUDI_SYSTEM_DIRECTORY=/home/travis/local/lib/python2.7/site-packages/pyaudi
elif [[ "${BUILD_TYPE}" == "Python34" ]]; then
    export PYAUDI_SYSTEM_DIRECTORY=/home/travis/local/lib/python3.4/site-packages/pyaudi
elif [[ "${BUILD_TYPE}" == "Python35" ]]; then
    export PYAUDI_SYSTEM_DIRECTORY=/home/travis/local/lib/python3.5/site-packages/pyaudi
fi

if [[ "${BUILD_TYPE}" == "Python27" || "${BUILD_TYPE}" == "Python34" || "${BUILD_TYPE}" == "Python35" ]]; then
    cmake -DBUILD_MAIN=no -DCMAKE_PREFIX_PATH=$deps_dir -DCMAKE_INSTALL_PREFIX=$deps_dir -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=no -DBUILD_PYAUDI=yes ../;
    make install VERBOSE=1;
    python -c "import pyaudi.test; pyaudi.test.run_test_suite()";
fi



set +e
set +x
