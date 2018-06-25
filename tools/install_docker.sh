#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

CMAKE_VERSION="3.11.1"
EIGEN3_VERSION="3.3.4"
BOOST_VERSION="1.67.0"
NLOPT_VERSION="2.4.2"

if [[ ${AUDI_BUILD} == *36 ]]; then
	PYTHON_DIR="cp36-cp36m"
elif [[ ${AUDI_BUILD} == *35 ]]; then
	PYTHON_DIR="cp35-cp35m"
elif [[ ${AUDI_BUILD} == *34 ]]; then
	PYTHON_DIR="cp34-cp34m"
elif [[ ${AUDI_BUILD} == *27 ]]; then
	PYTHON_DIR="cp27-cp27mu"
else
	echo "Invalid build type: ${AUDI_BUILD}"
	exit 1
fi

# HACK: for python 3.x, the include directory
# is called 'python3.xm' rather than just 'python3.x'.
# This confuses the build system of Boost.Python, thus
# we create a symlink to 'python3.x'.
cd /opt/python/${PYTHON_DIR}/include
PY_INCLUDE_DIR_NAME=`ls`
# If the include dir ends with 'm', create a symlink
# without the 'm'.
if [[ $PY_INCLUDE_DIR_NAME == *m ]]; then
	ln -s $PY_INCLUDE_DIR_NAME `echo $PY_INCLUDE_DIR_NAME|sed 's/.$//'`
fi

cd
mkdir install
cd install

# Install Boost
curl -L http://dl.bintray.com/boostorg/release/${BOOST_VERSION}/source/boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2 > boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2
tar xjf boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2
cd boost_`echo ${BOOST_VERSION}|tr "." "_"`
sh bootstrap.sh --with-python=/opt/python/${PYTHON_DIR}/bin/python > /dev/null
./bjam --toolset=gcc link=shared threading=multi cxxflags="-std=c++11" variant=release --with-python --with-serialization --with-iostreams --with-regex --with-chrono --with-timer --with-test --with-system -j2 install > /dev/null
cd ..

# Install gmp (before mpfr as its used by it)
curl -L https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2 > gmp-6.1.2.tar.bz2
tar xvf gmp-6.1.2.tar.bz2  > /dev/null 2>&1
cd gmp-6.1.2 > /dev/null 2>&1
./configure --enable-fat > /dev/null 2>&1
make > /dev/null 2>&1
make install > /dev/null 2>&1
cd ..

# Install mpfr
curl -L http://www.mpfr.org/mpfr-3.1.6/mpfr-3.1.6.tar.gz > mpfr-3.1.6.tar.gz
tar xvf mpfr-3.1.6.tar.gz > /dev/null 2>&1
cd mpfr-3.1.6
./configure > /dev/null 2>&1
make > /dev/null 2>&1
make install > /dev/null 2>&1
cd ..

# Install CMake
curl -L https://github.com/Kitware/CMake/archive/v${CMAKE_VERSION}.tar.gz > v${CMAKE_VERSION}
tar xzf v${CMAKE_VERSION} > /dev/null 2>&1
cd CMake-${CMAKE_VERSION}/
./configure > /dev/null
gmake -j2 > /dev/null
gmake install > /dev/null
cd ..

# Install Eigen
curl -L http://bitbucket.org/eigen/eigen/get/${EIGEN3_VERSION}.tar.gz > ${EIGEN3_VERSION}
tar xzf ${EIGEN3_VERSION} > /dev/null 2>&1
cd eigen*
mkdir build
cd build
cmake ../ > /dev/null
make install > /dev/null
cd ..
cd ..

# Install piranha
curl -L https://github.com/bluescarni/piranha/archive/v0.10.tar.gz > v0.10
tar xvf v0.10 > /dev/null 2>&1
cd piranha-0.10
mkdir build
cd build
cmake ../ > /dev/null
make install > /dev/null 2>&1
cd ..

# Install audi headers
cd /audi
mkdir build_audi
cd build_audi
cmake -DAUDI_BUILD_AUDI=yes -DAUDI_BUILD_TESTS=no -DCMAKE_BUILD_TYPE=Release ../
make install
cd ..

# Compile and install pyaudi (build directory is created by .travis.yml)
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DAUDI_BUILD_AUDI=no -DAUDI_BUILD_PYAUDI=yes -DPYTHON_EXECUTABLE=/opt/python/${PYTHON_DIR}/bin/python ../;
make -j2 install

# Compile wheels
cd wheel
# Copy the installed pyaudi files, wherever they might be in /usr/local,
# into the current dir.
cp -a `find /usr/local/lib -type d -iname 'pyaudi'` ./
# Create the wheel and repair it.
/opt/python/${PYTHON_DIR}/bin/python setup.py bdist_wheel
auditwheel repair dist/pyaudi* -w ./dist2
# Try to install it and run the tests.
cd /
/opt/python/${PYTHON_DIR}/bin/pip install /audi/build/wheel/dist2/pyaudi*
/opt/python/${PYTHON_DIR}/bin/python -c "from pyaudi import test; test.run_test_suite()"

# Upload in PyPi
# This variable will contain something if this is a tagged build (vx.y.z), otherwise it will be empty.
export AUDI_RELEASE_VERSION=`echo "${TRAVIS_TAG}"|grep -E 'v[0-9]+\.[0-9]+.*'|cut -c 2-`
if [[ "${AUDI_RELEASE_VERSION}" != "" ]]; then
    echo "Release build detected, uploading to PyPi."
    /opt/python/${PYTHON_DIR}/bin/pip install twine
    /opt/python/${PYTHON_DIR}/bin/twine upload -u darioizzo /audi/build/wheel/dist2/pyaudi*
fi
