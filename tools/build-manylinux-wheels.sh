#!/bin/bash
set -e -x

yum install -y gmp
yum install -y mpfr

cd /io
mkdir libs
# Install boost
wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.bz2 --no-check-certificate
tar --bzip2 -xf /io/boost_1_62_0.tar.bz2
cd boost_1_62_0
./bootstrap.sh --prefix=/io/libs

#
# Install cmake
#
# Install piranha

# Install and compile pyaudi



# Compile wheels
for PYBIN in /opt/python/*/bin; do
    ${PYBIN}/pip wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair $whl -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    ${PYBIN}/pip install python-manylinux-demo --no-index -f /io/wheelhouse
    (cd $HOME; ${PYBIN}/nosetests pymanylinuxdemo)
done
