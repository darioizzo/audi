#!/bin/bash
set -e -x

# Install a system package required by our library
yum install -y gmp
yum install -y mpfr
yum install -y boost-devel
yum install -y cmake

# Install piranha release 0.8
wget https://github.com/bluescarni/piranha/archive/v0.8.zip
unzip v0.8
cd piranha-0.8
mkdir build
cd build
cmake ../
make install

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
    ${PYBIN}/pip install pyaudi --no-index -f /io/wheelhouse
    (cd $HOME; ${PYBIN}/nosetests pymanylinuxdemo)
done
