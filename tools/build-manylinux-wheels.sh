#!/bin/bash
set -e -x

yum install -y gmp
yum install -y mpfr

# Install boost
#
# Install cmake
#
# Install piranha


# Compile wheels
for PYBIN in /opt/python/*/bin; do
    ${PYBIN}/pip install -r numpy 
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
