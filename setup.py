#!/usr/bin/env python
from distutils.core import setup, Extension

extension_module = Extension(
    'pyaudi/_core',
     sources=['pyaudi/core.cpp'],
     libraries=['boost_python3', 'boost_serialization', 'gmp', 'mpfr'],
     include_dirs=['/io/include']
)

setup(
    name = 'pyaudi',
    version = '1.0',
    description = 'This is a demo package with a compiled C extension.',
    ext_modules = [extension_module],
    packages=['pyaudi']
)
