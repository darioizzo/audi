from setuptools.dist import Distribution
from setuptools import setup

class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True

setup(name="pyaudi",
    version='1.0.3',
    packages=['pyaudi'],
    # Include pre-compiled extension
    package_data={
                'pyaudi': ['_core.so']
                },
    distclass=BinaryDistribution)
