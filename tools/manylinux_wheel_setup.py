from distutils.core import setup, Extension

extension_module = Extension(
    'dummy',
     sources=['dummy.cpp']
)

setup(name="pyaudi",
    version='1.0.3',
    packages=['pyaudi'],
    ext_modules = [extension_module],
    # Include pre-compiled extension
    package_data={
               	'pyaudi': ['_core.so']
               	},
)
