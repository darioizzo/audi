mkdir build
cd build
cmake `
    -G "Visual Studio 17 2022" `
    -A x64 `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\audi_devel\Library `
    -DAUDI_BUILD_AUDI=yes `
    -DAUDI_BUILD_TESTS=yes `
    ..

cmake --build . --config Release --target install -- /m 
ctest -j4 -V -C Release --output-on-failure

cd ..
mkdir build_python
cd build_python

cmake `
    -G "Visual Studio 17 2022" `
    -A x64 `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\audi_devel\Library `
    -DAUDI_BUILD_AUDI=no `
    -DAUDI_BUILD_PYAUDI=yes `
    ..

cmake --build . --config Release --target install -- /m 

python -c "from pyaudi import test; test.run_test_suite();"
