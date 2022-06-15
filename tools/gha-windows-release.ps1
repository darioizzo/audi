# Powershell script

# Install conda environment
conda config --set always_yes yes
conda create --name audi cmake eigen obake-devel mppp boost boost-cpp pybind11 python=3.10
conda activate audi

# Define environment variables for clang ...
# ... and make them persistent 
cmd.exe /c "call `"C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars64.bat`" && set > %temp%\vcvars.txt"

Get-Content "$env:temp\vcvars.txt" | Foreach-Object {
  if ($_ -match "^(.*?)=(.*)$") {
    Set-Content "env:\$($matches[1])" $matches[2]
  }
}

#build and run the dcgp tests
mkdir build
cd build
cmake `
    -G "Ninja" `
    -DCMAKE_C_COMPILER=clang-cl `
    -DCMAKE_CXX_COMPILER=clang-cl `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\audi `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\audi `
    -DAUDI_BUILD_AUDI=yes `
    -DAUDI_BUILD_TESTS=yes ..

cmake --build . --target install --config Release
ctest -j4 -V -C Release

cd ..
mkdir build_python
cd build_python

cmake `
    -G "Ninja" `
    -DCMAKE_C_COMPILER=clang-cl `
    -DCMAKE_CXX_COMPILER=clang-cl `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\audi `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\audi `
    -DAUDI_BUILD_AUDI=no `
    -DAUDI_BUILD_PYAUDI=yes ..

cmake --build . --target install --config Release

python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()"
