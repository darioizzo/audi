import os
import re
import sys
import shutil

def wget(url, out):
    import urllib.request
    print('Downloading "' + url + '" as "' + out + '"')
    urllib.request.urlretrieve(url, out)


def rm_fr(path):
    import os
    import shutil
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def run_command(raw_command, directory=None, verbose=True):
    # Helper function to run a command and display optionally its output
    # unbuffered.
    import shlex
    import sys
    from subprocess import Popen, PIPE, STDOUT
    print(raw_command)
    proc = Popen(shlex.split(raw_command), cwd=directory,
                 stdout=PIPE, stderr=STDOUT)
    if verbose:
        output = ''
        while True:
            line = proc.stdout.readline()
            if not line:
                break
            line = str(line, 'utf-8')
            # Don't print the newline character.
            print(line[:-1])
            sys.stdout.flush()
            output += line
        proc.communicate()
    else:
        output = str(proc.communicate()[0], 'utf-8')
    if proc.returncode:
        raise RuntimeError(output)
    return output


# Build type setup.
BUILD_TYPE = os.environ['BUILD_TYPE']
is_release_build = (os.environ['APPVEYOR_REPO_TAG'] == 'true') and bool(
    re.match(r'v[0-9]+\.[0-9]+.*', os.environ['APPVEYOR_REPO_TAG_NAME']))
if is_release_build:
    print("Release build detected, tag is '" +
          os.environ['APPVEYOR_REPO_TAG_NAME'] + "'")
is_python_build = 'Python' in BUILD_TYPE


# Get mingw and set the path.
wget(r'https://github.com/bluescarni/binary_deps/raw/master/x86_64-6.2.0-release-posix-seh-rt_v5-rev1.7z', 'mw64.7z')
run_command(r'7z x -oC:\\ mw64.7z', verbose=False)
ORIGINAL_PATH = os.environ['PATH']
os.environ['PATH'] = r'C:\\mingw64\\bin;' + os.environ['PATH']

# Download common deps.
wget(r'https://github.com/bluescarni/binary_deps/raw/master/gmp_mingw_64.7z', 'gmp.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/mpfr_mingw_64.7z', 'mpfr.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/boost_mingw_64.7z', 'boost.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/eigen3.7z', 'eigen3.7z')

# Extract them.
run_command(r'7z x -aoa -oC:\\ gmp.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ mpfr.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ boost.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ eigen3.7z', verbose=False)


# Download piranha 0.10 https://github.com/bluescarni/piranha/archive/v0.10.zip
wget(r'https://github.com/bluescarni/piranha/archive/v0.10.zip', 'piranhav10.zip')
run_command(r'unzip piranhav10.zip', verbose=False)
# Move to the directory created and make piranha install its headers
os.chdir('piranha-0.10')
os.makedirs('build')
os.chdir('build')
print("Installing piranha")
run_command(
    r'cmake -G "MinGW Makefiles" .. -DCMAKE_INSTALL_PREFIX=c:\\local ', verbose=False)
run_command(r'mingw32-make install VERBOSE=1', verbose=False)
os.chdir('../../')
print("Piranha sucessfully installed .. continuing")

# Download mppp 0.11 https://github.com/bluescarni/mppp/archive/v0.11.zip
wget(r'https://github.com/bluescarni/mppp/archive/v0.11.zip', 'mpppv10.zip')
run_command(r'unzip mpppv10.zip', verbose=False)
# Move to the directory created and make piranha install its headers
os.chdir('mppp-0.11')
os.makedirs('build')
os.chdir('build')
print("Installing mppp")
run_command(
    r'cmake -G "MinGW Makefiles" .. -DMPPP_WITH_QUADMATH=yes -DCMAKE_INSTALL_PREFIX=c:\\local ', verbose=False)
run_command(r'mingw32-make install VERBOSE=1', verbose=False)
os.chdir('../../')
print("mppp sucessfully installed .. continuing")

# Setup of the dependencies for a Python build.
if is_python_build:
    if BUILD_TYPE == 'Python35':
        python_version = '35'
    elif BUILD_TYPE == 'Python36':
        python_version = '36'
    elif BUILD_TYPE == 'Python27':
        python_version = '27'
    else:
        raise RuntimeError('Unsupported Python build: ' + BUILD_TYPE)
    python_package = r'python' + python_version + r'_mingw_64.7z'
    boost_python_package = r'boost_python_' + python_version + r'_mingw_64.7z'
    # Remove any existing Python installation.
    rm_fr(r'c:\\Python' + python_version)
    # Set paths.
    pinterp = r'c:\\Python' + python_version + r'\\python.exe'
    pip = r'c:\\Python' + python_version + r'\\scripts\\pip'
    twine = r'c:\\Python' + python_version + r'\\scripts\\twine'
    pyaudi_install_path = r'C:\\Python' + \
        python_version + r'\\Lib\\site-packages\\pyaudi'

    # Get Python.
    wget(r'https://github.com/bluescarni/binary_deps/raw/master/' +
         python_package, 'python.7z')
    run_command(r'7z x -aoa -oC:\\ python.7z', verbose=False)
        # Install pip and deps.
    wget(r'https://bootstrap.pypa.io/get-pip.py', 'get-pip.py')
    run_command(pinterp + ' get-pip.py')
    if is_release_build:
        # call pip via python, workaround to avoid path issues when calling pip from win
        # (https://github.com/pypa/pip/issues/1997)
        run_command(pinterp + r' -m pip install twine')

    # Download the python lib https://github.com/mitsuba-renderer/dependencies_win64/raw/master/lib/python36.lib
    wget(r'https://github.com/mitsuba-renderer/dependencies_win64/raw/master/lib/python'+python_version+r'.lib', r'python'+python_version+r'.lib')
    shutil.move(r'python' + python_version + r'.lib', r'C:\\Python' + python_version + r'\\libs\\')
    # Download pybind11 https://github.com/pybind/pybind11/archive/v2.2.4.zip
    wget(r'https://github.com/pybind/pybind11/archive/v2.2.4.zip', 'pybind11_v224.zip')
    run_command(r'unzip pybind11_v224.zip', verbose=False)
    # Move to the directory created and make piranha install its headers
    os.chdir('pybind11-2.2.4')
    os.makedirs('build')
    os.chdir('build')
    print("Installing pybind11")
    run_command(
       r'cmake -G "MinGW Makefiles" ..' + ' ' +
       r'-DPYBIND11_TEST=OFF' + ' ' +
       r'-DPYTHON_PREFIX=C:\\Python' + python_version + ' ' +
       r'-DPYTHON_EXECUTABLE=C:\\Python' + python_version + r'\\python.exe', verbose=True)
    run_command(r'mingw32-make install VERBOSE=1', verbose=False)
    os.chdir('../../')
    print("pybind11 sucessfully installed .. continuing")
    shutil.move(r'C:\\Python' + python_version + r'\\libs\\python' + python_version + r'.lib', r'C:\\Python' + python_version + r'\\libs\\python' + python_version + r'OLD.lib')

# Set the path so that the precompiled libs can be found.
os.environ['PATH'] = os.environ['PATH'] + r';c:\\local\\lib'

# Proceed to the build.
common_cmake_opts = r'-DCMAKE_PREFIX_PATH=c:\\local -DCMAKE_INSTALL_PREFIX=c:\\local -DAUDI_WITH_MPPP=yes'

# Configuration step.
if is_python_build:
    os.makedirs('build_audi')
    os.chdir('build_audi')
    run_command(
        r'cmake -G "MinGW Makefiles" ..  -DCMAKE_BUILD_TYPE=Release -DAUDI_BUILD_TESTS=no -DAUDI_BUILD_AUDI=yes -DAUDI_BUILD_PYAUDI=no' + ' ' + common_cmake_opts)
    run_command(r'mingw32-make install VERBOSE=1 -j2')
    os.chdir('..')
    os.makedirs('build_pyaudi')
    os.chdir('build_pyaudi')
    run_command(r'cmake -G "MinGW Makefiles" ..  -DPYAUDI_INSTALL_PATH=c:\\local -DAUDI_BUILD_AUDI=no -DAUDI_BUILD_PYAUDI=yes -DCMAKE_BUILD_TYPE=Release ' + common_cmake_opts + ' ' +
                r'-DPYTHON_EXECUTABLE=C:\\Python' + python_version + r'\\python.exe -DPYTHON_LIBRARY=C:\\Python' + python_version + r'\\libs\\python' + python_version + r'.dll')
    run_command(r'mingw32-make install VERBOSE=1 -j2')
elif BUILD_TYPE in ['Release', 'Debug']:
    os.makedirs('build_audi')
    os.chdir('build_audi')
    cmake_opts = r'-DCMAKE_BUILD_TYPE=' + BUILD_TYPE + \
        r' -DAUDI_BUILD_TESTS=yes ' + common_cmake_opts
    run_command(r'cmake -G "MinGW Makefiles" .. ' + cmake_opts)
    run_command(r'mingw32-make install VERBOSE=1 -j2')
    run_command(r'ctest')
else:
    raise RuntimeError('Unsupported build type: ' + BUILD_TYPE)


# Packaging.
if is_python_build:
    # Run the Python tests.
    run_command(
        pinterp + r' -c "from pyaudi import test; test.run_test_suite()"')
    # Build the wheel.
    os.chdir('wheel')
    shutil.move(pyaudi_install_path, r'.')
    wheel_libs = 'mingw_wheel_libs_python{}.txt'.format(python_version[0])
    DLL_LIST = [_[:-1] for _ in open(wheel_libs, 'r').readlines()]
    for _ in DLL_LIST:
        shutil.copy(_, 'pyaudi')
    run_command(pinterp + r' setup.py bdist_wheel')
    os.environ['PATH'] = ORIGINAL_PATH
    # workaround necessary to be able to call pip install via python (next line) and
    # not find a pyaudi already in the pythonpath
    os.makedirs('garbage')
    shutil.move('pyaudi', r'garbage')
    shutil.move('pyaudi.egg-info', r'garbage')
    # call pip via python, workaround to avoid path issues when calling pip from win
    # (https://github.com/pypa/pip/issues/1997)
    run_command(pinterp + r' -m pip install dist\\' + os.listdir('dist')[0])
    run_command(
        pinterp + r' -c "from pyaudi import test; test.run_test_suite()"', directory=r'c:\\')
    if is_release_build:
        run_command(twine + r' upload -u darioizzo dist\\' +
                    os.listdir('dist')[0])
