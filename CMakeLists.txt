# CMake version check.
cmake_minimum_required(VERSION 3.3)

# Set default build type to "Release".
# NOTE: this should be done before the project command since the latter can set
# CMAKE_BUILD_TYPE itself (it does so for nmake).
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
	    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
endif()

# Main audi/pyaudi project version.
set(AUDI_PROJECT_VERSION 1.7)

# Module path setup.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules" "${CMAKE_SOURCE_DIR}/cmake_modules/yacma")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Main build options: build audi or pyaudi. They cannot be on at the same time,
# and only one must be chosen.
option(AUDI_BUILD_AUDI "Build audi." ON)
option(AUDI_BUILD_PYAUDI "Build pyaudi." OFF)

# Check consistency.
if(AUDI_BUILD_AUDI AND AUDI_BUILD_PYAUDI)
    message(FATAL_ERROR "Please select whether to build audi or pyaudi: you cannot build them both at the same time.")
endif()

if((NOT AUDI_BUILD_AUDI) AND (NOT AUDI_BUILD_PYAUDI))
    message(FATAL_ERROR "Please select if you want to build audi or pyudi.")
endif()

if(AUDI_BUILD_AUDI)
    # Initial setup of a audi build.
    project(audi VERSION ${AUDI_PROJECT_VERSION} LANGUAGES CXX C)
    enable_testing()

    # Build option: enable test set.
    option(AUDI_BUILD_TESTS "Build test set." ON)

    # Build Option: when active the file main.cpp is built.
    option(AUDI_BUILD_MAIN "Build 'main.cpp'." OFF)

else()
    # Initial setup of a pyaudi build.
    project(pyaudi VERSION ${AUDI_PROJECT_VERSION} LANGUAGES CXX C)
endif()

# Common general bits.

# Initial setup of compiler flags.
include(YACMACompilerLinkerSettings)

# Threading setup.
include(YACMAThreadingSetup)

# Assemble the flags.
set(AUDI_CXX_FLAGS_DEBUG ${YACMA_CXX_FLAGS} ${YACMA_CXX_FLAGS_DEBUG} ${YACMA_THREADING_CXX_FLAGS})
set(AUDI_CXX_FLAGS_RELEASE ${YACMA_CXX_FLAGS} ${YACMA_THREADING_CXX_FLAGS})
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND YACMA_COMPILER_IS_CLANGXX)
  message(STATUS "Clang compiler on OSX detected, setting the standard library to 'libc++'.")
  list(APPEND AUDI_CXX_FLAGS_DEBUG "-stdlib=libc++")
  list(APPEND AUDI_CXX_FLAGS_RELEASE "-stdlib=libc++")
endif()
if(YACMA_COMPILER_IS_MSVC)
  include(CheckCXXCompilerFlag)
  # Disable the idiotic minmax macros on MSVC, some annoying warnings,
  # and enable the bigobj option.
  # Also, enable the WIN32_LEAN_AND_MEAN definition.
  list(APPEND AUDI_CXX_FLAGS_DEBUG "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj" "-DWIN32_LEAN_AND_MEAN")
  list(APPEND AUDI_CXX_FLAGS_RELEASE "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj" "-DWIN32_LEAN_AND_MEAN")
  # Enable strict conformance mode, if supported.
  set(CMAKE_REQUIRED_QUIET TRUE)
  check_cxx_compiler_flag("/permissive-" _AUDI_MSVC_SUPPORTS_STRICT_CONFORMANCE)
  unset(CMAKE_REQUIRED_QUIET)
  if(_AUDI_MSVC_SUPPORTS_STRICT_CONFORMANCE)
    message(STATUS "The '/permissive-' flag is supported, enabling it.")
    list(APPEND AUDI_CXX_FLAGS_DEBUG "/permissive-")
    list(APPEND AUDI_CXX_FLAGS_RELEASE "/permissive-")
  endif()
  unset(_AUDI_MSVC_SUPPORTS_STRICT_CONFORMANCE)
  if(YACMA_COMPILER_IS_CLANGXX)
    # clang-cl emits various warnings from GMP/MPFR, let's just silence them.
    # NOTE: at one point in the recent past, MSVC added an options similar to GCC's isystem:
    # https://blogs.msdn.microsoft.com/vcblog/2017/12/13/broken-warnings-theory/
    # We probably just need to wait for this to be picked up by CMake/clang-cl. Let's
    # revisit the issue in the future.
    list(APPEND _AUDI_CLANG_CL_DISABLED_WARNINGS
      "-Wno-unused-variable"
      "-Wno-inconsistent-dllimport"
      "-Wno-unknown-pragmas"
      "-Wno-unused-parameter"
      "-Wno-sign-compare"
      "-Wno-deprecated-declarations"
      "-Wno-deprecated-dynamic-exception-spec"
      "-Wno-old-style-cast"
      "-Wno-sign-conversion"
      "-Wno-non-virtual-dtor"
      "-Wno-deprecated"
      "-Wno-shadow"
      "-Wno-shorten-64-to-32"
      "-Wno-reserved-id-macro"
      "-Wno-undef"
      "-Wno-c++98-compat-pedantic"
      "-Wno-documentation-unknown-command"
      "-Wno-zero-as-null-pointer-constant"
      "-Wno-language-extension-token"
      "-Wno-gnu-anonymous-struct"
      "-Wno-nested-anon-types"
      "-Wno-documentation"
      "-Wno-comma"
      "-Wno-nonportable-system-include-path"
      "-Wno-global-constructors"
      "-Wno-redundant-parens"
      "-Wno-exit-time-destructors"
      "-Wno-missing-noreturn"
      "-Wno-switch-enum"
      "-Wno-covered-switch-default"
      "-Wno-float-equal"
      "-Wno-double-promotion"
      "-Wno-microsoft-enum-value"
      "-Wno-missing-prototypes"
      "-Wno-implicit-fallthrough"
      "-Wno-format-nonliteral"
      "-Wno-cast-qual"
      "-Wno-disabled-macro-expansion"
      "-Wno-unused-private-field"
      "-Wno-unused-template"
      "-Wno-unused-macros"
      "-Wno-extra-semi-stmt"
      "-Wno-c++98-compat")
    list(APPEND AUDI_CXX_FLAGS_DEBUG ${_AUDI_CLANG_CL_DISABLED_WARNINGS})
    list(APPEND AUDI_CXX_FLAGS_RELEASE ${_AUDI_CLANG_CL_DISABLED_WARNINGS})
    unset(_AUDI_CLANG_CL_DISABLED_WARNINGS)
  endif()
endif()
if(YACMA_COMPILER_IS_INTELXX)
  # NOTE: on MSVC we use the push/pop pragmas, but they do not seem to work on Intel (the pragmas
  # in icc influence the behaviour at instantiation point, not at definition point).
  # These warnings are useful in principle, but they are generated a lot from cereal and we have no
  # way of disabling them selectively. Just rely on the other compilers to provde good diagnostic.
  list(APPEND AUDI_CXX_FLAGS_DEBUG "-diag-disable" "2259,1682,68")
  list(APPEND AUDI_CXX_FLAGS_RELEASE "-diag-disable" "2259,1682,68")
endif()
if(MINGW)
	# Flag needed to deal with big binaries in MinGW.
	message(STATUS "Enabling the '-Wa,-mbig-obj' flag in MinGW builds.")
	list(APPEND AUDI_CXX_FLAGS_DEBUG "-Wa,-mbig-obj")
	list(APPEND AUDI_CXX_FLAGS_RELEASE "-Wa,-mbig-obj")
endif()

# Some flags that generate warnings due to Eigen obsolete versions
list(REMOVE_ITEM AUDI_CXX_FLAGS_DEBUG "-Wduplicated-branches")
list(REMOVE_ITEM AUDI_CXX_FLAGS_DEBUG "-Wold-style-cast")

# Configure config.hpp.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/include/audi/config.hpp" @ONLY)

if(AUDI_BUILD_PYAUDI)
    # audi dependencies.
    include(YACMAPythonSetup)
    # Python version check.
    if(${PYTHON_VERSION_MAJOR} LESS 3)
        message(FATAL_ERROR "Minimum supported Python version is 3.0")
    endif()
endif()

# Boost setup (common to audi/pyaudi).
if(AUDI_BUILD_TESTS)
    # Internal variable that will be used to tell PagmoFindBoost to locate the
    # Boost unit test framework, if tests are required.
    set(_AUDI_FIND_BOOST_UNIT_TEST_FRAMEWORK TRUE)
endif()
include(AudiFindBoost)

if(AUDI_BUILD_AUDI)
    # audi dependencies.

    # Eigen setup
    find_package(Eigen3 REQUIRED)
    message(STATUS "Eigen include directory: ${EIGEN3_INCLUDE_DIR}")
    message(STATUS "Eigen version detected: ${EIGEN3_VERSION}")
    
    # obake setup.
    find_package(obake REQUIRED)
    message(STATUS "obake library found.")

    # mp++ setup.
    find_package(mp++ REQUIRED)
    message(STATUS "mp++ library found.")

    # Setup of the header-only audi library.
    add_library(audi INTERFACE)
    target_link_libraries(audi INTERFACE Threads::Threads Boost::boost Boost::serialization obake::obake Eigen3::eigen3 mp++::mp++)
    
    # This sets up the include directory to be different if we build
    target_include_directories(audi INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)


    # Builds the main file
    IF(AUDI_BUILD_MAIN)
        add_executable(main main.cpp)
        target_link_libraries(main audi)
        set_property(TARGET main PROPERTY CXX_STANDARD 17)
        set_property(TARGET main PROPERTY CXX_STANDARD_REQUIRED YES)
        set_property(TARGET main PROPERTY CXX_EXTENSIONS NO)
    ENDIF()

    # Builds the tests (performance in Release, the others in Debug)
    if(AUDI_BUILD_TESTS)
        # Internal variable that will be used to tell PagmoFindBoost to locate the
        # Boost unit test framework, if tests are required.
        set(_AUDI_FIND_BOOST_UNIT_TEST_FRAMEWORK TRUE)
        add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
    endif()

    # Configure the doc files.
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py" @ONLY)

    # Setup of the export.
    install(TARGETS audi EXPORT audi_export)
    # Setup of the optional deps.
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/audi-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/audi-config.cmake" @ONLY)

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/audi-config.cmake" DESTINATION "lib/cmake/audi")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/AudiFindBoost.cmake"
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindEigen3.cmake" DESTINATION "lib/cmake/audi")
    install(EXPORT audi_export NAMESPACE Audi:: DESTINATION lib/cmake/audi)
    # Take care of versioning.
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/audi-config-version.cmake" VERSION ${audi_VERSION}
            COMPATIBILITY ExactVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/audi-config-version.cmake" DESTINATION "lib/cmake/audi")

    # Do the actual library installation.
    install(DIRECTORY include/ DESTINATION include)
endif()

if(AUDI_BUILD_PYAUDI)
    # Add the pyaudi subdirectory.
    add_subdirectory("${CMAKE_SOURCE_DIR}/pyaudi")
    if(MINGW)
        message(STATUS "Creating the files for the generation of a binary wheel.")
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py" @ONLY)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" @ONLY)
    endif()
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        # NOTE: this is necessary on linux but harmful on mingw.
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.cfg" @ONLY)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py" @ONLY)
    endif()
endif()
