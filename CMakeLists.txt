# CMake version check.
cmake_minimum_required(VERSION 3.2)

# Module path setup.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules" "${CMAKE_SOURCE_DIR}/cmake_modules/yacma")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Set default build type to "Release".
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

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

# Main audi/pyaudi project version.
set(AUDI_PROJECT_VERSION 1.1)

if(AUDI_BUILD_AUDI)
    # Initial setup of a audi build.
    project(audi VERSION ${AUDI_PROJECT_VERSION})
    enable_testing()

    # Build option: enable test set.
    option(AUDI_BUILD_TESTS "Build test set." ON)

    # Build Option: when active the file main.cpp is built.
    option(AUDI_BUILD_MAIN "Build 'main.cpp'." OFF)

else()
    # Initial setup of a pyaudi build.
    project(pyaudi VERSION ${AUDI_PROJECT_VERSION})
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
  # Disable the idiotic minmax macros on MSVC, some annoying warnings,
  # and enable the bigobj option.
  list(APPEND AUDI_CXX_FLAGS_DEBUG "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj")
  list(APPEND AUDI_CXX_FLAGS_RELEASE "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj")
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

# Some flags that generate warnings due to Eigen or Piranha obsolete versions
list(REMOVE_ITEM AUDI_CXX_FLAGS_DEBUG "-Wduplicated-branches")
list(REMOVE_ITEM AUDI_CXX_FLAGS_DEBUG "-Wold-style-cast")

if(AUDI_BUILD_PYAUDI)
    # audi dependencies.
    include(YACMAPythonSetup)
    # Python version check.
    if(${PYTHON_VERSION_MAJOR} LESS 2 OR (${PYTHON_VERSION_MAJOR} EQUAL 2 AND ${PYTHON_VERSION_MINOR} LESS 7))
        message(FATAL_ERROR "Minimum supported Python version is 2.7.")
    endif()

    # NOTE: for the time being, require that audi/pyaudi versions are matching exactly.
    find_package(audi ${AUDI_PROJECT_VERSION} EXACT REQUIRED)

    # Internal variable that will be used to tell AudiFindBoost to locate Boost.Python.
    set(_AUDI_FIND_BOOST_PYTHON TRUE)
endif()

# Boost setup (common to audi/pyaudi).
include(AudiFindBoost)

if(AUDI_BUILD_AUDI)
    # audi dependencies.

    # Eigen setup
    find_package(Eigen3 REQUIRED)
    message(STATUS "Eigen include directory: ${EIGEN3_INCLUDE_DIR}")
    message(STATUS "Eigen version detected: ${EIGEN3_VERSION}")
    include_directories(${EIGEN3_INCLUDE_DIR})
    
    # Piranha setup.
    find_package(Piranha REQUIRED)
    message(STATUS "Piranha library found.")
    message(STATUS "Piranha include dir is: ${Piranha_INCLUDE_DIR}")
    include_directories(${Piranha_INCLUDE_DIR})

    # GMP setup.
    find_package(GMP REQUIRED)
    message(STATUS "GMP library found.")
    message(STATUS "GMP include dir is: ${GMP_INCLUDE_DIR}")
    message(STATUS "GMP library is: ${GMP_LIBRARIES}")
    include_directories(${GMP_INCLUDE_DIR})
    
    # MPFR setup.
    FIND_PACKAGE(MPFR REQUIRED)
    message(STATUS "MPFR library found.")
    message(STATUS "MPFR include dir is: ${MPFR_INCLUDE_DIR}")
    message(STATUS "MPFR library is: ${MPFR_LIBRARIES}")
    message(STATUS "Checking MPFR version.")
    try_compile(MPFR_VERSION_CHECK ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake_modules/mpfr_check_version.cpp"
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MPFR_INCLUDE_DIR};${GMP_INCLUDE_DIR}")
    if(NOT MPFR_VERSION_CHECK)
        message(FATAL_ERROR "Unsupported MPFR version, please upgrade.")
    endif()
    message(STATUS "MPFR version is ok.")
    include_directories(${MPFR_INCLUDE_DIR})

    # Setup of the header-only audi library.
    add_library(audi INTERFACE)
    target_link_libraries(audi INTERFACE Threads::Threads Boost::boost Boost::serialization)    
    target_link_libraries(audi INTERFACE Eigen3::eigen3)
    
    # This sets up the include directory to be different if we build
    target_include_directories(audi INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

    # These lines are temporary with piranha v0.10, will be simplified when this is updated
    target_include_directories(audi INTERFACE ${Piranha_INCLUDE_DIR})
    target_include_directories(audi INTERFACE ${MPFR_INCLUDE_DIR})
    target_include_directories(audi INTERFACE ${GMP_INCLUDE_DIR})
    # First mpfr then gmp otherwise mingw is unhappy
    target_link_libraries(audi INTERFACE ${MPFR_LIBRARIES})
    target_link_libraries(audi INTERFACE ${GMP_LIBRARIES})
    

    # Builds the main file
    IF(AUDI_BUILD_MAIN)
        add_executable(main main.cpp)
        target_link_libraries(main audi)
        set_property(TARGET main PROPERTY CXX_STANDARD 11)
        set_property(TARGET main PROPERTY CXX_STANDARD_REQUIRED YES)
        set_property(TARGET main PROPERTY CXX_EXTENSIONS NO)
    ENDIF()

    # Builds the tests (performance in Release, the others in Debug)
    if(AUDI_BUILD_TESTS)
        add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
    endif()

    # Setup of the export.
    install(TARGETS audi EXPORT audi_export)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/audi-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/audi-config.cmake" @ONLY)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/audi-config.cmake" DESTINATION "lib/cmake/audi")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindGMP.cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindMPFR.cmake"
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindPiranha.cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/AudiFindBoost.cmake"
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
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt" @ONLY)
    endif()
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        # NOTE: this is necessary on linux but harmful on mingw.
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.cfg" @ONLY)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py" @ONLY)
    endif()
endif()



