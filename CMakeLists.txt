cmake_minimum_required(VERSION 3.2.0)

project(audi VERSION 1.0.3)

enable_testing()

# Here we augment the cmake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/yacma")

# Set default build type to "Release".
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
endif(NOT CMAKE_BUILD_TYPE)

include(YACMACompilerLinkerSettings)

# Setup the CXX flags from YACMA and show final result
YACMA_SETUP_CXX_FLAGS()

# Build Option: when active the file main.cpp is built.
option(BUILD_MAIN "Build 'main.cpp'." OFF)

# Build option: enable test set.
option(BUILD_TESTS "Build test set." ON)

# Build option: enable test set.
option(BUILD_PYAUDI "Build the Python bindings." OFF)

# Initialise (empty) list of libraries to link against.
SET(MANDATORY_LIBRARIES "")

# Threading setup.
include(YACMAThreadingSetup)
set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${YACMA_THREADING_CXX_FLAGS}")

# Report the CXX flags.
message(STATUS "Current CXX flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Current CXX debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")

set(REQUIRED_BOOST_LIBS "")

if(BUILD_PYAUDI)
	include(YACMAPythonSetup)
    message(STATUS "Linking boost python.")
	list(APPEND REQUIRED_BOOST_LIBS python)
endif()

if(BUILD_PYAUDI OR BUILD_TESTS)
	message(STATUS "Linking boost serialization.")
	list(APPEND REQUIRED_BOOST_LIBS serialization)
endif()

if(BUILD_TESTS)
	message(STATUS "Linking unit tests to Boost.Test.")
	list(APPEND REQUIRED_BOOST_LIBS unit_test_framework)
	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		message(STATUS "Linking performance tests to Boost Timer and Boost.Chrono.")
		list(APPEND REQUIRED_BOOST_LIBS timer chrono system)
	endif()
endif()

message(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")
# Timers are available since 1.48.0.
FIND_PACKAGE(Boost 1.48.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
message(STATUS "Detected Boost version: ${Boost_VERSION}")
# Include system Boost headers.
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})

FIND_PACKAGE(Piranha REQUIRED)
message(STATUS "Piranha library found.")
message(STATUS "Piranha include dir is: ${Piranha_INCLUDE_DIR}")
INCLUDE_DIRECTORIES(${Piranha_INCLUDE_DIR})

# GMP setup.
FIND_PACKAGE(GMP REQUIRED)
message(STATUS "GMP library found.")
message(STATUS "GMP include dir is: ${GMP_INCLUDE_DIR}")
message(STATUS "GMP library is: ${GMP_LIBRARIES}")
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})

# MPFR setup.
FIND_PACKAGE(MPFR REQUIRED)
message(STATUS "MPFR library found.")
message(STATUS "MPFR include dir is: ${MPFR_INCLUDE_DIR}")
message(STATUS "MPFR library is: ${MPFR_LIBRARIES}")
message(STATUS "Checking MPFR version.")
TRY_COMPILE(MPFR_VERSION_CHECK ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake/mpfr_check_version.cpp"
	CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MPFR_INCLUDE_DIR};${GMP_INCLUDE_DIR}")
IF(NOT MPFR_VERSION_CHECK)
	message(FATAL_ERROR "Unsupported MPFR version, please upgrade.")
ENDIF(NOT MPFR_VERSION_CHECK)
message(STATUS "MPFR version is ok.")
INCLUDE_DIRECTORIES(${MPFR_INCLUDE_DIR})

# NOTE: MPFR should be linked in before GMP, that's why we link GMP here.
set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

# Add the subdirectory for the library / headers
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

# Link main to mandatory libraries.
IF(BUILD_MAIN)
	add_executable(main main.cpp)
	TARGET_LINK_LIBRARIES(main ${MANDATORY_LIBRARIES})
	set_property(TARGET main PROPERTY CXX_STANDARD 11)
    set_property(TARGET main PROPERTY CXX_STANDARD_REQUIRED YES)
ENDIF(BUILD_MAIN)

IF(BUILD_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
ENDIF(BUILD_TESTS)

if(BUILD_PYAUDI)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/pyaudi")
    if(MINGW)
		message(STATUS "Creating the files for the generation of a binary wheel for MinGW.")
		configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py")
		configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}.txt")
	endif()
endif()
