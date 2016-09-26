PROJECT(audi)

ENABLE_TESTING()

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules" "${CMAKE_SOURCE_DIR}/cmake_modules/yacma")

MESSAGE(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Set default build type to "Release".
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Build Option: when active the file main.cpp is built.
OPTION(BUILD_MAIN "Build 'main.cpp'." OFF)

# Build option: enable test set.
OPTION(BUILD_TESTS "Build test set." ON)

# Build option: enable test set.
OPTION(BUILD_PYAUDI "Build the Python bindings." OFF)

# Initialise (empty) list of libraries to link against.
SET(MANDATORY_LIBRARIES "")

# Initial setup of the compilation flags.
include(AudiCompilerLinkerSettings)

# Additional platform-specific setup.
include(AudiPlatformSettings)

# Threading setup.
include(YACMAThreadingSetup)
set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${YACMA_THREADING_CXX_FLAGS}")

# Report the CXX flags.
message(STATUS "Current CXX flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Current CXX debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")

if(BUILD_PYAUDI)
	include(YACMAPythonSetup)
endif()

set(REQUIRED_BOOST_LIBS "")

if(BUILD_PYAUDI OR BUILD_TESTS)
	message(STATUS "Linking boost serialization.")
	list(APPEND REQUIRED_BOOST_LIBS serialization)
endif()

IF(BUILD_TESTS)
	MESSAGE(STATUS "Linking unit tests to Boost.Test.")
	LIST(APPEND REQUIRED_BOOST_LIBS unit_test_framework)
	IF(CMAKE_BUILD_TYPE STREQUAL "Release")
		MESSAGE(STATUS "Linking performance tests to Boost Timer and Boost.Chrono.")
		LIST(APPEND REQUIRED_BOOST_LIBS timer chrono system)
	ENDIF()
ENDIF()

MESSAGE(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")
# Timers are available since 1.48.0.
FIND_PACKAGE(Boost 1.48.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
MESSAGE(STATUS "Detected Boost version: ${Boost_VERSION}")
# Include system Boost headers.
MESSAGE(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})

FIND_PACKAGE(Piranha REQUIRED)
MESSAGE(STATUS "Piranha library found.")
MESSAGE(STATUS "Piranha include dir is: ${Piranha_INCLUDE_DIR}")
INCLUDE_DIRECTORIES(${Piranha_INCLUDE_DIR})

# GMP setup.
FIND_PACKAGE(GMP REQUIRED)
MESSAGE(STATUS "GMP library found.")
MESSAGE(STATUS "GMP include dir is: ${GMP_INCLUDE_DIR}")
MESSAGE(STATUS "GMP library is: ${GMP_LIBRARIES}")
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})

# MPFR setup.
FIND_PACKAGE(MPFR REQUIRED)
MESSAGE(STATUS "MPFR library found.")
MESSAGE(STATUS "MPFR include dir is: ${MPFR_INCLUDE_DIR}")
MESSAGE(STATUS "MPFR library is: ${MPFR_LIBRARIES}")
MESSAGE(STATUS "Checking MPFR version.")
TRY_COMPILE(MPFR_VERSION_CHECK ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake_modules/mpfr_check_version.cpp"
	CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MPFR_INCLUDE_DIR};${GMP_INCLUDE_DIR}")
IF(NOT MPFR_VERSION_CHECK)
	MESSAGE(FATAL_ERROR "Unsupported MPFR version, please upgrade.")
ENDIF(NOT MPFR_VERSION_CHECK)
MESSAGE(STATUS "MPFR version is ok.")
INCLUDE_DIRECTORIES(${MPFR_INCLUDE_DIR})

# NOTE: MPFR should be linked in before GMP, that's why we link GMP here.
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

# Add the subdirectory for the library / headers
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

# Link main to mandatory libraries.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)

IF(BUILD_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
ENDIF(BUILD_TESTS)

if(BUILD_PYAUDI)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/pyaudi")
endif()
