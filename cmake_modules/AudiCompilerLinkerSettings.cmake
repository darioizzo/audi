INCLUDE(CheckCXXCompilerFlag)

message(STATUS "The C++ compiler ID is: ${CMAKE_CXX_COMPILER_ID}")

# Clang detection:
# http://stackoverflow.com/questions/10046114/in-cmake-how-can-i-test-if-the-compiler-is-clang
# http://www.cmake.org/cmake/help/v2.8.10/cmake.html#variable:CMAKE_LANG_COMPILER_ID
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	SET(CMAKE_COMPILER_IS_CLANGXX 1)
ENDIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

macro(AUDI_CHECK_ENABLE_CXX_FLAG flag)
	set(AUDI_CHECK_CXX_FLAG)
	check_cxx_compiler_flag("${flag}" AUDI_CHECK_CXX_FLAG)
	if(AUDI_CHECK_CXX_FLAG)
		message(STATUS "Enabling the '${flag}' compiler flag.")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
	else()
		message(STATUS "Disabling the '${flag}' compiler flag.")
	endif()
	unset(AUDI_CHECK_CXX_FLAG CACHE)
endmacro()

macro(AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG flag)
	set(AUDI_CHECK_DEBUG_CXX_FLAG)
	check_cxx_compiler_flag("${flag}" AUDI_CHECK_DEBUG_CXX_FLAG)
	if(AUDI_CHECK_DEBUG_CXX_FLAG)
		message(STATUS "Enabling the '${flag}' debug compiler flag.")
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${flag}")
	else()
		message(STATUS "Disabling the '${flag}' debug compiler flag.")
	endif()
	unset(AUDI_CHECK_DEBUG_CXX_FLAG CACHE)
endmacro()

# Configuration for GCC.
IF(CMAKE_COMPILER_IS_GNUCXX)
	MESSAGE(STATUS "GNU compiler detected, checking version.")
	# The trouble here is that -g (which implies -g2) results in ICE in some tests and in
	# some pyranha exposition cases. We just append -g1 here, which overrides the default -g.
	message(STATUS "Forcing the debug flag to -g1 for GCC.")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g1")
	# Set the standard flag.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	# Color diagnostic available since GCC 4.9.
	AUDI_CHECK_ENABLE_CXX_FLAG(-fdiagnostics-color=auto)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(CMAKE_COMPILER_IS_CLANGXX)
	MESSAGE(STATUS "Clang compiler detected, checking version.")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	# This used to be necessary with earlier versions of Clang which
	# were not completely compatible with GCC's stdlib. Nowadays it seems
	# unnecessary on most platforms.
	# AUDI_CHECK_ENABLE_CXX_FLAG(-stdlib=libc++)
	# For now it seems like -Wshadow from clang behaves better than GCC's, just enable it here
	# for the time being.
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wshadow)
ENDIF(CMAKE_COMPILER_IS_CLANGXX)

# Common configuration for GCC, Clang and Intel.
if(CMAKE_COMPILER_IS_CLANGXX OR CMAKE_COMPILER_IS_GNUCXX)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wall)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wextra)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wnon-virtual-dtor)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wnoexcept)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wlogical-op)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wconversion)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wdeprecated)
	# In the serialization work, we started hitting the template recursive instantiation
	# limit on clang. This limit is supposed to be at least 1024 in C++11, but for some reason
	# clang sets this to 256, and gcc to 900.
	AUDI_CHECK_ENABLE_CXX_FLAG(-ftemplate-depth=1024)
	# NOTE: this can be useful, but at the moment it triggers lots of warnings in type traits.
	# Keep it in mind for the next time we touch type traits.
	# AUDI_CHECK_ENABLE_CXX_FLAG(-Wold-style-cast)
	# NOTE: disable this for now, as it results in a lot of clutter from Boost.
	# AUDI_CHECK_ENABLE_CXX_FLAG(-Wzero-as-null-pointer-constant)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-pedantic-errors)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-Wdisabled-optimization)
	AUDI_CHECK_ENABLE_CXX_FLAG(-fvisibility-inlines-hidden)
	AUDI_CHECK_ENABLE_CXX_FLAG(-fvisibility=hidden)
	# This is useful when the compiler decides the template backtrace is too verbose.
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-ftemplate-backtrace-limit=0)
	AUDI_CHECK_ENABLE_DEBUG_CXX_FLAG(-fstack-protector-all)
endif()
