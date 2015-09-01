# POSIX thread setup. Intended both for UNIX and Windows (the latter when using some sort of
# pthread emulation/wrapper like pthreads-win32).
IF(CMAKE_USE_PTHREADS_INIT)
	MESSAGE(STATUS "POSIX threads detected.")
	# For POSIX threads, we try to see if the compiler accepts the -pthread flag. It is a bit of a kludge,
	# but I do not have any better idea at the moment. The situation is hairy, e.g.,
	# different systems require different GCC flags:
	# http://gcc.gnu.org/onlinedocs/libstdc++/manual/using_concurrency.html
	CHECK_CXX_COMPILER_FLAG(-pthread AUDI_PTHREAD_COMPILER_FLAG)
	# NOTE: we do not enable the -pthread flag on OS X as it is apparently ignored.
	IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND AUDI_PTHREAD_COMPILER_FLAG)
		MESSAGE(STATUS "Enabling the -pthread compiler flag.")
		# NOTE: according to GCC docs, this sets the flag for both compiler and linker. This should
		# work similarly for clang as well.
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
	ENDIF()
ENDIF()

IF(UNIX)
	# Install path for libraries.
	SET(LIB_INSTALL_PATH "lib")
ENDIF(UNIX)

# OS X setup.
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	IF(CMAKE_COMPILER_IS_CLANGXX)
		# On OS X with clang we need to use libc++.
		MESSAGE(STATUS "Clang compiler on OS X detected, using libc++ as standard library.")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
	ENDIF()
ENDIF()