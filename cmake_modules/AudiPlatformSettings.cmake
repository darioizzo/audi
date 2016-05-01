if(UNIX)
	# Install path for libraries.
	set(LIB_INSTALL_PATH "lib")
endif()

# OS X setup.
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(CMAKE_MACOSX_RPATH OFF)
	if(CMAKE_COMPILER_IS_CLANGXX)
		# On OS X with clang we need to use libc++.
		message(STATUS "Clang compiler on OS X detected, using libc++ as standard library.")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
	endif()
endif()
