INCLUDE(YACMACompilerLinkerSettings)

# Setup the C++11 flag.
if(YACMA_COMPILER_IS_CLANGXX OR YACMA_COMPILER_IS_INTELXX OR YACMA_COMPILER_IS_GNUCXX)
	message(STATUS "Enabling the '-std=c++11' flag.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
