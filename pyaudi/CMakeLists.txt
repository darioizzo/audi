YACMA_PYTHON_MODULE(_core SHARED
	core.cpp
	# Python files.
	__init__.py
	test.py
)

set_property(TARGET _core PROPERTY CXX_STANDARD 11)
set_property(TARGET _core PROPERTY CXX_STANDARD_REQUIRED YES)

# Boost.Python version 1.64 (and possibly later?) has issues with default hidden visibility.
if(${Boost_MAJOR_VERSION} GREATER 1 OR (${Boost_MAJOR_VERSION} EQUAL 1 AND ${Boost_MINOR_VERSION} GREATER 63))
string(FIND ${CMAKE_CXX_FLAGS} "-fvisibility=hidden" _PYAUDI_HAS_VHIDDEN)
if (NOT ${_PYAUDI_HAS_VHIDDEN} EQUAL -1)
	message(STATUS "Removing the '-fvisibility=hidden' flag for Boost.Python > 1.63. (release)")
	string(REPLACE "-fvisibility=hidden" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif()
string(FIND ${CMAKE_CXX_FLAGS} "-fvisibility=hidden" _PYAUDI_HAS_VHIDDEN_DEBUG)
if (NOT ${_PYAUDI_HAS_VHIDDEN_DEBUG} EQUAL -1)
	message(STATUS "Removing the '-fvisibility=hidden' flag for Boost.Python > 1.63. (debug)")
	string(REPLACE "-fvisibility=hidden" "" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
	endif()
endif()

# Setup the installation path.
set(PYAUDI_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pyaudi")

# Module needs the Python libs, plus the mandatory libraries.
if(WIN32)
    message(STATUS "Windows detected, linking also to the python library")
    set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
endif()

TARGET_LINK_LIBRARIES(_core ${MANDATORY_LIBRARIES} ${Boost_SERIALIZATION_LIBRARY} ${Boost_PYTHON_LIBRARY})
INSTALL(TARGETS _core
	RUNTIME DESTINATION ${PYAUDI_INSTALL_PATH}
	LIBRARY DESTINATION ${PYAUDI_INSTALL_PATH}
)

# Add the global init file and the additional submodules.
INSTALL(FILES __init__.py test.py DESTINATION ${PYAUDI_INSTALL_PATH})
