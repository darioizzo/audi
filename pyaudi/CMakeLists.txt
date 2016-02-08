include(AudiPythonSetup)

include_directories(${PYTHON_INCLUDE_PATH})

set(PYAUDI_INSTALL_PATH "${PYTHON_MODULES_PATH}/pyaudi")

ADD_LIBRARY(_core SHARED
	core.cpp
	# Python files.
	__init__.py
)

# Initial definition of extra compile flags.
set(PYAUDI_EXTRA_COMPILE_FLAGS "")

# Prefix and suffix setup.
SET_TARGET_PROPERTIES(_core PROPERTIES PREFIX "")
IF(PYDEXTENSION)
	SET_TARGET_PROPERTIES(_core PROPERTIES SUFFIX ".pyd")
ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Apparently on OS X Python expects the .so extension for compiled modules.
	MESSAGE(STATUS "OS X detected, pyaudi's module suffix will be '.so'.")
	SET_TARGET_PROPERTIES(_core PROPERTIES SUFFIX ".so")
ENDIF()

# Strict aliasing setup.
# NOTE: Python < 3 is not strict C, as per:
# http://www.python.org/dev/peps/pep-3123/
# When using GCC or Clang, we add -fno-strict-aliasing here and when compiling core.cpp to prevent possible problems.
# Other compilers should probably use similar flags.
if(${PYTHON_LIBRARY_VERSION} LESS 30 AND (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX))
	message(STATUS "Python < 3 version detected: pyaudi will be compiled with the '-fno-strict-aliasing' flag.")
	set(PYAUDI_EXTRA_COMPILE_FLAGS "${PYAUDI_EXTRA_COMPILE_FLAGS} -fno-strict-aliasing")
	set_target_properties(_core PROPERTIES LINK_FLAGS "-fno-strict-aliasing")
endif()

# Setup the extra compile flags.
message(STATUS "Extra compile flags for pyaudi: ${PYAUDI_EXTRA_COMPILE_FLAGS}")
set_target_properties(_core PROPERTIES COMPILE_FLAGS "${PYAUDI_EXTRA_COMPILE_FLAGS}")

# Module needs the Python libs, plus the mandatory libraries.
TARGET_LINK_LIBRARIES(_core ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS _core
	RUNTIME DESTINATION ${PYAUDI_INSTALL_PATH}
	LIBRARY DESTINATION ${PYAUDI_INSTALL_PATH}
)

# Add the global init file and the additional submodules.
INSTALL(FILES __init__.py DESTINATION ${PYAUDI_INSTALL_PATH})
