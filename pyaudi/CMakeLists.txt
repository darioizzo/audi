YACMA_PYTHON_MODULE(_core SHARED
	core.cpp
	# Python files.
	__init__.py
	test.py
)

set_property(TARGET _core PROPERTY CXX_STANDARD 11)
set_property(TARGET _core PROPERTY CXX_STANDARD_REQUIRED YES)

# Setup the installation path.
set(PYAUDI_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pyaudi")

# Module needs the Python libs, plus the mandatory libraries.
TARGET_LINK_LIBRARIES(_core ${MANDATORY_LIBRARIES} ${Boost_SERIALIZATION_LIBRARY} ${Boost_PYTHON_LIBRARY})
INSTALL(TARGETS _core
	RUNTIME DESTINATION ${PYAUDI_INSTALL_PATH}
	LIBRARY DESTINATION ${PYAUDI_INSTALL_PATH}
)

# Add the global init file and the additional submodules.
INSTALL(FILES __init__.py test.py DESTINATION ${PYAUDI_INSTALL_PATH})
