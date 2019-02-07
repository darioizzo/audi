if(UNIX)
# On unix platforms, we set the rpath of the pyaudi libraries
# to the lib subdirectory of the install prefix. This allows
# to install pyaudi in nonstandard locations (e.g., ~/.local)
# without having to muck around with library search paths.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Setup of the header-only pyaudi library.
add_library(pyaudi INTERFACE)
target_link_libraries(pyaudi INTERFACE Audi::audi)
target_include_directories(pyaudi INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/pyaudi> $<INSTALL_INTERFACE:include>)

# Setup of the pyaudi core module.
pybind11_add_module(core core.cpp docstrings.cpp)
target_link_libraries(core PRIVATE pyaudi)
target_compile_options(core PRIVATE "$<$<CONFIG:DEBUG>:${AUDI_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${AUDI_CXX_FLAGS_RELEASE}>")
set_property(TARGET core PROPERTY CXX_STANDARD 11)
set_property(TARGET core PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)

# Setup the installation path.
set(PYAUDI_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pyaudi")
install(TARGETS core
RUNTIME DESTINATION ${PYAUDI_INSTALL_PATH}
LIBRARY DESTINATION ${PYAUDI_INSTALL_PATH}
)

# Add the Python files.
install(FILES __init__.py _version.py test.py DESTINATION ${PYAUDI_INSTALL_PATH})




































