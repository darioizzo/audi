# Python version check.
if(${PYTHON_VERSION_MAJOR} LESS 3)
    message(FATAL_ERROR "Minimum supported Python version is 3.0")
endif()

# Setting the rpath
if(UNIX)
    # On unix platforms, we set the rpath of the pygmo libraries
    # to the lib subdirectory of the install prefix. This allows
    # to install pygmo in nonstandard locations (e.g., ~/.local)
    # without having to muck around with library search paths.
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Create the module
YACMA_PYTHON_MODULE(core core.cpp docstrings.cpp)

# Finding audi 
#NOTE: for the time being, require that audi/pyaudi versions are matching exactly.
find_package(audi ${AUDI_PROJECT_VERSION} EXACT REQUIRED)

# Finding pybind11
find_package(pybind11 REQUIRED)

# Setup of the header-only pyaudi library.
add_library(pyaudi INTERFACE)
target_link_libraries(pyaudi INTERFACE ${PYAUDI_BP_TARGET} Boost::disable_autolinking Audi::audi )
target_include_directories(pyaudi INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/pyaudi> $<INSTALL_INTERFACE:include>)

# Setup of the pyaudi core module.
target_link_libraries(core PRIVATE ${PYAUDI_BP_TARGET} pyaudi)
target_compile_options(core PRIVATE "$<$<CONFIG:DEBUG>:${AUDI_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${AUDI_CXX_FLAGS_RELEASE}>")
target_include_directories(core SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
target_compile_definitions(core PRIVATE "${pybind11_DEFINITIONS}")

set_property(TARGET core PROPERTY CXX_STANDARD 17)
set_property(TARGET core PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)

# Helper file with version number.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_version.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/_version.py" @ONLY)   

# Setup the installation path.
set(PYAUDI_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/pyaudi")
install(TARGETS core
RUNTIME DESTINATION ${PYAUDI_INSTALL_PATH}
LIBRARY DESTINATION ${PYAUDI_INSTALL_PATH}
)

# Add the Python files.
install(FILES __init__.py _version.py test.py DESTINATION ${PYAUDI_INSTALL_PATH})




































