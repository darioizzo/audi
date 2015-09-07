# Add the unit test library.
# NOTE: probably we should have a check in the build system to determine
# whether Boost libraries are static or dynamic.
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

MACRO(ADD_AUDI_TESTCASE arg1)
	ADD_EXECUTABLE(${arg1} helpers.hpp ${arg1}.cpp)
	TARGET_LINK_LIBRARIES(${arg1} ${MANDATORY_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
	ADD_TEST(${arg1} ${arg1})
ENDMACRO(ADD_AUDI_TESTCASE)

ADD_AUDI_TESTCASE(gdual)
ADD_AUDI_TESTCASE(functions)
