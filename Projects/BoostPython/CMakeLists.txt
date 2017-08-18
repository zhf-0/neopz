file(COPY testPython.py DESTINATION .)

if (BUILD_BOOST_PYTHON)
	# Creating lib with c++ functions to be parsed to python
  source_group(BoostPythonFunc FILES PZPythonFunc.cpp PZPythonFunc.h)
 	add_library(PZPythonFunc PZPythonFunc.cpp PZPythonFunc.h)
 	target_link_libraries(PZPythonFunc pz)

 	# Python module where the parsing is defined
	PYTHON_ADD_MODULE(pyz pyz.cpp)
	target_link_libraries(pyz ${PYTHON_LIBRARIES} PZPythonFunc  ${Boost_LIBRARIES})

endif()