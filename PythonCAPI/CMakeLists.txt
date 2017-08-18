# @file neopz/PythonCAPI/CMakeLists.txt  -- CMake file for the PythonCAPI module

file(GLOB headers *.h)
file(GLOB sources *.cpp)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})
set(pzpythoncapi_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
