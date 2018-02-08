# @file neopz/Pre/CMakeLists.txt  -- CMake file for the Pre module

file(GLOB headers *.h std_cxx14/*.h)
file(GLOB sources *.cpp)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

set(headersP ${headers} PARENT_SCOPE)
set(sourcesP ${sources} PARENT_SCOPE)