# @file neopz/Pre/CMakeLists.txt  -- CMake file for the Pre module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if (USING_BOOST)
add_subdirectory(ParameterHandler)
endif()#USING_BOOST

if (BUILD_COMPLEX_PROJECTS)
file (GLOB to_remove_from_headers pzreadmeshhr.h)
list (REMOVE_ITEM headers ${to_remove_from_headers})
file (GLOB to_remove_from_sources pzreadmeshhr.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})
endif()

install(FILES ${headers} ${headersP} DESTINATION ${NEOPZ_INCLUDE})
set(pzpre_STAT_SRCS ${headers} ${sources} ${headersP} ${sourcesP} PARENT_SCOPE)
