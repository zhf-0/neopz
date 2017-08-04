# @file neopz/Material/REAL/Plasticity/CMakeLists.txt  -- CMake file for Plasticity Materials

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

add_subdirectory(PrincipalValues)

set(pzmaterial_Plasticity_SRCS ${headers} ${sources} ${pzmaterial_Plasticity_PrincipalValues_SRCS} PARENT_SCOPE)