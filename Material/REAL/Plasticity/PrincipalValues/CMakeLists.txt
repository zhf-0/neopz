# @file neopz/Material/REAL/Plasticity/PrincipalValues/CMakeLists.txt  -- CMake file for Plasticity Materials Implemented in Principal Values

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzmaterial_Plasticity_PrincipalValues_SRCS ${headers} ${sources} PARENT_SCOPE)