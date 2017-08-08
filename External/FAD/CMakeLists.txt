# @file neopz/External/FAD/CMakeLists.txt  -- CMake file for the FAD external module

file(GLOB headers *.h)

INSTALL(DIRECTORY ./
  DESTINATION ${NEOPZ_INCLUDE}
  FILES_MATCHING PATTERN "*.h"
  )
INSTALL(DIRECTORY utils
  DESTINATION ${NEOPZ_INCLUDE}
  FILES_MATCHING PATTERN "*.h"
  )
INSTALL(DIRECTORY Fad
  DESTINATION ${NEOPZ_INCLUDE}
  FILES_MATCHING PATTERN "*.h"
  )
INSTALL(DIRECTORY TinyFad
  DESTINATION ${NEOPZ_INCLUDE}
  FILES_MATCHING PATTERN "*.h"
  )
INSTALL(DIRECTORY TinyFadET
  DESTINATION ${NEOPZ_INCLUDE}
  FILES_MATCHING PATTERN "*.h"
  )
file(GLOB sources *.cpp)

add_subdirectory(Fad)
add_subdirectory(TinyFadET)
add_subdirectory(TinyFad)

set(FAD_fad_SRCS ${FAD_fad_SRCS} PARENT_SCOPE)

set(FAD_Tiny_SRCS ${FAD_Tiny_SRCS} PARENT_SCOPE)

set(FAD_TinyET_SRCS ${FAD_TinyET_SRCS} PARENT_SCOPE)

set(FAD_SRCS ${headers} ${sources} PARENT_SCOPE)

