# @file neopz/External/FAD/TinyFad/CMakeLists.txt

file(GLOB headers *.h)

INSTALL(DIRECTORY Specializations
  DESTINATION ${NEOPZ_INCLUDE}/TinyFad/
  FILES_MATCHING PATTERN "*.h"
  )

file(GLOB sources *.cpp)

set(FAD_Tiny_SRCS ${headers} ${sources} PARENT_SCOPE)
