# @file neopz/External/FAD/TinyFadET/CMakeLists.txt

file(GLOB headers *.h)

file(GLOB sources *.cpp)

set(FAD_TinyET_SRCS ${headers} ${sources} PARENT_SCOPE)
