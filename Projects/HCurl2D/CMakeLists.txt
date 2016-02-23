########### next target ###############
file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})
add_executable(HCurl2D ${main_SRCS})
target_link_libraries(HCurl2D pz)