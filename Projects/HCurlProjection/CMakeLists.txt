########### next target ##############
file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})
add_executable(HCurlProjection ${main_SRCS})
target_link_libraries(HCurlProjection pz)