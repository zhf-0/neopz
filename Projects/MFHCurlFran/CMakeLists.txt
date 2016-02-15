########### next target ###############
file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})
add_executable(MFHCurlFran ${main_SRCS})
target_link_libraries(MFHCurlFran pz)