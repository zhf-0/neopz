########### next target ###############

file(GLOB headers *.h)
file(GLOB sources *.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(slepc-integration ${main_SRCS})
target_link_libraries(slepc-integration pz)