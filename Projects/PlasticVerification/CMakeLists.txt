file(GLOB header *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(PlasticVerification ${main_SRCS} ${header})

target_link_libraries(PlasticVerification pz)

