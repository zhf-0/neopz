file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(MixedPoroElasticity ${Mymain})

target_link_libraries(MixedPoroElasticity pz)

