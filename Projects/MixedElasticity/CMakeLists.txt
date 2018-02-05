file(GLOB headers *.h)
file(GLOB sources meshgen.cpp pzmixedelasmat.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(MixedElasticity ${Mymain})

target_link_libraries(MixedElasticity pz)

