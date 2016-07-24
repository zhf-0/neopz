file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain Main.cpp ${sources} ${headers})

add_executable(Elasticity2DTensor ${Mymain})

target_link_libraries(Elasticity2DTensor pz)

