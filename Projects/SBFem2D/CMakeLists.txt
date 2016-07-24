file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain Main.cpp ${sources} ${headers})

add_executable(SBFem2D ${Mymain})

target_link_libraries(SBFem2D pz)

