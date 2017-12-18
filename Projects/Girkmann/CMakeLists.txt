file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(Girkmann ${Mymain})

target_link_libraries(Girkmann pz)

