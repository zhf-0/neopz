file(GLOB headers *.h)
file(GLOB sources meshgeng.cpp pzmixedelasmatg.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(GirktestSquare ${Mymain})

target_link_libraries(GirktestSquare pz)

set(MymainG mainGirkmann.cpp ${sources} ${headers})

add_executable(GirktestGirkmann ${MymainG})

target_link_libraries(GirktestGirkmann pz)

configure_file(Girkmann.msh Girkmann.msh COPYONLY)
