########### next target ##############
include_directories("${PROJECT_SOURCE_DIR}/Projects/HCurlProjection")
file(GLOB headers *.h ../HCurlProjection/TPZMatHCurlProjection.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})
add_executable(Helmholtz2D ${main_SRCS})
target_link_libraries(Helmholtz2D pz)