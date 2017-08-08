########### next target ###############
include_directories(
	"${PROJECT_SOURCE_DIR}/Projects/HCurlProjection"
)

file(GLOB headers *.h ../HCurlProjection/TPZMatHCurlProjection.h)
file(GLOB sources *.cpp ../HCurlProjection/TPZMatHCurlProjection.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(ModalAnalysis ${main_SRCS})
target_link_libraries(ModalAnalysis pz)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/circle.geo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/circle.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
