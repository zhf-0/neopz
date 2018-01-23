########## next target ##############
include_directories(
	"${PROJECT_SOURCE_DIR}/Projects/ModalAnalysis"
)

file(GLOB headers *.h ../ModalAnalysis/TPZMatModalAnalysis.h)
file(GLOB sources *.cpp ../ModalAnalysis/TPZMatModalAnalysis.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(NewModalAnalysis ${main_SRCS})
target_link_libraries(NewModalAnalysis pz)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/coarseMesh.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fineMesh.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/refinedInterfaceMesh.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
