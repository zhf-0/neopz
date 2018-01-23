########## next target ##############

file(GLOB headers *.h)
file(GLOB sources *.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(ModalAnalysis ${main_SRCS})
target_link_libraries(ModalAnalysis pz)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/circle.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/circlequad.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
