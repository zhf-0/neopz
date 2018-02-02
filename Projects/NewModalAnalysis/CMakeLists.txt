########## next target ##############

file(GLOB headers *.h)
file(GLOB sources *.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(NewModalAnalysis ${main_SRCS})
target_link_libraries(NewModalAnalysis pz)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/veryCoarseMesh.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/coarseMesh.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fineMesh.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/refinedInterface.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/wr90.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/refinedWR90.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
