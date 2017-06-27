########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(GPU_poisson ${sources} ${header})

target_link_libraries(GPU_poisson pz ${CBLAS_LIB} ${LAPACK_LIB}) 
