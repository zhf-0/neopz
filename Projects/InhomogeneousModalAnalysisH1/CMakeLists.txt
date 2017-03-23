########### next target ###############
include_directories(
"${PROJECT_SOURCE_DIR}/Projects/HCurlProjection"
)

file(GLOB headers *.h ../HCurlProjection/TPZMatHCurlProjection.h)
file(GLOB sources *.cpp ../HCurlProjection/TPZMatHCurlProjection.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(InhomogeneousModalAnalysisH1 ${main_SRCS})
target_link_libraries(InhomogeneousModalAnalysisH1 pz)
