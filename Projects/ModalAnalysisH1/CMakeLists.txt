########### next target ###############
include_directories(
	"${PROJECT_SOURCE_DIR}/Projects/HCurlProjection"
)

file(GLOB headers *.h ../HCurlProjection/TPZMatHCurlProjection.h)
file(GLOB sources *.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(ModalAnalysisH1 ${main_SRCS})
target_link_libraries(ModalAnalysisH1 pz)