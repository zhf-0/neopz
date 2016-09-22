########### next target ###############
include_directories(
	"${PROJECT_SOURCE_DIR}/Projects/HCurlProjection"
)

file(GLOB headers *.h ../HCurlProjection/TPZMatHCurlProjection.h)
file(GLOB sources *.cpp ../HCurlProjection/TPZMatHCurlProjection.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(MFHCurlH1Lag ${main_SRCS})
target_link_libraries(MFHCurlH1Lag pz)