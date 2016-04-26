########### next target ###############
include_directories(
	"${PROJECT_SOURCE_DIR}/Projects/HCurl2D"
)

file(GLOB headers *.h ../HCurl2D/TPZMatHCurl2D.h)
file(GLOB sources *.cpp ../HCurl2D/TPZMatHCurl2D.cpp)


set(main_SRCS ${sources} ${headers})
add_executable(MFHCurlH1 ${main_SRCS})
target_link_libraries(MFHCurlH1 pz)