# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.7 FATAL_ERROR)

project (PZ)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#this is the folder where additional findMYLIB.cmake modules must be put
#in order to use find_package(MYLIB)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(PzDependentOption)#macro for dependent options
#Specify where the headers files must be finded
include_directories(
	"${PROJECT_SOURCE_DIR}/Python"
	"${PROJECT_SOURCE_DIR}/Util"
	"${PROJECT_SOURCE_DIR}/PerfUtil"
	"${PROJECT_SOURCE_DIR}/Common"
	"${PROJECT_SOURCE_DIR}/Save"
	"${PROJECT_SOURCE_DIR}/Matrix"
	"${PROJECT_SOURCE_DIR}/Topology"
	"${PROJECT_SOURCE_DIR}/Geom"
	"${PROJECT_SOURCE_DIR}/SpecialMaps"
	"${PROJECT_SOURCE_DIR}/Refine"
	"${PROJECT_SOURCE_DIR}/Shape"
	"${PROJECT_SOURCE_DIR}/Material"
	"${PROJECT_SOURCE_DIR}/Material/REAL"
	"${PROJECT_SOURCE_DIR}/Material/REAL/Plasticity"
	"${PROJECT_SOURCE_DIR}/Material/Complex"
	"${PROJECT_SOURCE_DIR}/Multigrid"
	"${PROJECT_SOURCE_DIR}/Mesh"
	"${PROJECT_SOURCE_DIR}/Mesh/REAL"
	"${PROJECT_SOURCE_DIR}/StrMatrix"
	"${PROJECT_SOURCE_DIR}/Integral"
	"${PROJECT_SOURCE_DIR}/Frontal"
	"${PROJECT_SOURCE_DIR}/Pre"
	"${PROJECT_SOURCE_DIR}/Post"
	"${PROJECT_SOURCE_DIR}/Analysis"
	"${PROJECT_SOURCE_DIR}/SubStruct"
	"${PROJECT_SOURCE_DIR}/LinearSolvers"
	"${PROJECT_SOURCE_DIR}/External"
	"${PROJECT_SOURCE_DIR}/External/sloan"
	"${PROJECT_SOURCE_DIR}/Publications"
)

#Defining variables to cmake

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include /usr/lib/atlas-base/ /usr/include/atlas /usr/local/atlas/include /usr/local/atlas/lib )

# INSTALL part
#To install neopz/include
set(NEOPZ_INCLUDE pzlib/include)
#To install neopz library
set(NEOPZ_LIB pzlib/lib)
############################################
############################################
############################################
##############  All options   ##############
#To enable Python bindings and create the library
option(BUILD_PYTHON_BINDING "To build the Python bindings" OFF)
#enabling Boost library
option(USING_BOOST "Whether the Boost package will be linked with the library" OFF)
#To enable unit test using boost UnitTest framework
option(BUILD_UNITTESTING "Whether the Boost UnitTest will be linked in" OFF)
#To use differentiation with FAD
option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
#enabling compiler interprocedural optimizations (IPO)
option(USING_IPO "Whether link and assembly support for compiler IPO will be activated or not (may cause linking to be very slow)." OFF)
#enabling METIS library
option(USING_METIS "Whether the METIS library will be linked in" OFF)
#enabling OPENSSL library (used for MD5 computation)
option(USING_OPENSSL "Whether the OPENSSL library will be linked in" OFF)
#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
#enabling OpenMP runtime and library
option(USING_OPENMP "Whether the OpenMP library and runtime will be linked in" OFF)
#enabling LIKWID library
option(USING_LIKWID "Whether the likwid library will be linked in" OFF)
#enabling NUMA library
option(USING_LIBNUMA "Whether the Lib NUMA library will be linked in" OFF)
#enabling MatLab library
option(USING_MATLAB_ENGINE "Whether MatLab Engine will be linked in" OFF)
#enabling BLAS library
option(USING_BLAS "Whether the Basic Linear Algebra Subprograms (BLAS) library will be linked in" OFF)
pz_dependent_option(BLA_VENDOR STRING "Which BLAS to use" ""
                       USING_BLAS "")#will be used in order to specify BLAS's vendor. not ready yet.
#enabling QWT library
option(USING_QWT "Whether the QWT library will be linked in" OFF)
#enabling MKL library
option(USING_MKL "Whether the MKL library will be linked in" OFF)
#enabling HWLOC library
option(USING_HWLOC "Whether the HWLOC library will be linked in" OFF)
#enabling PAPI library
option(USING_PAPI "Whether the PAPI library will be linked in" OFF)
#enabling LAPACK library
option(USING_LAPACK "Whether the LAPACK library will be linked in" OFF)
##############PZ configuration##############
#enabling the new skyline matrix format
option(USING_NEW_SKYLMAT "Enable the new skyline matrix format" OFF)
#To building several examples for training and experimentation PZ
option(BUILD_TUTORIAL "To build examples for training or teaching PZ" OFF)
#To building several examples for training and experimentation PZ
option(BUILD_BROKEN "To build projects which are broken" OFF)
#To building all the specific projects developed with PZ environment
option(BUILD_PROJECTS "To build specific projects." OFF)
#To building all the performance projects developed with PZ environment
option(BUILD_PERF_TESTS "To build the performance test." OFF)
#To building all the plasticity materials
option(BUILD_PLASTICITY_MATERIALS "To build plasticity materials." OFF)
############################################
############################################
############################################
##########option-dependent actions##########
#To enable Python bindings and create the library
if (BUILD_PYTHON_BINDING)
	find_package(SWIG "2.0" REQUIRED)
	find_package(PythonLibs REQUIRED)
endif()

#ENABLING EXTERNAL LIBRARIES

#enabling pthread library
#pthread's location for windows
if(WIN32)
	set(PTHREADS_ROOT_DIR ./../externallibs)
	find_package(PThreads-win32 REQUIRED)
endif(WIN32)

set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

IF(APPLE)
   # Mac OS X specific code
   find_library(ACCELERATE_LIB Accelerate)
#   SET(OperatingSystem "Mac OS X")
ENDIF(APPLE)


#To enable unit test using boost UnitTest framework
if(BUILD_UNITTESTING)
	set(USING_BOOST ON CACHE BOOL "Whether the Boost package will be linked with the library" FORCE)
endif()

if(USING_BOOST)
	IF (WIN32)
		set(Boost_USE_STATIC_LIBS   ON)
		set(Boost_USE_MULTITHREADED ON)
		set(BOOST_ROOT "../externallibs/lib/boost_1_65_0")	
	ENDIF()
	find_package( Boost COMPONENTS unit_test_framework REQUIRED )
   	include_directories(${Boost_INCLUDE_DIR})
endif(USING_BOOST)


if(USING_FAD)
	include_directories(External/FAD External/FAD/Fad External/FAD/TinyFad External/FAD/TinyFadET)
endif(USING_FAD)

#enabling to use LOG4CXX library
if(USING_LOG4CXX)
	find_path(LOG4CXX_INCLUDE log4cxx PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../log4cxx ${PROJECT_SOURCE_DIR}/../log4cxx_src/include ${PROJECT_SOURCE_DIR}/../externallibs/log4cxx_src ${PROJECT_SOURCE_DIR}/../externallibs/log4cxx_src/include)
	find_library(LOG4CXX_LIB_RELEASE
		         NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a
				 PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib ./../externallibs/lib/Release
				)
	find_library(LOG4CXX_LIB_DEBUG
				 NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a       # Change name here if debug lib has different name.
				 PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib ./../externallibs/lib/Debug
				)
	set(LOG4CXX_LIB
	    debug ${LOG4CXX_LIB_DEBUG}
	    optimized ${LOG4CXX_LIB_RELEASE}
	   )
	if(NOT LOG4CXX_INCLUDE)
		set (LOG4CXX_INCLUDE "" CACHE PATH "Directory where log4cxx header files can be found")
		message(FATAL_ERROR "log4cxx's headers not found")
	else()
		include_directories(${LOG4CXX_INCLUDE})
	endif()

	if(NOT LOG4CXX_LIB)
		set (LOG4CXX_LIB "" CACHE PATH "Directory where the log4cxx library can be found")
		message(FATAL_ERROR "log4cxx lib not found")
	endif()
endif(USING_LOG4CXX)

#enabling METIS library
if(USING_METIS)
	find_path(METIS_INCLUDE metis.h PATHS ${SEARCH_DIRS})
	find_library(METIS_LIB NAMES libmetis.dylib metis.lib libmetis.so libmetis.a PATHS ${SEARCH_DIRS})
	if(NOT METIS_INCLUDE)
		set (METIS_INCLUDE "" CACHE PATH "Directory where metis header files can be found")
		message(FATAL_ERROR "metis's headers not found")
	else()
		include_directories(${METIS_INCLUDE})
	endif()
	if(NOT METIS_LIB)
		set (METIS_LIB "" CACHE PATH "Directory where the metis library can be found")
		message(FATAL_ERROR "metis lib not found")
	endif()
endif(USING_METIS)

#enabling TBB library
if(USING_TBB)
	find_path(TBB_INCLUDE tbb.h PATHS ${SEARCH_DIRS}
	/opt/intel/tbb/include/tbb /usr/include/tbb /opt/local/include/tbb)
	find_library(TBB_LIB NAMES libtbb.dylib libtbb.so PATHS ${SEARCH_DIRS}  /opt/intel/tbb/lib )
	if(NOT TBB_INCLUDE)
		set (TBB_INCLUDE "" CACHE PATH "Directory where tbb header files can be found")
		message(FATAL_ERROR "tbb's headers not found")
	else()
		include_directories(/opt/local/include)
		include_directories(${TBB_INCLUDE})
	endif()
	if(NOT TBB_LIB)
		set (TBB_LIB "" CACHE PATH "Directory where the tbb library can be found")
		message(FATAL_ERROR "tbb lib not found")
	endif()
endif(USING_TBB)

#enabling OpenMP runtime and library
if(USING_OPENMP)
	find_path(OPENMP_INCLUDE omp.h PATHS ${SEARCH_DIRS} /opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.2.0/include)
	find_library(OPENMP_LIB NAMES libgomp.dylib PATHS ${SEARCH_DIRS}  /opt/local/lib/gcc5)
	if(NOT OPENMP_INCLUDE)
		set (OPENMP_INCLUDE "" CACHE PATH "Directory where omp header files can be found")
	else()
		include_directories(${OPENMP_INCLUDE})
	endif()
	if(NOT OPENMP_LIB)
		set (TBB_OPENMP "" CACHE PATH "Directory where the OpenMP library can be found")
		message(FATAL_ERROR "openmp lib not found")
	endif()
endif(USING_OPENMP)

#enabling LIKWID library
if(USING_LIKWID)
	find_path(LIKWID_INCLUDE likwid.h PATHS ${SEARCH_DIRS})
	find_library(LIKWID_LIB NAMES liblikwid.a PATHS ${SEARCH_DIRS})
	if(NOT LIKWID_INCLUDE)
		set (LIKWID_INCLUDE "" CACHE PATH "Directory where likwid header files can be found")
		message(FATAL_ERROR "likwid's headers not found")
	else()
		include_directories(${LIKWID_INCLUDE})
	endif()
	if(NOT LIKWID_LIB)
		set (LIKWID_LIB "" CACHE PATH "Directory where the likwid library can be found")
		message(FATAL_ERROR "likwid lib not found")
	endif()
endif(USING_LIKWID)

#enabling numa library
if(USING_LIBNUMA)
	find_path(LIBNUMA_INCLUDE numa.h PATHS ${SEARCH_DIRS})
	find_library(LIBNUMA_LIB NAMES libnuma.so PATHS ${SEARCH_DIRS})
	if(NOT LIBNUMA_INCLUDE)
		set (LIBNUMA_INCLUDE "" CACHE PATH "Directory where libnuma header files can be found (consider installing libnuma-dev)")
		message(FATAL_ERROR "Lib NUMA's headers not found")
	else()
		include_directories(${LIBNUMA_INCLUDE})
	endif()
	if(NOT LIBNUMA_LIB)
		set (LIBNUMA_LIB "" CACHE PATH "Directory where the libnuma library can be found")
		message(FATAL_ERROR "Lib NUMA lib not found")
	endif()
endif(USING_LIBNUMA)

#enabling MatLab library
if(USING_MATLAB_ENGINE)
	find_package(MATLAB REQUIRED)
	include_directories(${Matlab_ROOT_DIR}/extern/include)
endif(USING_MATLAB_ENGINE)

#enabling BLAS library
if(USING_BLAS)
	find_package(BLAS REQUIRED)
endif(USING_BLAS)

#enabling LAPACK library
if(USING_LAPACK)
	if(APPLE)
	else()
		find_path(LAPACK_INCLUDE clapack.h PATHS ${SEARCH_DIRS})
		find_library(LAPACK_LIB NAMES libclapack.dylib PATHS ${SEARCH_DIRS})
		if(NOT LAPACK_INCLUDE)
			set (LAPACK_INCLUDE "" CACHE PATH "Directory where LAPACK header files can be found")
		else()
			include_directories(${LAPACK_INCLUDE})
		endif()
		if(NOT LAPACK_LIB)
			set (LAPACK_LIB "" CACHE PATH "Directory where the LAPACK library can be found")
		endif()
	endif(APPLE)
endif(USING_LAPACK)

#enabling PAPI library
if(USING_PAPI)
		find_path(PAPI_INCLUDE NAMES papi.h PATHS ${SEARCH_DIRS})
		find_library(PAPI_LIB NAMES libpapi.a libpapi.so PATHS ${SEARCH_DIRS})
		if(NOT PAPI_INCLUDE)
			set (PAPI_INCLUDE "" CACHE PATH "Directory where papi.h can be found")
		endif()
		if(NOT PAPI_LIB)
			set (PAPI_LIB "" CACHE PATH "Directory where the papi library can be found")
		endif()
endif(USING_PAPI)

#enabling HWLOC library
if(USING_HWLOC)
		find_path(HWLOC_INCLUDE NAMES hwloc.h PATHS ${SEARCH_DIRS})
		find_library(HWLOC_LIB NAMES hwloc.a libhwloc.so libhwloc.dylib PATHS ${SEARCH_DIRS})
		if(NOT HWLOC_INCLUDE)
			set (HWLOC_INCLUDE "" CACHE PATH "Directory where hwloc.h can be found")
		endif()
		if(NOT HWLOC_LIB)
			set (HWLOC_LIB "" CACHE PATH "Directory where the hwloc library can be found")
		endif()
endif(USING_HWLOC)



#enabling MKL library
if(USING_MKL)

		if (WIN32)
			set(MKL_ROOT "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows")
			find_path(MKL_INCLUDE NAMES mkl.h PATHS ${MKL_ROOT}/mkl/include)
			find_library(MKL_LIB_INTEL NAMES mkl_rt.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
			find_library(MKL_LIB_CORE NAMES mkl_core.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
			find_library(MKL_LIB_THREAD NAMES mkl_intel_thread.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
			find_library(COMPOSER_OMP NAMES libiomp5md.lib PATHS ${MKL_ROOT}/compiler/lib/intel64_win)
		else()
			find_path(MKL_INCLUDE NAMES mkl.h PATHS ${SEARCH_DIRS} /opt/intel/mkl/include /softwares/intel/mkl/include)
			find_library(MKL_LIB_INTEL NAMES libmkl_intel_lp64.so libmkl_intel_lp64.dylib PATHS
					${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64)
			find_library(MKL_LIB_CORE NAMES libmkl_core.so libmkl_core.dylib PATHS ${SEARCH_DIRS}
					/opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
			find_library(MKL_LIB_THREAD NAMES libmkl_intel_thread.so libmkl_intel_thread.dylib PATHS
					${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
			find_library(COMPOSER_OMP NAMES libiomp5.so libiomp5.dylib PATHS ${SEARCH_DIRS}
					/opt/intel/composer_xe/compiler/lib
					/opt/intel/composerxe/lib/intel64
					/opt/intel/lib /softwares/intel/lib/intel64)
			#find_path(OMP_INCLUDE NAMES omp.h PATHS ${SEARCH_DIRS} /opt/intel/composer_xe_2013_sp1.1.103/compiler/include)
			#include_directories(${OMP_INCLUDE})
		endif()

		if(NOT MKL_INCLUDE)
			set (MKL_INCLUDE "" CACHE PATH "Directory where mkl.h can be found")
		else()
			include_directories(${MKL_INCLUDE})
		endif()
		if(NOT MKL_LIB)
			set (MKL_LIB_INTEL "" CACHE PATH "Directory where the mkl library can be found")
		endif()
endif(USING_MKL)


#enabling QWT library
if(USING_QWT)
	if (WIN32)
		set(QWT_ROOT "C:/Qwt-6.1.3")
	else()
		set(QWT_ROOT "" CACHE PATH "Root directory where the QWT library is istalled")
	endif()
endif(USING_QWT)

##############PZ configuration##############
#To select type of REAL
set(REALTypesAre "Types are: double, float, long double, pzfpcounter")
set(REAL_TYPE "double" CACHE STRING "${REALTypesAre}")



# Begin of config.h part
# This file contains all definitions used to compile the library pz
# this will be saved at include diretory in the install target
set (CONFIG_H "${PROJECT_BINARY_DIR}/Common/config.h")
configure_file ("${PROJECT_SOURCE_DIR}/Common/config.h.in" ${CONFIG_H})
# add it to include search path
include_directories("${PROJECT_BINARY_DIR}/Common")
# add the .h file in the install target
install (FILES ${CONFIG_H} DESTINATION ${NEOPZ_INCLUDE})
# End of config.h part

#increment all the directories where the sources files is
if(NOT BUILD_COMPLEX_PROJECTS)
    if(BUILD_PERF_TESTS)
        add_subdirectory(PerfTests)
    endif()
    if (BUILD_PYTHON_BINDING)
        add_subdirectory(Python)
    endif()
endif()

add_subdirectory(PerfUtil)
add_subdirectory(Util)
add_subdirectory(Common)
add_subdirectory(Save)
add_subdirectory(Integral)
add_subdirectory(LinearSolvers)
add_subdirectory(Matrix)
add_subdirectory(Topology)
add_subdirectory(Geom)
add_subdirectory(SpecialMaps)
add_subdirectory(Shape)
add_subdirectory(Refine)
add_subdirectory(External)
add_subdirectory(Material)
add_subdirectory(Mesh)
add_subdirectory(Analysis)
add_subdirectory(Multigrid)
add_subdirectory(Post)
add_subdirectory(Frontal)
add_subdirectory(StrMatrix)
add_subdirectory(Pre)
add_subdirectory(SubStruct)
add_subdirectory(Publications)


add_subdirectory(lib)

if(BUILD_PROJECTS OR BUILD_TUTORIAL)
	add_subdirectory(Projects)
endif()


if(BUILD_UNITTESTING)
	add_subdirectory(UnitTest_PZ)
endif()

#Enabling tests
ENABLE_TESTING()
include(CTest)

#To insert cmakelist files as a group
if(BUILD_TUTORIAL)
	#Grouping CMakeList.txt for pz Projects - UnitTest_PZ
	file(GLOB_RECURSE makefiles *Lists.txt)

	#creating text files group
	source_group(CMakefilesPZ FILES ${makefilespz})
	source_group(CMakefilesProjects FILES ${makefilesproject})
	source_group(CMakefilesUnitTest FILES ${makefilesunittest})
	add_custom_target(CMAKE_Files SOURCES ${makefiles})
endif(BUILD_TUTORIAL)

#To insert doxygen documentation files
file(GLOB_RECURSE doxyfiles Documentation/*.h ./dox*.h)
add_custom_target(DOXYGEN_Files SOURCES  ${doxyfiles})
