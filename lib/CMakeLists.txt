# @file neopz/lib/CMakeList.txt  -- CMake file grouping pz modules and creating pz library. 

#Creating groups for organized display of files into PZ modules
source_group(Common FILES ${pzcommon_STAT_SRCS})
source_group(Save FILES ${pzsave_STAT_SRCS})
source_group(Util FILES ${pzutil_STAT_SRCS})
source_group(PerfUtil FILES ${pzperfutil_STAT_SRCS})
source_group(Matrix FILES ${pzmatrix_STAT_SRCS})
source_group("Linear Solvers" FILES ${pzlinearsolvers_STAT_SRCS})
source_group(Topology FILES ${pztopology_STAT_SRCS})
source_group(Geom FILES ${pzgeom_STAT_SRCS})
source_group("Special Maps" FILES ${pzspecialmaps_STAT_SRCS})
source_group(Mesh FILES ${pzmesh_STAT_SRCS})
source_group(Shape FILES ${pzshape_STAT_SRCS})
source_group(Refine FILES ${pzrefine_STAT_SRCS})
source_group(Material FILES ${pzmaterial_STAT_SRCS})
source_group(Analisys FILES ${pzanalysis_STAT_SRCS})
source_group(Multigrid FILES ${pzmultigrid_STAT_SRCS})
source_group(Post FILES ${pzpost_STAT_SRCS})
source_group(Frontal FILES ${pzfrontal_STAT_SRCS})
source_group("StructMatrix" FILES ${pzstrmatrix_STAT_SRCS})
source_group(Pre FILES ${pzpre_STAT_SRCS})
source_group("Substruct" FILES ${pzsubstruct_STAT_SRCS})
source_group(Integral FILES ${pzintegral_STAT_SRCS})
source_group(External FILES ${pzexternal_STAT_SRCS})
source_group(Publications FILES ${pzpublications_STAT_SRCS})
source_group(External\\Sloan FILES ${sloan_SRCS})
if(USING_FAD)
	source_group(External\\FAD FILES ${FAD_SRCS})
	source_group(External\\FAD\\Fad FILES ${FAD_fad_SRCS})
	source_group(External\\FAD\\TinyFad FILES ${FAD_Tiny_SRCS})
	source_group(External\\FAD\\TinyFadET FILES ${FAD_TinyET_SRCS})
endif(USING_FAD)

#creating pz lib and groups of the sources files
set(pzsources 
	${pzanalysis_STAT_SRCS} 
	${pzcommon_STAT_SRCS} 
	${pzexternal_STAT_SRCS} ${sloan_SRCS} ${FAD_SRCS} ${FAD_fad_SRCS} ${FAD_Tiny_SRCS} ${FAD_TinyET_SRCS}
	${pzfrontal_STAT_SRCS} 
	${pzgeom_STAT_SRCS} 
	${pzintegral_STAT_SRCS} 
        ${pzlinearsolvers_STAT_SRCS}
	${pzmaterial_STAT_SRCS} 
	${pzmatrix_STAT_SRCS} 
	${pzmesh_STAT_SRCS} 
	${pzmultigrid_STAT_SRCS} 
	${pzperfutil_STAT_SRCS} 
	${pzpre_STAT_SRCS} 
	${pzpost_STAT_SRCS}
	${pzrefine_STAT_SRCS} 
	${pzsave_STAT_SRCS} 
	${pzshape_STAT_SRCS} 
	${pzspecialmaps_STAT_SRCS} 
	${pzstrmatrix_STAT_SRCS} 
	${pzpublications_STAT_SRCS}
	${pzsubstruct_STAT_SRCS}
        ${pztopology_STAT_SRCS} 
	${pzutil_STAT_SRCS} 
	)

#creating target pz library and installing it
if (BUILD_PYTHON_BINDING AND NOT WIN32)
  add_library(pz SHARED ${pzsources})
else(BUILD_PYTHON_BINDING AND NOT WIN32)
  add_library(pz STATIC ${pzsources})
endif()

if(USING_FAD)
	target_compile_definitions(pz PRIVATE _AUTODIFF)
endif(USING_FAD)


#target_link_libraries(pz Threads::Threads)

if(USING_BOOST)
	#target_link_libraries(pz ${Boost_LIBRARIES})
	target_compile_definitions(pz PRIVATE USING_BOOST)
	if(BUILD_UNITTESTING)
		target_compile_definitions(pz PRIVATE BUILD_UNITTESTING)
	endif(BUILD_UNITTESTING)
endif(USING_BOOST)

if (USING_IPO)
  set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
endif(USING_IPO)

if (USING_LOG4CXX)
  #target_link_libraries(pz ${LOG4CXX_LIB})
#  if(WIN32)
#		target_link_libraries(pz odbc32.lib ws2_32.lib mswsock.lib)
#  endif(WIN32)
  target_compile_definitions(pz PRIVATE LOG4CXX)#should it be public?
endif(USING_LOG4CXX)

if(USING_TBB)
	#target_link_libraries(pz ${TBB_LIB})
	target_compile_definitions(pz PRIVATE USING_TBB)
endif(USING_TBB)

if (USING_OPENMP)
	#target_link_libraries(pz ${OPENMP_LIB})
	target_compile_definitions(pz PRIVATE USING_OPENMP)
endif(USING_OPENMP)

if (USING_LIKWID)
	#target_link_libraries(pz ${LIKWID_LIB})
	target_compile_definitions(pz PRIVATE USING_LIKWID)
endif(USING_LIKWID)

if (USING_LIBNUMA)
	#target_link_libraries(pz ${LIBNUMA_LIB})
	target_compile_definitions(pz PRIVATE USING_LIBNUMA)
endif(USING_LIBNUMA)

if(USING_MATLAB_ENGINE)
	target_compile_definitions(pz PRIVATE USING_MATLAB)
endif(USING_MATLAB_ENGINE)

if(USING_PAPI)
	#target_link_libraries(pz ${PAPI_LIB})
	target_compile_definitions(pz PRIVATE USING_PAPI)
endif(USING_PAPI)

if(USING_HWLOC)
	#target_link_libraries(pz ${HWLOC_LIB})
	target_compile_definitions(pz PRIVATE USING_HWLOC)
endif(USING_HWLOC)

if (USING_METIS)
	#target_link_libraries(pz ${METIS_LIB})
	target_compile_definitions(pz PRIVATE USING_METIS)
endif(USING_METIS)

if(USING_QWT)
	target_compile_definitions(pz PRIVATE USING_QWT)
endif(USING_QWT)

if(USING_OPENSSL)
	target_link_libraries(pz crypto)###what?
	target_compile_definitions(pz PRIVATE USING_OPENSSL)
endif(USING_OPENSSL)

if(USING_MKL)
	#target_link_libraries(pz ${MKL_LIBRARIES})
	target_compile_definitions(pz PUBLIC USING_MKL)
endif(USING_MKL)

if(USING_LAPACK)
	#target_link_libraries(pz ${LAPACK_LIBRARIES})
	target_compile_definitions(pz PRIVATE USING_LAPACK)
endif(USING_LAPACK)

if(USING_BLAS)
	#target_link_libraries(pz ${BLAS_LIBRARIES})
	target_compile_definitions(pz PRIVATE USING_BLAS)
endif(USING_BLAS)


if(USING_PETSC)
	#target_link_libraries(pz PETSC::petsc)
	target_compile_definitions(pz PRIVATE USING_PETSC)
endif(USING_PETSC)

if(USING_SLEPC)
	#target_link_libraries(pz SLEPC::slepc)
	target_compile_definitions(pz PRIVATE USING_SLEPC)
endif(USING_SLEPC)


target_link_libraries(pz ${PZ_LIBRARIES})
## inserting definitions
#specifying the path to neopz source code
set(PZSOURCEDIR ${PROJECT_SOURCE_DIR})
#specify where the refinement patterns can be found (default behaviour when running the code)
set(REFPATTERNDIR ${PROJECT_SOURCE_DIR}/Refine/RefPatterns)

# first the common flags, then release flags,
# and finally the debug flags
#####################################################################
#This should not be public, but private.							#
#We must urgently remove PZSOURCEDIR and REFPATTERNDIR usages from  #
#our projects.														#
#####################################################################

target_compile_definitions(pz PUBLIC
	PZSOURCEDIR=\"${PZSOURCEDIR}\"
	REFPATTERNDIR=\"${REFPATTERNDIR}\"
	HAVE_STRUCT_TIMESPEC
    $<$<CONFIG:Debug>:PZDEBUG ZERO_INTERNAL_RESIDU MAKEINTERNAL>
    $<$<CONFIG:Release>:NODEBUG >
)

if(APPLE)
	target_compile_definitions(pz PRIVATE
		MACOSX	
	)
endif(APPLE)

if (NOT WIN32)
	target_compile_options(pz PRIVATE "-fPIC" "-Wno-narrowing")
endif(NOT WIN32)

install(TARGETS pz DESTINATION ${NEOPZ_LIB})

##############PZ configuration##############

if(USING_NEW_SKYLMAT)
	target_compile_definitions(pz PRIVATE USING_NEW_SKYLMAT)
endif(USING_NEW_SKYLMAT)

#To building several examples for training and experimentation PZ
if(BUILD_TUTORIAL)
	target_compile_definitions(pz PRIVATE BUILD_TUTORIAL)
endif()

#To building several examples for training and experimentation PZ

if(BUILD_TUTORIAL)
	target_compile_definitions(pz PRIVATE BUILD_TUTORIAL)
endif()


#defining STATE and REAL types


if (REAL_TYPE STREQUAL "double")
  target_compile_definitions(pz PUBLIC REALdouble)
  set(REAL_TYPE_DEF "REALdouble")
elseif (REAL_TYPE STREQUAL "float")
  target_compile_definitions(pz PUBLIC REALfloat)
  set(REAL_TYPE_DEF "REALfloat")
elseif (REAL_TYPE STREQUAL "long double")
  target_compile_definitions(pz PUBLIC REALlongdouble)
  set(REAL_TYPE_DEF "REALlongdouble")
elseif (REAL_TYPE STREQUAL "pzfpcounter")
  target_compile_definitions(pz PUBLIC REALpzfpcounter)
  set(REAL_TYPE_DEF "REALpzfpcounter")
else()
  message (FATAL_ERROR "ERROR: Please specify a valid type for REAL. ${REALTypesAre}")
endif()
message (STATUS "Using REAL type = ${REAL_TYPE}")

if (STATE_TYPE STREQUAL "double")
  target_compile_definitions(pz PUBLIC STATEdouble)
  set(STATE_TYPE_DEF "STATEdouble")
elseif (STATE_TYPE STREQUAL "float")
  target_compile_definitions(pz PUBLIC STATEfloat)
  set(STATE_TYPE_DEF "STATEfloat")
elseif (STATE_TYPE STREQUAL "long double")
  target_compile_definitions(pz PUBLIC STATElongdouble)
  set(STATE_TYPE_DEF "STATElongdouble")
elseif (STATE_TYPE STREQUAL "complex<float>")
  target_compile_definitions(pz PUBLIC STATEcomplexf)
  target_compile_definitions(pz PUBLIC STATE_COMPLEX)
  set(STATE_TYPE_DEF "STATEcomplexf")
elseif (STATE_TYPE STREQUAL "complex<double>")
  target_compile_definitions(pz PUBLIC STATEcomplexd)
  target_compile_definitions(pz PUBLIC STATE_COMPLEX)
  set(STATE_TYPE_DEF "STATEcomplexd")
elseif (STATE_TYPE STREQUAL "complex<long double>")
  target_compile_definitions(pz PUBLIC STATEcomplexld)
  target_compile_definitions(pz PUBLIC STATE_COMPLEX)
  set(STATE_TYPE_DEF "STATEcomplexld")
else()
  message (FATAL_ERROR "ERROR: Please specify a valid type for STATE. ${STATETypesAre}")
endif()
message (STATUS "Using STATE type = ${STATE_TYPE}")

# Some definitions for windows
if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on config.h")
	target_compile_definitions(pz PUBLIC VC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	target_compile_definitions(pz PUBLIC _USE_MATH_DEFINES)
	#disabling VC warnings
	target_compile_definitions(pz PUBLIC _SCL_SECURE_NO_WARNINGS)
	target_compile_definitions(pz PUBLIC _CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	#disabling STL warning for std::_Vector_alloc when LOG4CXX is enabled
	if (USING_LOG4CXX)
		target_compile_definitions(pz PUBLIC "/wd4251")
	endif()
	#define use of pthread static lib.
	target_compile_definitions(pz PUBLIC PTW32_STATIC_LIB)
    #enabling /bigobj
    target_compile_definitions(pz PUBLIC "/bigobj")

    FOREACH(FLAG_TYPE EXE MODULE SHARED)
        STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
        STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
        STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
        STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
        SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
        MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
    ENDFOREACH ()

endif()


# Begin of config.h part
# This file contains all definitions used to compile the library pz
# this will be saved at include diretory in the install target
set (CONFIG_H "${PROJECT_BINARY_DIR}/Common/config.h")
configure_file ("${PROJECT_SOURCE_DIR}/Common/config.h.in" ${CONFIG_H})
# add it to include search path
target_include_directories(pz PUBLIC "${PROJECT_BINARY_DIR}/Common")
# add the .h file in the install target
install (FILES ${CONFIG_H} DESTINATION ${NEOPZ_INCLUDE})
# End of config.h part


#creating pzmatrix lib -- It is comment now, because it is not needed when the pzlib is installed.
#set(pzmatrixsources 
#	${pzcommon_STAT_SRCS} ${pzutil_STAT_SRCS} ${pzperfutil_STAT_SRCS} ${pzsave_STAT_SRCS} 
#	${pzmatrix_STAT_SRCS} ${pzlinearsolvers_STAT_SRCS}
#	)

#add_library(pzmatrix ${pzmatrixsources})
#install(TARGETS pzmatrix ARCHIVE DESTINATION lib)

#Verifying the content of the strings
#file(APPEND "${LogFile}" "${refpatterns_STAT_SRCS}   ..  \n\n")
