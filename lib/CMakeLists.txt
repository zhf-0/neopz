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
source_group(External\\FAD FILES ${FAD_SRCS})
source_group(External\\FAD\\Fad FILES ${FAD_fad_SRCS})
source_group(External\\FAD\\TinyFad FILES ${FAD_Tiny_SRCS})
source_group(External\\FAD\\TinyFadET FILES ${FAD_TinyET_SRCS})

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

target_link_libraries(pz Threads::Threads)

if(APPLE)
	target_link_libraries(pz ${ACCELERATE_LIB})
endif(APPLE)

if(USING_BOOST)
	target_link_libraries(pz ${Boost_LIBRARIES})
endif(USING_BOOST)

if (USING_IPO)
  set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
endif(USING_IPO)

if (USING_LOG4CXX)
  target_link_libraries(pz ${LOG4CXX_LIB})
  if(WIN32)
		target_link_libraries(pz odbc32.lib ws2_32.lib mswsock.lib)
  endif(WIN32)
endif(USING_LOG4CXX)

if (USING_METIS)
  target_link_libraries(pz ${METIS_LIB})
endif(USING_METIS)

if(USING_OPENSSL)
	target_link_libraries(crypto)
endif(USING_OPENSSL)

if(USING_TBB)
	target_link_libraries(pz ${TBB_LIB})
endif(USING_TBB)

if (USING_OPENMP)
	target_link_libraries(pz ${OPENMP_LIB})
endif(USING_OPENMP)

if (USING_LIKWID)
	target_link_libraries(pz ${LIKWID_LIB})
endif(USING_LIKWID)

if (USING_LIBNUMA)
	target_link_libraries(pz ${LIBNUMA_LIB})
endif(USING_LIBNUMA)

if(USING_BLAS)
  target_link_libraries(pz ${CBLAS_LIB} ${ATLAS_LIB})
endif(USING_BLAS)

if(USING_PAPI)
  target_link_libraries(pz ${PAPI_LIB})
endif(USING_PAPI)

if(USING_HWLOC)
	target_link_libraries(pz ${HWLOC_LIB})
endif(USING_HWLOC)

if(USING_MKL)
	target_link_libraries(pz ${MKL_LIB_INTEL} ${MKL_LIB_CORE} ${MKL_LIB_THREAD} ${COMPOSER_OMP})
endif(USING_MKL)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

install(TARGETS pz DESTINATION ${NEOPZ_LIB})

#creating pzmatrix lib -- It is comment now, because it is not needed when the pzlib is installed.
#set(pzmatrixsources 
#	${pzcommon_STAT_SRCS} ${pzutil_STAT_SRCS} ${pzperfutil_STAT_SRCS} ${pzsave_STAT_SRCS} 
#	${pzmatrix_STAT_SRCS} ${pzlinearsolvers_STAT_SRCS}
#	)

#add_library(pzmatrix ${pzmatrixsources})
#install(TARGETS pzmatrix ARCHIVE DESTINATION lib)

#Verifying the content of the strings
#file(APPEND "${LogFile}" "${refpatterns_STAT_SRCS}   ..  \n\n")
