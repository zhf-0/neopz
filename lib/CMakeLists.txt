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
source_group("Random" FILES ${pzrandom_STAT_SRCS})
source_group("Optimization" FILES ${pzoptimization_STAT_SRCS})

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
	${pzrandom_STAT_SRCS} 
	${pzoptimization_STAT_SRCS} 
	)

#creating target pz library and installing it
if (BUILD_PYTHON_BINDING)
	IF (NOT WIN32)
		#Build a shared library to be used with python binding
		add_library(pz SHARED ${pzsources})
	ELSE()
		add_library(pz STATIC ${pzsources})
	ENDIF()
else()
  #Build a static library (default)
  add_library(pz ${pzsources})
endif()

target_link_libraries(pz ${PTHREAD_LIB})

IF (USING_IPO)
  set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
ENDIF()

IF (NOT WIN32)
	target_link_libraries(pz pthread)
ENDIF()

if (USING_BOOST)
   #link_directories(${Boost_LIBRARY_DIRS})
   target_link_libraries(pz ${Boost_LIBRARIES})
   target_include_directories(pz PUBLIC ${Boost_INCLUDE_DIR})
endif()

if (USING_LOG4CXX)
  target_link_libraries(pz ${LOG4CXX_LIB})
  target_include_directories(pz PUBLIC ${LOG4CXX_INCLUDE})
  if(WIN32)
	target_link_libraries(pz odbc32.lib ws2_32.lib mswsock.lib)
  endif()
endif()

if (USING_OPENMP)
  target_link_libraries(pz ${OPENMP_LIB})
  target_include_directories(pz PUBLIC ${OPENMP_INCLUDE})
endif()

if (USING_METIS)
  target_link_libraries(pz ${METIS_LIB})
  target_include_directories(pz PUBLIC ${METIS_INCLUDE})
endif()

if (USING_LAPACK)
  target_link_libraries(pz ${LAPACK_LIB} )
  target_include_directories(pz PUBLIC ${LAPACK_INCLUDE})
endif()

if (USING_BLAS)
  target_link_libraries(pz ${CBLAS_LIB} ) # lapack BLAS
endif()

if (USING_PAPI)
  target_link_libraries(pz ${PAPI_LIB})
endif()

if (USING_MKL)
  target_link_libraries(pz ${MKL_LIB_INTEL} ${MKL_LIB_CORE} ${MKL_LIB_THREAD} ${COMPOSER_OMP})
  target_include_directories(pz PUBLIC ${MKL_INCLUDE})
endif()

if (USING_HWLOC)
  target_link_libraries(pz ${HWLOC_LIB})
endif()


#install(TARGETS pz DESTINATION ${NEOPZ_LIB})


install(TARGETS pz 
		EXPORT PZTargets
		ARCHIVE
		DESTINATION ${NEOPZ_LIB})
	
#creating pzmatrix lib -- It is comment now, because it is not needed when the pzlib is installed.
#set(pzmatrixsources 
#	${pzcommon_STAT_SRCS} ${pzutil_STAT_SRCS} ${pzperfutil_STAT_SRCS} ${pzsave_STAT_SRCS} 
#	${pzmatrix_STAT_SRCS} ${pzlinearsolvers_STAT_SRCS}
#	)

#add_library(pzmatrix ${pzmatrixsources})
#install(TARGETS pzmatrix ARCHIVE DESTINATION lib)

#Verifying the content of the strings
#file(APPEND "${LogFile}" "${refpatterns_STAT_SRCS}   ..  \n\n")
