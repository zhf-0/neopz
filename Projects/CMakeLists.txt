# @file neopz/Projects/CMakeList.txt  -- First CMake file for Projects directory 

#
# COMPLEX NUMBER PROJECTS
#
# If your project use Complex numbers (STATE type = complex any type), please add it on the list below
if(BUILD_COMPLEX_PROJECTS)

	#Example to work with scalar as complex
	add_subdirectory(Helmholtz1D)

	#Examples used as tutorial for FEEC group (Prof Hugo)
	add_subdirectory(PDElinear1D)
	add_subdirectory(PDElinear2D)

else()

	if(BUILD_TUTORIAL)

		#Example to reading and using geometric mesh
		add_subdirectory(testgeom)

		#Testing and using shape functions
		add_subdirectory(testshape)

		#Curso PZ:
		add_subdirectory(CursoPZ)
		#add_subdirectory(PDElinear1D)
		#add_subdirectory(PDElinear2D)

		#Generic Example 2D with elasticity, using matrices phi, phi, etc
		add_subdirectory(testmat)

		#Using visual matrix of the DX
		add_subdirectory(visualmatrix)
		
		#Using discontinuous approximation functions
		add_subdirectory(descontinuo)

		#Doing h-refinement
	#	add_subdirectory(adaptivity)

		#source_group(CursoPZ FILES ${pzcurso})
		#source_group(Shape FILES ${pztutorshape})
		#source_group(Geometry FILES ${pztutorgeom})
		#source_group(Matrix FILES ${pztutormat})
		#source_group(VisualMatrixDX FILES ${visualmatrix})

		#set(pztutorial ${pztutorial} ${pzcurso} ${pztutorshape} ${pztutorgeom} ${pztutormat} ${visualmatrix})

		#add_custom_target(PZTutorial SOURCES ${pztutorial})


	endif(BUILD_TUTORIAL)

	if(BUILD_PROJECTS)
		
		add_subdirectory(${PROJECT_SOURCE_DIR}/../ProdIndex-Source ${PROJECT_SOURCE_DIR}/../ProdIndex-Source)
		
	endif(BUILD_PROJECTS)


endif(BUILD_COMPLEX_PROJECTS)