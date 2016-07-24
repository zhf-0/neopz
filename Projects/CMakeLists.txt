# @file neopz/Projects/CMakeList.txt  -- First CMake file for Projects directory 


	if(BUILD_PROJECTS)
		
		#
		# COMPLEX NUMBER PROJECTS
		#
		# If your project works with Complex numbers (STATE type = complex any type), please add it on the list below
		add_subdirectory(Elasticity2DTensor)
		add_subdirectory(SBFem2D)
		
	endif(BUILD_PROJECTS)


   
