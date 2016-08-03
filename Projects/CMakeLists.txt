# @file neopz/Projects/CMakeList.txt  -- First CMake file for Projects directory 


	if(BUILD_TUTORIAL)

		if(BUILD_COMPLEX_PROJECTS)
			#add_subdirectory(PDElinear1D)
			#add_subdirectory(PDElinear2D)
		else(BUILD_COMPLEX_PROJECTS)

			# project to save and load a computational mesh
			add_subdirectory(SaveAndLoad)

			#Using discontinuous approximation functions
			add_subdirectory(descontinuo)

                	#Testing gradient reconstruction
			add_subdirectory(Gradient_Reconstruction)
                	#Testing generating ellipse from a vector of points by least square method
			add_subdirectory(EllipseGenerator)
		endif(BUILD_COMPLEX_PROJECTS)
	endif(BUILD_TUTORIAL)
	
	if(BUILD_PROJECTS)
		
		#
		# COMPLEX NUMBER PROJECTS
		#
		# If your project works with Complex numbers (STATE type = complex any type), please add it on the list below
		if(BUILD_COMPLEX_PROJECTS)

			#Example to work with scalar as complex
			#add_subdirectory(Helmholtz1D)
			add_subdirectory(TesteComplexo1D)
			add_subdirectory(StriplineFran)
			add_subdirectory(ValidacaoHCurlFran1)
			add_subdirectory(ValidacaoHCurlFran2)
			
			#add_subdirectory(HDiv_demo)

			#Examples used as tutorial for FEEC group (Prof Hugo)
			#add_subdirectory(PDElinear1D)
			#add_subdirectory(PDElinear2D)

		else(BUILD_COMPLEX_PROJECTS)
			#
			# REGULAR PROJECTS
			#
			# If your project DOES NOT use Complex numbers, please add it on the list below
			#
			add_subdirectory(Elasticity2DTensor)
			add_subdirectory(SBFem2D)
			
			
			########################## Nao aplicavel para Windows
			if (NOT WIN32)
				#Projeto Teo Refinamento Paralelo
				add_subdirectory(RefThread)
				# Projeto Caju: precisamos consertar para rodar no windows! problema com SWX components
				#add_subdirectory(HydraulicFracturePropagation)
			endif()
			##########################

			#projeto Diogo
			add_subdirectory(Plasticity)

			add_subdirectory(Elasticity2D)
			add_subdirectory(Gridgeneration)
			#add_subdirectory(Poisson2DLabyrinth)
			add_subdirectory(Andrade)
			add_subdirectory(MatrixMarket)
			add_subdirectory(Frac1DHdiv)
			add_subdirectory(SubStruct)
			add_subdirectory(MultiphysicsSimulation)
			add_subdirectory(ProblemGirkmann)
			
			add_subdirectory(multplaca)

			#Projeto Lua 2012
			add_subdirectory(Lua01)

			#teste nathan
			add_subdirectory(PlasticityPV)

			#Projeto Denise
			#add_subdirectory(HdivPressureTeste)
			add_subdirectory(HpAdaptivity)
			add_subdirectory(HDiv)
			add_subdirectory(Hdiv3D)
			add_subdirectory(AcoplamentoH1Hdiv)
			add_subdirectory(Validacao)

			#add_subdirectory(prodindex)
			#Projeto Agnaldo
			add_subdirectory(Testes_Agnaldo)
			add_subdirectory(ReducedSpace)
			add_subdirectory(ReducedSpacePlasticity)
			add_subdirectory(HydraulicFracturePropagation)
			add_subdirectory(StabilizedHybridElasticity)
			add_subdirectory(GlobalPostPro)
			# validacaoFranHCurl
			add_subdirectory(ValidacaoHCurlFran1)

			#Projetos de Thiago Dias dos Santos
			#add_subdirectory(SantosProjects)

			add_subdirectory(Bima1)
#			add_subdirectory(EstabSUPG)
#			add_subdirectory(MetodosEstab)
			#
			# Colaboration with Laurent Stainier
			
			add_subdirectory(MatLib)
			add_subdirectory(TrabalhosLNCC15)

			if(USING_MKL)
				add_subdirectory(PardisoSample)
			endif(USING_MKL)
			

		endif(BUILD_COMPLEX_PROJECTS)
		
		
	endif(BUILD_PROJECTS)


   
add_subdirectory(artigoPira)