/**
 * @file
 * @brief Contains the TPZVTKGeoMesh class which implements the graphical mesh to VTK environment to geometric mesh.
 */
/*
 *  TPZVTKGeoMesh.h
 *  Crack
 *
 *  Created by Cesar Lucci on 16/08/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 */

#ifndef TPZVTKGEOMESHH
#define TPZVTKGEOMESHH 

#include <set>
#include "pzgeoel.h"

/**
 * @ingroup post
 * @author Cesar Lucci
 * @since 16/08/10
 * @brief To export a graphical mesh to VTK environment to geometric mesh. \ref post "Post processing"
 */
class TPZVTKGeoMesh
{
	
public:
	
	TPZVTKGeoMesh();
	~TPZVTKGeoMesh();
	
	/** @brief Generate an output of all geomesh to VTK */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, bool matColor = false);
	
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, TPZVec<int> &elData);
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, char *filename, TPZChunkVector<int> &elData);
	
	/**
	 * @brief Based on a given geomesh, just the elements that have an neighbour with a given material id will be exported to an VTK file
	 */
	static void PrintGMeshVTKneighbour_material(TPZGeoMesh *gmesh, std::ofstream &file, int neighMaterial, bool matColor = false);
	
	/**
	 * @brief Based on a given geomesh, just the elements that have the given material id will be exported to an VTK file
	 */
	static void PrintGMeshVTKmy_material(TPZGeoMesh *gmesh, std::ofstream &file, std::set<int> myMaterial, bool matColor = false);
	
	static int GetVTK_ElType(TPZGeoEl *gel);
};

#endif
