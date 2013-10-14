/*
 *  TPZMultiphyisicsCompElMHM.h
 *  PZ
 *
 *  Created by Frederico on 01/02/2013.
 *  Copyright 2013 LNCC. All rights reserved.
 *
 */

#ifndef PZMULTIPHYSICCOMPELMHMH
#define PZMULTIPHYSICCOMPELMHMH

#include <iostream>

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzmaterialdata.h"

#include "pzelctemp.h"

class TPZTransform;

/** @brief class to create a compute element multiphysics */
template <class TSHAPE, class TGeometry>
class TPZMultiphysicsCompElMHM : public TPZMultiphysicsCompEl<TGeometry>
{

public:
	/**
	 * @brief Creates a multiphysic computational element within mesh. 
	 * @param mesh mesh multiphysic where will be created the element
	 * @param gel geometric element for which the computational element will be created
	 * @param index new elemen index
	 */
	TPZMultiphysicsCompElMHM(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

	TPZMultiphysicsCompElMHM();

	virtual ~TPZMultiphysicsCompElMHM();

	/**
	 * @brief Post processing method which computes the solution for the var post processed variable.
	 * @param qsi coordinate of the point in master element space where the solution will be evaluated
	 * @param var variable which will be computed
	 * @param sol (output) solution computed at the given point
	 * @see TPZMaterial::VariableIndex
	 * @see TPZMaterial::NSolutionVariables
	 * @see TPZMaterial::Solution
	 */
	/** The var index is obtained by calling the TPZMaterial::VariableIndex method with a post processing name */
	virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol);

	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix &axes);

	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
	 * This method will function for both volumetric and interface elements
	 * @param qsi master element coordinate of the interface element
	 * @param normal vector
	 * @param leftsol finite element solution
	 * @param dleftsol solution derivatives
	 * @param leftaxes axes associated with the left solution
	 * @param rightsol finite element solution
	 * @param drightsol solution derivatives
	 * @param rightaxes axes associated with the right solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZVec<REAL> &normal,
                    TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix &leftaxes,
					TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix &rightaxes);
	
	/**
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
	 * @param axes [in] axes indicating the direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
								 const TPZFMatrix &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
	
	/**
	* @brief Computes the element stiffness matrix and right hand side
	* @param ek element matrix
	* @param ef element right hand side
	*/
	virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

	/** @brief Initialize element matrix in which is computed CalcStiff */
	void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);

	/** 
	 * @brief Initialize a material data vector and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	void InitMaterialData(TPZVec<TPZMaterialData > &dataVec);

    virtual void EvaluateError(void (* /*fp*/)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv), TPZVec<REAL> &/*errors*/,TPZBlock * /*flux*/);

};

/** @brief Creates computational point element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsMHMPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational linear element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsMHMLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational quadrilateral element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsMHMQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational triangular element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsMHMTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational cube element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsMHMCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational prismal element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsMHMPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational pyramidal element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsMHMPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @brief Creates computational tetrahedral element for Multiphysics approximate space */
TPZCompEl *CreateMultiphysicsMHMTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif
