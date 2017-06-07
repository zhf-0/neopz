/**
 * @file
 * @brief Contains declaration of TPZNedFirstKindEl class which implements
 * on triangular topology the first family of computational Hcurl-conforming
 * elements proposed by J.C. Nédélec.
 */

#ifndef PZNEDFIRSTKINDEL_H
#define PZNEDFIRSTKINDEL_H

#include "pzinterpolationspace.h"
class TPZMaterialData;

/**
 * @brief Implements the Nédélec HCurl-conforming element of the first kind,
 * as in Nedelec, J.C. Numer. Math. (1980) 35: 315. doi:10.1007/BF01396415,
 * on a triangular topology.
 * \ref CompElement "Computational element"
 * @since June 02, 2017
 * @ingroup CompElement
 */
class TPZNedFirstKindEl : public TPZInterpolationSpace
{
public:
	
	/** @brief Default constructor */
	TPZNedFirstKindEl();
	
	/** @brief Default destructor */
	virtual ~TPZNedFirstKindEl();
	/** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] intpoint point in master element coordinates 
	 * @param[in] data stores all input data
	 */
	void ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data);
	/**
	 * @name Computational methods
	 * @brief Methods used to perform computations on the interpolated element
	 * @{
	 */
    
	/**
	 <#description#>

	 @param intpoint <#intpoint description#>
	 @param X <#X description#>
	 @param jacobian <#jacobian description#>
	 @param axes <#axes description#>
	 @param detjac <#detjac description#>
	 @param jacinv <#jacinv description#>
	 @param phi <#phi description#>
	 @param dphi <#dphi description#>
	 @param dphidx <#dphidx description#>
	 */
	void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx);

	/** 
	 * @brief Computes the shape function set at the point x. 
	 * @param qsi point in master element coordinates
	 * @param phi vector of values of shapefunctions, dimension (numshape,1)
	 * @param dphi matrix of derivatives of shapefunctions in master element coordinates, dimension (dim,numshape)
	 */
	/**
	 * This method uses the order of interpolation
	 * of the element along the sides to compute the number of shapefunctions
	 */
	void Shape(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

	void ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi);

	void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data);
	
};

#endif
