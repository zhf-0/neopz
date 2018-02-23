/**
 * @file
 * @brief Contains declaration of TPZHCurlNedFLinEl class which implements a
 * computational HCurl-conforming element for 1D boundaries.
 */

#ifndef TPZHCURLNEDFLINEL_H_
#define TPZHCURLNEDFLINEL_H_

#include "pzintel.h"
#include "pzshapelinear.h"
#include "TPZHCurlNedFTriEl.h"
/** \addtogroup CompElement */
/** @{ */
/**
 * @brief Implements a computational HCurl-conforming element for 1D boundaries.
 * \ref CompElement "Computational Element"
 * @author Francisco Orlandini
 * @since Jun 19, 2017.
 */

class TPZHCurlNedFLinEl : public TPZInterpolatedElement {

  public:
	/**
	 * @brief Constructor with a mesh and geometric element as arguments
	 * @param mesh mesh object into which the element will insert itself
	 * @param reference geometric element to which this element will refer
	 * @param index index in the vector of elements of mesh where this element
	 * was inserted
	 */
	TPZHCurlNedFLinEl(TPZCompMesh &mesh, TPZGeoEl *reference, long &index);
	
	/**
	 * @brief Constructor aimed at creating a copy of an interpolated element
	 * within a new mesh
	 */
	TPZHCurlNedFLinEl(TPZCompMesh &mesh, const TPZHCurlNedFLinEl &copy);
	
	/**
	 * @brief Copy the given element into a new patch mesh
	 * @param mesh patch mesh
	 * @param copy element to be copied
	 * @param gl2lcElMap map the indexes of the orginal mesh to the patched mesh
	 */
	TPZHCurlNedFLinEl(TPZCompMesh &mesh, const TPZHCurlNedFLinEl &copy,
					  std::map<long, long> &gl2lcElMap);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of
	 * connect index from global mesh to clone mesh
	 */
	TPZHCurlNedFLinEl(TPZCompMesh &mesh, const TPZHCurlNedFLinEl &copy,
					  std::map<long, long> &gl2lcConMap,
					  std::map<long, long> &gl2lcElMap);
	
	TPZHCurlNedFLinEl();
	/** @brief Destructor, does nothing */
	virtual ~TPZHCurlNedFLinEl();
	
	virtual TPZHCurlNedFLinEl *Clone(TPZCompMesh &mesh) const;
	
	virtual TPZHCurlNedFLinEl *
	ClonePatchEl(TPZCompMesh &mesh, std::map<long, long> &gl2lcConMap,
				 std::map<long, long> &gl2lcElMap) const;
	
	/**
	 Dimension of geometric element
	 
	 @return dimension
	 */
	virtual int Dimension() const;
	
	/**
	 Return the number of connects of
	 the element associated with its vertices.
	 
	 @return number of corner connects
	 */
	virtual int NCornerConnects() const;
	
	/**
	 Returns the number of connects of
	 the element.
	 
	 @return number of connects
	 */
	virtual int NConnects() const;
	
	/**
	 Returns the global index of a connect
	 described by its id (local connectivity index)
	 
	 @param i connect id
	 @return global connect index
	 */
	virtual long ConnectIndex(int i) const;
	
	/**
	 Associates the ith connect of the element
	 to a global index.
	 
	 @param i connect id
	 @param connectindex global connect index
	 */
	virtual void SetConnectIndex(int i, long connectindex);
	
	/**
	 Returns the number of connects of the element
	 associated with a given side.
	 
	 @param side side index
	 @return number of connects associated with side.
	 */
	virtual int NSideConnects(int side) const;
	
	/**
	 Returns the con-ith connect associated with
	 the side is
	 
	 @param con connect number (regarding side is)
	 @param is side index
	 @return local connect id
	 */
	virtual int SideConnectLocId(int con, int is) const;
	
	
	/**
	 Returns the number of shape functions associated with
	 the connect con when it has polynomial order order
	 
	 @param con local connect id
	 @param order polynomial order of the connect
	 @return number of shape functions
	 */
	virtual int NConnectShapeF(int con, int order) const;
	
	
	virtual void SideShapeFunction(int side, TPZVec<REAL> &point,
								   TPZFMatrix<REAL> &phi,
								   TPZFMatrix<REAL> &curlPhiHat);
	
	
	/**
	 Sets integration rule order ord
	 
	 @param ord integration rule order
	 */
	virtual void SetIntegrationRule(int ord);
	
	
	/**
	 Gets current integration rule (points and weights)
	 
	 @return integration rule
	 */
	virtual const TPZIntPoints &GetIntegrationRule() const;
	
	/**
	 Gets current integration rule (points and weights)
	 
	 @return integration rule
	 */
	virtual TPZIntPoints &GetIntegrationRule();
	
	
	/**
	 Sets computational element preferred order.
	 It acts as the maximum order of the element,
	 since it may have a connect with lower order
	 because of a neighbour.
	 
	 @param order preferred approximation order
	 */
	virtual void SetPreferredOrder(int order);
	
	virtual void GetInterpolationOrder(TPZVec<int> &ord);
	
	virtual int PreferredSideOrder(int side);
	
	
	/**
	 Given a local connect id, returns
	 the polynomial order of its highest order
	 function
	 
	 @param connect local connect id
	 @return polynomial order
	 */
	int ConnectOrder(int connect) const;
	
	
	/**
	 Returns effective polynomial order of the
	 functions associated with side side
	 
	 @param side side id
	 @return functions polynomial order
	 */
	virtual int EffectiveSideOrder(int side) const;
	
	/**
	 Sets polynomial order of the functions
	 associated with side side
	 
	 @param side side id
	 @param order polynomial order
	 */
	virtual void SetSideOrder(int side, int order);
	
	virtual TPZTransform<> TransformSideToElement(int side);
	
	
	/**
	 Calculates shape functions and their curls.
	 Given an integration point, it calculates the shape functions
	 (and their curls) in the reference element, the geometric mapping
	 and the Piola covariant transformation for calculating the shape
	 functions along with their curls in the deformed element.
	 
	 @param intpoint integration point @ reference el
	 @param X integration point @ deformed el (DONT use it)
	 @param jacobian jacobian of the geometric mapping (2D->2D)
	 @param axes places deformed element in the cartesian 3D universe
	 @param detjac determinant of the jacobian
	 @param jacinv inverse of the jacobian
	 @param phi shape functions values
	 @param curlPhiHat curl of the shape functions @ reference el
	 @param curlPhi curl of the shape functions @ deformed el
	 */
	void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
					  TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
					  REAL &detjac, TPZFMatrix<REAL> &jacinv,
					  TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &curlPhiHat,
					  TPZFMatrix<REAL> &curlPhi);

	static void CalcShape(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi,
	                                  TPZFMatrix<REAL> &curlPhiHat, TPZVec<int> &order, TPZVec<int> nShapeF);
	
	/**
	 Applies Piola covariant transformation for proper calculation
	 of the shape functions in the deformed element.
	 
	 @param phiHat shape functions @ reference el
	 @param jacinv inverse of the jacobian of the element geometric mapping
	 @param phi shape functions @ deformed el
	 */
	void ShapeTransform(const TPZFMatrix<REAL> &phiHat,
						const TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi);
	
	/**
	 It does nothing. This is a boundary element, the curl is not calculated in it.
	 
	 @param curlPhiHat curl values @ reference el
	 @param jacinv inverse of the jacobian of the element geometric mapping
	 @param curlPhi curl values @ deformed el
	 */
	void CurlTransform(const TPZFMatrix<REAL> &curlPhiHat,
					   const TPZFMatrix<REAL> &jacinv,
					   TPZFMatrix<REAL> &curlPhi);
	
	/**
	 Compute the value of the shape functions
	 evaluated at the integration point.
	 Implementation provided in TPZHCurlNedFLinElShape[,Scaled,N].cpp,
	 along with proper description. See TPZHCurlNedFTriElShape.cpp
	 for details.
	 
	 @param qsi integration point (in reference geometric element)
	 @param phi vector of shape functions values
	 @param curlPhiHat vector of shape functions curl values
	 */
	virtual void Shape(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi,
					   TPZFMatrix<REAL> &curlPhiHat);
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
	
	virtual void SetCreateFunctions(TPZCompMesh *mesh) {
		mesh->SetAllCreateFunctionsHCurl();
        }

      protected:
	//! Stores global indexes of the connects of the element.
	TPZManVector<long, pzshape::TPZShapeLinear::NSides> fConnectIndexes;
	//! stores side orientation (1 or -1)
	/*!
	 It guarantees tangential continuity between elements
	 */
	int	fSideOrient;
	//! stores integration rule points
	pzshape::TPZShapeLinear::IntruleType fIntRule;
};

#endif /* TPZHCURLNEDFLINEL_H_ */
