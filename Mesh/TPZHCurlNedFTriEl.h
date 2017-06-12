/**
 * @file
 * @brief Contains declaration of TPZHCurlNedFTriEl class which implements
 * on triangular topology the first family of computational Hcurl-conforming
 * elements proposed by J.C. Nédélec.
 */

#ifndef TPZHCURLNEDFTRIEL_H
#define TPZHCURLNEDFTRIEL_H

#include "pzintel.h"
#include "pzshapetriang.h"

/**
 * @brief Implements the Nédélec HCurl-conforming element of the first kind,
 * as in Nedelec, J.C. Numer. Math. (1980) 35: 315. doi:10.1007/BF01396415,
 * on a triangular topology. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
class TPZHCurlNedFTriEl : public TPZInterpolatedElement {
	
public:
	/**
	 * @brief Constructor with a mesh and geometric element as arguments
	 * @param mesh mesh object into which the element will insert itself
	 * @param reference geometric element to which this element will refer
	 * @param index index in the vector of elements of mesh where this element was inserted
	 */
	TPZHCurlNedFTriEl(TPZCompMesh &mesh, TPZGeoEl *reference, long &index);
	
	/**
	 * @brief Constructor aimed at creating a copy of an interpolated element within a new mesh
	 */
	TPZHCurlNedFTriEl(TPZCompMesh &mesh, const TPZHCurlNedFTriEl &copy);
	
	/**
	 * @brief Copy the given element into a new patch mesh
	 * @param mesh patch mesh
	 * @param copy element to be copied
	 * @param gl2lcElMap map the indexes of the orginal mesh to the patched mesh
	 */
	TPZHCurlNedFTriEl ( TPZCompMesh &mesh,
							const TPZHCurlNedFTriEl &copy,
							std::map<long,long> & gl2lcElMap);
	
    /**
     * @brief Constructor used to generate patch mesh... generates a map of connect index from
     * global mesh to clone mesh
     */
    TPZHCurlNedFTriEl (TPZCompMesh &mesh,
                  const TPZHCurlNedFTriEl &copy,
                  std::map<long,long> & gl2lcConMap,
                  std::map<long,long> & gl2lcElMap);
    
	TPZHCurlNedFTriEl();
	/** @brief Destructor, does nothing */
	virtual ~TPZHCurlNedFTriEl();
    
    virtual TPZHCurlNedFTriEl *Clone(TPZCompMesh &mesh) const;
    
    virtual TPZHCurlNedFTriEl *ClonePatchEl(TPZCompMesh &mesh, std::map<long,long> & gl2lcConMap, std::map<long,long> & gl2lcElMap) const;
    
    virtual int Dimension() const;
    
    virtual int NCornerConnects() const;
    
    virtual int NConnects() const;
    
    virtual long ConnectIndex(int i) const;
    
    virtual void SetConnectIndex(int i, long connectindex);
    
    virtual int NSideConnects(int side) const;
    
    virtual int SideConnectLocId(int con,int is) const;
    
    virtual int NConnectShapeF(int con) const;
    
    virtual int NConnectShapeF(int con, int order) const;
    
    void CreateDofVec();
    
    void EvaluateShapeF(const TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &curlPhiHat);
    
    virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &curlPhiHat);
    
    virtual void SideShapeFunction(int side, TPZVec<REAL> &point, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &curlPhiHat);
    
    virtual void SetIntegrationRule(int ord);
    
    virtual const TPZIntPoints &GetIntegrationRule() const;
    
    virtual TPZIntPoints &GetIntegrationRule();
    
    virtual void SetPreferredOrder(int order);
    
    virtual void GetInterpolationOrder(TPZVec<int> &ord);
    
    virtual int PreferredSideOrder(int side);
    
    int ConnectOrder(int connect) const;
    
    virtual int EffectiveSideOrder(int side) const;
    
    virtual void SetSideOrder(int side, int order);
    
    virtual TPZTransform TransformSideToElement(int side);
    
protected:
	
	static bool fHaveShapeFBeenCreated;
	bool fHaveDofVecBeenCreated;
    TPZManVector< TPZManVector<TPZManVector<REAL,31>,255 >, 4> *fPhiVecDofs;
    TPZManVector<long,pzshape::TPZShapeTriang::NSides> fConnectIndexes;
    
    TPZManVector<int, pzshape::TPZShapeTriang::NFaces> fSideOrient;//TODO: TRANSFER TO LINEAR EL
    
    pzshape::TPZShapeTriang::IntruleType fIntRule;
	
};

TPZCompEl * CreateHCurlNedFLinEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
TPZCompEl * CreateHCurlNedFTriEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
#endif
