/**
 * @file
 * @brief Contains declaration of TPZHCurlNedFLinEl class which implements a computational HCurl-conforming element for 1D boundaries.
 */

#ifndef TPZHCURLNEDFLINEL_H_
#define TPZHCURLNEDFLINEL_H_

#include "pzintel.h"
#include "pzshapelinear.h"

/** \addtogroup CompElement */
/** @{ */
/**
 * @brief Implements a computational HCurl-conforming element for 1D boundaries. \ref CompElement "Computational Element"
 * @author Francisco Orlandini
 * @since Jun 19, 2017.
 */
 

class TPZHCurlNedFLinEl : public TPZInterpolatedElement{
    
public:
    TPZHCurlNedFLinEl(TPZCompMesh &mesh, TPZGeoEl *reference, long &index);
    
    TPZHCurlNedFLinEl(TPZCompMesh &mesh, const TPZHCurlNedFLinEl &copy);
    
    TPZHCurlNedFLinEl ( TPZCompMesh &mesh,
                       const TPZHCurlNedFLinEl &copy,
                       std::map<long,long> & gl2lcElMap);
    
    TPZHCurlNedFLinEl (TPZCompMesh &mesh,
                       const TPZHCurlNedFLinEl &copy,
                       std::map<long,long> & gl2lcConMap,
                       std::map<long,long> & gl2lcElMap);
    
    TPZHCurlNedFLinEl();
    /** @brief Destructor, does nothing */
    virtual ~TPZHCurlNedFLinEl();
    
    virtual TPZHCurlNedFLinEl *Clone(TPZCompMesh &mesh) const;
    
    virtual TPZHCurlNedFLinEl *ClonePatchEl(TPZCompMesh &mesh, std::map<long,long> & gl2lcConMap, std::map<long,long> & gl2lcElMap) const;
    
    virtual int Dimension() const;
    
    virtual int NCornerConnects() const;
    
    virtual int NConnects() const;
    
    virtual long ConnectIndex(int i) const;
    
    virtual void SetConnectIndex(int i, long connectindex);
    
    virtual int NSideConnects(int side) const;
    
    virtual int SideConnectLocId(int con,int is) const;
    
    virtual int NConnectShapeF(int con, int order) const;
    
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
    
    void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                      TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                      REAL &detjac, TPZFMatrix<REAL> &jacinv,
                      TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &curlPhiHat, TPZFMatrix<REAL> &curlPhi);
    
    void ShapeTransform(const TPZFMatrix<REAL> &phiHat, const TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi);
    
    void CurlTransform(const TPZFMatrix<REAL> &curlPhiHat, const TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &curlPhi);
    
    virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &curlPhiHat);
    
    virtual void SetCreateFunctions(TPZCompMesh *mesh){
        mesh->SetAllCreateFunctionsHCurl();
    }//TODO: is this necessary?
    
protected:
    
    TPZManVector<long,pzshape::TPZShapeLinear::NSides> fConnectIndexes;
    
    int fSideOrient;
    
    pzshape::TPZShapeLinear::IntruleType fIntRule;
};

#endif /* TPZHCURLNEDFLINEL_H_ */
