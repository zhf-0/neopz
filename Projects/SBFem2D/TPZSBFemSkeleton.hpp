//
//  TPZSBFemSkeleton.hpp
//  PZ
//
//  Created by Philippe Devloo on 3/27/16.
//
//

#ifndef TPZSBFemSkeleton_hpp
#define TPZSBFemSkeleton_hpp


#include <stdio.h>
#include "pzcompel.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzelctemp.h"
#include "pzshapelinear.h"

/// class which computes the contribution to K11 K12 and K22 for an SBFem pyramid element
// The derivation of TPZCompEl is to reimplement ComputeSolution for post processing
class TPZSBFemSkeleton : public TPZIntelGen<pzshape::TPZShapeLinear>
{
    
    /// index of the 2D reference element
    long fGeoReference;
    
    
public:
    
    TPZSBFemSkeleton(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) : TPZIntelGen<pzshape::TPZShapeLinear>(mesh,gel,index)
    {
        
    }
        
    
    void SetReferenceGeometric(long refindex)
    {
        fGeoReference = refindex;
    }
    
    /** @brief Method for creating a copy of the element */
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const
    {
        // till I remember how this works
        DebugStop();
    }
    
    /**
     * @brief Method for creating a copy of the element in a patch mesh
     * @param mesh Patch clone mesh
     * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
     * @param gl2lcElMap map the computational elements
     */
    /**
     * Otherwise of the previous clone function, this method don't
     * copy entire mesh. Therefore it needs to map the connect index
     * from the both meshes - original and patch
     */
    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                    std::map<long,long> & gl2lcConMap,
                                    std::map<long,long> & gl2lcElMap) const
    {
        // till I remember how this works
        DebugStop();
        return 0;
    }
    
    virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                         TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                         REAL &detjac, TPZFMatrix<REAL> &jacinv,
                                                 TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidaxes);
    
    void ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data);



};

#endif /* TPZSBFemSkeleton_hpp */
