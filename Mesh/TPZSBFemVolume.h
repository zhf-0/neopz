//
//  TPZSBFemVolume.hpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#ifndef TPZSBFemVolume_hpp
#define TPZSBFemVolume_hpp

#include <stdio.h>
#include "pzcompel.h"
#include "pzelmat.h"

class TPZSBFemVolume : public TPZCompEl
{
    
    /// index of element group
    long fElementGroupIndex;
    
    /// index of the skeleton element
    long fSkeleton;
    
    /// extend the border shape functions for SBFem computations
    void ExtendShapeFunctions(TPZMaterialData &data1d, TPZMaterialData &data2d);

public:
    
    TPZSBFemVolume(TPZCompMesh &mesh, TPZGeoEl *gel, long &index);
    
    /// Compute the E0, E1 and E2 matrices
    void ComputeKMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2);
    
    /// Data structure initialization
    void SetSkeleton(long skeleton)
    {
        fSkeleton = skeleton;
    }
    
    void SetElementGroupIndex(long index)
    {
        fElementGroupIndex = index;
    }
    /** @brief Method for creating a copy of the element */
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const
    {
        // till I remember how this works
        DebugStop();
        return 0;
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

    /** @brief Returns the number of nodes of the element */
    virtual int NConnects() const
    {
        if (fSkeleton == -1) {
            return 0;
        }
        return Mesh()->Element(fSkeleton)->NConnects();
    }
    
    /**
     * @brief Returns the index of the ith connectivity of the element
     * @param i connectivity index who want knows
     */
    virtual long ConnectIndex(int i) const
    {
        if (fSkeleton == -1) {
            DebugStop();
        }
        return Mesh()->Element(fSkeleton)->ConnectIndex(i);
    }
    /** @brief Dimension of the element */
    virtual int Dimension() const
    {
        TPZGeoEl *reference = Reference();
        return reference->Dimension();
    }
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<long> &connectindexes) const
    {
        if (fSkeleton == -1) {
            DebugStop();
        }
        Mesh()->Element(fSkeleton)->BuildCornerConnectList(connectindexes);
    }
    
    /**
     * @brief Set the index i to node inode
     * @param inode node to set index
     * @param index index to be set
     */
    virtual void SetConnectIndex(int inode, long index)
    {
        if (fSkeleton == -1) {
            DebugStop();
        }
        Mesh()->Element(fSkeleton)->SetConnectIndex(inode,index);
    }
    

};


TPZCompEl * CreateSBFemCompEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);



#endif /* TPZSBFemVolume_hpp */
