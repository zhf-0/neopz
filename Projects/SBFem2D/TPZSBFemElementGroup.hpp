//
//  TPZSBFemElementGroup.hpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#ifndef TPZSBFemElementGroup_hpp
#define TPZSBFemElementGroup_hpp

#include <stdio.h>

#include "pzelementgroup.h"


class TPZSBFemElementGroup : public TPZElementGroup
{
    
    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFMatrix<STATE> fPhi;
    
    /// Inverse of the eigenvector matrix (transfers eigenvector coeficients to side shape coeficients)
    TPZFMatrix<STATE> fPhiInverse;
    
    /// Vector of eigenvalues of the SBFem analyis
    TPZManVector<STATE> fEigenvalues;
    
public:
    
    TPZSBFemElementGroup() : TPZElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemElementGroup(TPZCompMesh &mesh, long &index) : TPZElementGroup(mesh,index)
    {
        
    }
    
    /// Compute the SBFem matrices
    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2);
    
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
    
    
    /**
     * @brief Computes the element right hand side
     * @param ef element load vector(s)
     */
    virtual void CalcResidual(TPZElementMatrix &ef)
    {
        TPZElementMatrix ek(Mesh(),TPZElementMatrix::EK);
        CalcStiff(ek,ef);
    }
    

    /// method to assemble E0, E1, E2

    /// method to compute the stiffness
    /// method to compute the solution
};

#endif /* TPZSBFemElementGroup_hpp */
