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
    
    /// Multiplying coefficients of each eigenvector
    TPZFMatrix<STATE> fCoef;
    
    TPZFMatrix<STATE> fMassMatrix;
    
    /// Compute the mass matrix based on the value of M0 and the eigenvectors
    void ComputeMassMatrix(TPZElementMatrix &M0);
    
public:
    
    TPZSBFemElementGroup() : TPZElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemElementGroup(TPZCompMesh &mesh, long &index) : TPZElementGroup(mesh,index)
    {
        
    }
    
    /// Compute the SBFem matrices
    /// method to assemble E0, E1, E2
    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0);
    
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
    
    /**
     * @brief Prints element data
     * @param out Indicates the device where the data will be printed
     */
    virtual void Print(std::ostream &out = std::cout) const
    {
        out << __PRETTY_FUNCTION__ << std::endl;
        TPZCompEl::Print(out);
        out << "EigenVectors for displacement\n";
        fPhi.Print("Phi = ",out,EMathematicaInput);
        out << "Inverse EigenVectors\n";
        fPhiInverse.Print("PhiInv = ",out,EMathematicaInput);
        out << "EigenValues " << fEigenvalues << std::endl;
        out << "Mass Matrix\n";
        fMassMatrix.Print("Mass = ",out);
        out << "Solution Coeficients\n";
        fCoef.Print("Coef ",out);
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            fElGroup[el]->Print(out);
        }
        out << "End of " << __PRETTY_FUNCTION__ << std::endl;
    }
    

    
    /**
     * @brief Computes the element right hand side
     * @param ef element load vector(s)
     */
    virtual void CalcResidual(TPZElementMatrix &ef)
    {
        TPZElementMatrix ek(Mesh(),TPZElementMatrix::EK);
        CalcStiff(ek,ef);
    }
    

    /** @brief Loads the solution within the internal data structure of the element */
    /**
     * Is used to initialize the solution of connect objects with dependency. \n
     * Is also used to load the solution within SuperElements
     */
    virtual void LoadSolution();


    /// method to compute the stiffness
    /// method to compute the solution
};

#endif /* TPZSBFemElementGroup_hpp */
