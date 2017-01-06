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
    
    /// Section of the phi vector associated with this volume element
    TPZFMatrix<STATE> fPhi;
    
    /// Eigenvlues associated with the internal shape functions
    TPZManVector<STATE> fEigenvalues;
    
    /// Multiplier coeficients associated with the solution
    TPZFMatrix<STATE> fCoeficients;
    
    /// vector of local indices of multipliers in the group
    TPZManVector<long> fLocalIndices;
    
    /// extend the border shape functions for SBFem computations
    void ExtendShapeFunctions(TPZMaterialData &data1d, TPZMaterialData &data2d);
    
    /// Density associated with the mass matrix
    REAL fDensity;

public:
    
    TPZSBFemVolume(TPZCompMesh &mesh, TPZGeoEl *gel, long &index);
    
    /// Compute the E0, E1 and E2 matrices
    void ComputeKMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0);
    
    /// Data structure initialization
    void SetSkeleton(long skeleton)
    {
        fSkeleton = skeleton;
    }
    
    void SetElementGroupIndex(long index)
    {
        fElementGroupIndex = index;
        std::map<long,int> globtolocal;
        TPZCompEl *cel = Mesh()->Element(index);
        int nc = cel->NConnects();
        TPZManVector<int,10> firsteq(nc+1,0);
        for (int ic = 0; ic<nc; ic++) {
            globtolocal[cel->ConnectIndex(ic)] = ic;
            TPZConnect &c = cel->Connect(ic);
            firsteq[ic+1] = firsteq[ic]+c.NShape()*c.NState();
        }
        int neq = 0;
        nc = NConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = Connect(ic);
            neq += c.NShape()*c.NState();
        }
        fLocalIndices.Resize(neq);
        int count = 0;
        for (int ic=0; ic<nc; ic++) {
            long cindex = ConnectIndex(ic);
#ifdef PZDEBUG
            if (globtolocal.find(cindex) == globtolocal.end()) {
                DebugStop();
            }
#endif
            TPZConnect &c = Connect(ic);
            int neq = c.NShape()*c.NState();
            int locfirst = firsteq[globtolocal[cindex]];
            for (int eq = 0; eq<neq; eq++) {
                fLocalIndices[count++] = locfirst+eq;
            }
        }
#ifdef PZDEBUG
        if(count != neq) DebugStop();
#endif
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
    
    /** @brief return the density associated with the element */
    REAL Density()
    {
        return fDensity;
    }
    
    /** @brief assign a different density */
    void SetDensity(REAL density)
    {
#ifdef PZDEBUG
        if(density <= 0.) DebugStop();
#endif
        fDensity = density;
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
    

    /// initialize the data structures of the eigenvectors and eigenvalues associated with this volume element
    void SetPhiEigVal(TPZFMatrix<STATE> &phi, TPZManVector<STATE> &eigval);
    
    /** @brief Loads the solution within the internal data structure of the element */
    /**
     * Is used to initialize the solution of connect objects with dependency. \n
     * Is also used to load the solution within SuperElements
     */
    virtual void LoadCoef(TPZFMatrix<STATE> &coef);
    
    /**
     * @brief Computes solution and its derivatives in the local coordinate qsi.
     * @param qsi master element coordinate
     * @param sol finite element solution
     * @param dsol solution derivatives
     * @param axes axes associated with the derivative of the solution
     */
    virtual void ComputeSolution(TPZVec<REAL> &qsi,
                                 TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
    

    /**
     * @brief Calculates the solution - sol - for the variable var
     * at point qsi, where qsi is expressed in terms of the
     * master element coordinates
     * @param qsi master element coordinate
     * @param var variable name
     * @param sol vetor for the solution
     */
    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol);
    

    /**
     * @brief Prints element data
     * @param out Indicates the device where the data will be printed
     */
    virtual void Print(std::ostream &out = std::cout) const
    {
        out << "Printing " << __PRETTY_FUNCTION__ << std::endl;
        TPZCompEl::Print(out);
        out << "Group Element Index " << fElementGroupIndex << std::endl;
        out << "Skeleton Element Index " << fSkeleton << std::endl;
        out << "Local Indices " << fLocalIndices << std::endl;
        out << "Coeficients ";
        fCoeficients.Print("Coef =",out,EMathematicaInput);
        out << "Displacement Eigenvectors\n";
        fPhi.Print("Phi = ",out,EMathematicaInput);
        
    }
    
    void CreateGraphicalElement(TPZGraphMesh &, int);

};


TPZCompEl * CreateSBFemCompEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);



#endif /* TPZSBFemVolume_hpp */
