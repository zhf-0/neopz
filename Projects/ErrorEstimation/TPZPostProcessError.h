//
//  TPZPostProcessError.hpp
//  PZ
//
//  Created by Philippe Devloo on 6/30/16.
//
//

#ifndef TPZPostProcessError_hpp
#define TPZPostProcessError_hpp

#include <stdio.h>
#include <iterator>

#include "pzmanvector.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzfunction.h"

struct TPZPatch
{
    // connect index of the partition of unity mesh
    long fPartitionConnectIndex;
    // vector of element indices of HDiv elements
    TPZManVector<long,20> fElIndices;
    // vector of open set of connect indices that will be used for flux and pressure computations
    TPZManVector<long,25> fConnectIndices;
    
    // vector of closed set of connect indexes included in the elements
    TPZManVector<long,30> fAllConnectIndices;
    
    void ClosedSet(std::set<long> &closed)
    {
//        std::copy (bar.begin(),bar.end(),std::inserter(foo,it));
        std::copy(&(fConnectIndices[0]),(&(fConnectIndices[0])+fConnectIndices.size()),std::inserter(closed,closed.begin()));
    }
    
    TPZPatch() : fPartitionConnectIndex(-1)
    {
        
    }
    
    TPZPatch(const TPZPatch &copy) : fPartitionConnectIndex(copy.fPartitionConnectIndex), fElIndices(copy.fElIndices),
    fConnectIndices(copy.fConnectIndices), fAllConnectIndices(copy.fAllConnectIndices)
    {
        
    }
    TPZPatch &operator=(const TPZPatch &copy)
    {
        fPartitionConnectIndex = copy.fPartitionConnectIndex;
        fElIndices = copy.fElIndices;
        fConnectIndices = copy.fConnectIndices;
        fAllConnectIndices = copy.fAllConnectIndices;
        return *this;
    }
    
    void Print(std::ostream &out)
    {
        out << "The generating partitionindex = " << fPartitionConnectIndex << std::endl;
        out << "Element indices " << fElIndices << std::endl;
        out << "Open set connect indices " << fConnectIndices << std::endl;
        out << "Closed set connect indices " << fAllConnectIndices << std::endl;
    }
};

class TPZPostProcessError
{
public:
    TPZPostProcessError(TPZVec<TPZCompMesh *> &meshvec);
    
private:
    // mesh vector
    TPZManVector<TPZCompMesh *,6> fMeshVector;
    
    // vector of vector of patches
    // each vector of patches corresponds to one color
    TPZManVector<TPZStack<TPZPatch>, 10> fVecVecPatches;
    
    // build vector of patches of a same color
    void BuildPatchStructures();
    
    // print the relevant information of the patches
    void PrintPatchInformation(std::ostream &out);
    
    // original connect sequence numbers
    TPZVec<long> fConnectSeqNumbers;
    
    // multiplying coefficients of the reconstructed fluxes and pressures
    TPZFMatrix<STATE> fSolution;
    
    // block corresponding to the original connect sequence numbers
    TPZBlock<STATE> fBlock;
    
    // plot the reconstructed fluxes
    void PlotFluxes(const std::string &filename);
    
    // solve for the reconstructed fluxes of a given color. Add the flux coefficients
    void ComputePatchFluxes();
    
public:
    
    // Collect the connect indices and elements which will contribute to the patch caracterized by the set of nodes
    // generally each node will form a patch
    TPZPatch BuildPatch(TPZCompElSide &seed);

    // compute the estimated H1 seminorm errors
    void ComputeHDivSolution();
    
    // compute the estimated H1 seminorm errors
    void ComputeElementErrors(TPZVec<STATE> &elementerrors);
    
    // compute the exact element errors
    void ComputeExactH1SemiNormErrors(TPZFunction<STATE> &exact, TPZVec<STATE> &exacterror)
    {
        DebugStop();
    }
    
};

#endif /* TPZPostProcessError_hpp */
