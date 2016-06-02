//
//  TPZAssemble.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#ifndef TPZAssemble_hpp
#define TPZAssemble_hpp

#include <stdio.h>
#include "pzstrmatrix.h"

class TPZAssemble
{
    
protected:
    
    TPZStructMatrix::TPZAssembleConfig *fConfig;
    
public:
    
    TPZAssemble(TPZStructMatrix::TPZAssembleConfig &config) : fConfig(0)
    {
        fConfig = &config;
    }
    
    virtual ~TPZAssemble()
    {
        
    }
    
    TPZStructMatrix::TPZAssembleConfig *Config()
    {
        return fConfig;
    }
    
    virtual TPZAssemble *Clone(TPZStructMatrix::TPZAssembleConfig &config) const = 0;
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs);
    
    /** @brief Assemble the global right hand side */
    virtual void Assemble(TPZFMatrix<STATE> & rhs);
    
    /** @brief Find the order to assemble the elements */
    static void OrderElement(TPZCompMesh *cmesh, TPZVec<long> &ElementOrder);
    
    /** @brief Create blocks of elements to parallel processing */
    static void ElementColoring(TPZCompMesh *cmesh, TPZVec<long> &elSequence, TPZVec<long> &elSequenceColor, TPZVec<long> &elBlocked);
    

protected:
    
    /** @brief Establish whether the element should be computed */
    bool ShouldCompute(int matid) const
    {
        const unsigned int size = fConfig->fMaterialIds.size();
        return size == 0 || fConfig->fMaterialIds.find(matid) != fConfig->fMaterialIds.end();
    }

    /// filter out the equations which are out of the range
    void FilterEquations(TPZVec<long> &origindex, TPZVec<long> &destindex) const;

    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Serial_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs);
    
    /** @brief Assemble the global right hand side */
    virtual void Serial_Assemble(TPZFMatrix<STATE> & rhs);
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs) = 0;
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs) = 0;
    

};

#endif /* TPZAssemble_hpp */
