//
//  TPZAssembleCS.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#ifndef TPZAssembleCSTBB_hpp
#define TPZAssembleCSTBB_hpp

#include <stdio.h>
#include "TPZAssembleCS.h"

/// class that performs parallel assembly using thread workers and a mutex for accessing the global matrix
// there is no thread that does the assembly
// this class has a particular implementation when linking with TBB
class TPZAssembleCSTBB : public TPZAssembleCS
{
    
public:
    
    TPZAssembleCSTBB(TPZStructMatrix::TPZAssembleConfig &config) : TPZAssembleCS(config)
    {
        
    }
    
    TPZAssembleCSTBB(const TPZAssembleCSTBB &copy) : TPZAssembleCS(copy)
    {
        
    }
    
    virtual ~TPZAssembleCSTBB()
    {
        
    }
    
    virtual TPZAssemble *Clone(TPZStructMatrix::TPZAssembleConfig &config) const
    {
        TPZAssembleCSTBB *result = new TPZAssembleCSTBB(*this);
        result->fConfig = &config;
        return result;
    }
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs);
    

#ifdef USING_TBB
    struct AssembleTask {
        ThreadData *data;
        AssembleTask(ThreadData *dt) : data(dt) {};
        void operator()(const tbb::blocked_range<size_t>& range) const;
    };
#endif
    

    
};


#endif /* TPZAssembleCS_hpp */
