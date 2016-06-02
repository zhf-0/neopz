//
//  TPZAssemble.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#ifndef TPZAssembleGC_hpp
#define TPZAssembleGC_hpp

#include <stdio.h>
#include "pzstrmatrix.h"
#include "TPZAssemble.h"

/// order the elements acording to the equation numbering and then color the elements
// this scheme should allow to assemble the elements and/or right hand side without ever blocking
class TPZAssembleGlobalColor : public TPZAssemble
{
    
    
public:
    
    TPZAssembleGlobalColor(TPZStructMatrix::TPZAssembleConfig &config) : TPZAssemble(config)
    {
    }
    
    virtual ~TPZAssembleGlobalColor()
    {
        
    }
    
    
    virtual TPZAssemble *Clone(TPZStructMatrix::TPZAssembleConfig &config) const
    {
        TPZAssembleGlobalColor *result = new TPZAssembleGlobalColor(*this);
        result->fConfig = &config;
        return result;
    }
    
protected:
    
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs);
    
    /** @brief Structure to manipulate thread to solve system equations */
    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZAssembleGlobalColor *strmat,TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZAssembleGlobalColor *strmat, TPZFMatrix<STATE> &rhs);
        /** @brief Destructor: Destroy the mutex semaphores and others */
        ~ThreadData();
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief The function which will compute the assembly */
        bool ShouldCompute(int matid)
        {
            return fStruct->ShouldCompute(matid);
        }
        
        /** @brief Current structmatrix object */
        TPZAssembleGlobalColor *fStruct;
        /** @brief Global matrix */
        TPZMatrix<STATE> *fGlobMatrix;
        /** @brief Global rhs vector */
        TPZFMatrix<STATE> *fGlobRhs;
        /** @brief List of computed element matrices (autopointers?) */
        std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted;
        /** @brief Elements which are being processed */
        std::set<int> fProcessed;
        /** @brief  Current element */
        long fNextElement;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t fAccessElement;
        /** @brief Semaphore (to wake up assembly thread) */
        TPZSemaphore fAssembly;
        
        pthread_cond_t fCondition;
        bool fSleeping;
        
        // Vectors for mesh coloring
        std::map<int,int> felBlocked;
        /// Vector for mesh coloring
        TPZVec<long> *fnextBlocked, *felSequenceColor;
        
        static void *ThreadWorkResidual(void *datavoid);
    };
    
    friend struct ThreadData;
    
    /** @brief Vectors for mesh coloring */
    /// fnextBlocked : the lowest element index that the element i will block
    /// felSequenceColor : element sequence in which the element should be processed
    TPZVec<long> fnextBlocked, felSequenceColor;
   

};

#endif /* TPZAssemble_hpp */
