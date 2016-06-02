//
//  TPZAssembleCS.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#ifndef TPZAssembleCS_hpp
#define TPZAssembleCS_hpp

#include <stdio.h>
#include "TPZAssemble.h"

/// class that performs parallel assembly using thread workers and a mutex for accessing the global matrix
// there is no thread that does the assembly
// this class has a particular implementation when linking with TBB
class TPZAssembleCS : public TPZAssemble
{
    
public:
    
    TPZAssembleCS(TPZStructMatrix::TPZAssembleConfig &config) : TPZAssemble(config)
    {
        
    }
    
    TPZAssembleCS(const TPZAssembleCS &copy) : TPZAssemble(copy)
    {
        
    }
    
    virtual ~TPZAssembleCS()
    {
        
    }
    
    virtual TPZAssemble *Clone(TPZStructMatrix::TPZAssembleConfig &config) const
    {
        TPZAssembleCS *result = new TPZAssembleCS(*this);
        result->fConfig = &config;
        return result;
    }
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs);
    

    
    /** @brief Structure to manipulate thread to solve system equations */
    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZAssembleCS *strmat,TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZAssembleCS *strmat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds);
        /** @brief Destructor: Destroy the mutex semaphores and others */
        ~ThreadData();
        /** @brief Look for an element index which needs to be computed and put it on the stack */
        long NextElement();
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief Establish whether the element should be computed */
        bool ShouldCompute(int matid)
        {
            return fAssemble->ShouldCompute(matid);
        }
        
        /** @brief Current structmatrix object */
        TPZAssembleCS *fAssemble;
        /** @brief Gui interface object */
        TPZAutoPointer<TPZGuiInterface> fGuiInterface;
        /** @brief Global matrix */
        TPZMatrix<STATE> *fGlobMatrix;
        /** @brief Global rhs vector */
        TPZFMatrix<STATE> *fGlobRhs;
        /** @brief List of computed element matrices (autopointers?) */
        std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > > fSubmitted;
        /** @brief  Current element */
        long fNextElement;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t fAccessElement;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t fAccessElementK;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t fAccessElementF;
        /** @brief Semaphore (to wake up assembly thread) */
        TPZSemaphore fAssembly;
    };
    
    
};


#endif /* TPZAssembleCS_hpp */
