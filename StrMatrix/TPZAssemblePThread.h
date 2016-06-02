//
//  TPZAssemblePThread.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#ifndef TPZAssemblePThread_hpp
#define TPZAssemblePThread_hpp

#include <stdio.h>
#include "TPZAssemble.h"

/// class that performs parallel assembly using thread workers and a single thread assembler
class TPZAssemblePThread : public TPZAssemble
{
    
public:
    
    TPZAssemblePThread(TPZStructMatrix::TPZAssembleConfig &config) : TPZAssemble(config)
    {
        
    }
    
    TPZAssemblePThread(const TPZAssemblePThread &copy) : TPZAssemble(copy)
    {
        
    }
    
    virtual ~TPZAssemblePThread()
    {
        
    }

    virtual TPZAssemble *Clone(TPZStructMatrix::TPZAssembleConfig &config) const
    {
        TPZAssemblePThread *result = new TPZAssemblePThread(*this);
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
        ThreadData(TPZAssemblePThread *strmat,TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZAssemblePThread *strmat, TPZFMatrix<STATE> &rhs, std::set<int> &MaterialIds);
        /** @brief Destructor: Destroy the mutex semaphores and others */
        ~ThreadData();
        /** @brief Look for an element index which needs to be computed and put it on the stack */
        long NextElement();
        /** @brief Put the computed element matrices in the map */
        void ComputedElementMatrix(long iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef);
        /** @brief The function which will compute the matrices */
        static void *ThreadWork(void *threaddata);
        /** @brief The function which will compute the assembly */
        static void *ThreadAssembly(void *threaddata);
        /** @brief Establish whether the element should be computed */
        bool ShouldCompute(int matid)
        {
            return fAssemble->ShouldCompute(matid);
        }
        
        /** @brief Current structmatrix object */
        TPZAssemblePThread *fAssemble;
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
    };

};
#endif /* TPZAssemblePThread_hpp */
