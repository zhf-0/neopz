//
//  TPZAssemble.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#ifndef TPZAssembleOptimizedColor_hpp
#define TPZAssembleOptimizedColor_hpp

#include <stdio.h>
#include "pzstrmatrix.h"
#include "TPZAssemble.h"

class TPZAssembleOptimizedColor : public TPZAssemble
{
    
protected:
    
    
public:
    
    TPZAssembleOptimizedColor(TPZStructMatrix::TPZAssembleConfig &config) : TPZAssemble(config)
    {
    }
    
    virtual ~TPZAssembleOptimizedColor()
    {
        
    }
    
    TPZAssembleOptimizedColor(const TPZAssembleOptimizedColor &copy) : TPZAssemble(copy)
    {
        DebugStop();
    }
    
    virtual TPZAssemble *Clone(TPZStructMatrix::TPZAssembleConfig &config) const
    {
        TPZAssembleOptimizedColor *result = new TPZAssembleOptimizedColor(*this);
        result->fConfig = &config;
        return result;

    }
    
    // elSequence (input) element sequence acording to the connect sequence numbers
    // elSequenceColor (output) the coloured element sequence
    // elBlocked the element index which needs to have been computed before assembling the element
    // elColors (output) number of elements in each color
    static void ElementColoring(TPZCompMesh *cmesh, TPZVec<long> &elSequence, TPZVec<long> &elSequenceColor,
                                TPZVec<long> &elBlocked, TPZVec<long> &NumelColors);

protected:
    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs);
    

protected:
    
    /** @brief Structure to manipulate thread to solve system equations */
    struct ThreadData
    {
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZAssembleOptimizedColor *strmat,int seqnum, TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs);
        /** @brief Initialize the mutex semaphores and others */
        ThreadData(TPZAssembleOptimizedColor *strmat, int seqnum, TPZFMatrix<STATE> &rhs);
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
        TPZAssembleOptimizedColor *fStruct;
        /** @brief Global matrix */
        TPZMatrix<STATE> *fGlobMatrix;
        /** @brief Global rhs vector */
        TPZFMatrix<STATE> *fGlobRhs;
        
#ifdef USING_BOOST
        boost::atomic<long> *fCurrentIndex;
#else
#endif
        /** @brief sequence number of the thread */
        int fThreadSeqNum;
        
        /** @brief vector indicating whether an element has been computed */
        TPZVec<long> *fComputedElements;
        /** @brief Mutexes (to choose which element is next) */
        pthread_mutex_t *fAccessElement;
        
        pthread_cond_t *fCondition;
        
        int *fSomeoneIsSleeping;
        
        /// Vector for mesh coloring
        TPZVec<long> *fElBlocked, *fElSequenceColor;
        
        /// All elements below or equal this index have been computed
        long *fElementCompleted;
        
        static void *ThreadWorkResidual(void *datavoid);
    };
    
    friend struct ThreadData;
protected:
    
    /** @brief Vectors for mesh coloring */
    TPZVec<long> fElBlocked, fElSequenceColor;
    
    /// vector of the size of the elements containing 0 or 1 if the element has been computed (in the order of computation sequence)
    TPZVec<long> fElementsComputed;
    
    /// All elements below or equal this index have been computed
    long fElementCompleted;
    
    /// variable indicating if a thread is sleeping
    int fSomeoneIsSleeping;
    
#ifdef USING_BOOST
    boost::atomic<long> fCurrentIndex;
#endif
    
    
    /** @brief Mutexes (to choose which element is next) */
    pthread_mutex_t fAccessElement;
    
    pthread_cond_t fCondition;
    
    
    
protected:
    
};

#endif /* TPZAssemble_hpp */
