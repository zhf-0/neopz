/**
 * @file
 * @brief Contains the TPZStructMatrixOR class which responsible for a interface among Matrix and Finite Element classes.
 */

#ifndef TPZStructMatrixOR_H
#define TPZStructMatrixOR_H

#include <set>
#include <map>
#include <semaphore.h>
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"
#include "pzequationfilter.h"
#include "TPZGuiInterface.h"

class TPZAssemble;

class TPZCompMesh;
template<class TVar>
class TPZMatrix;
template<class TVar>
class TPZFMatrix;

/**
 * @brief Refines geometrical mesh (all the elements) num times
 * @ingroup geometry
 */
//void UniformRefine(int num, TPZGeoMesh &m);

/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrix {
    
public:
    
    enum MAssembleStyle {EPThread, ETBB, EColoring};
    
    TPZStructMatrix(TPZCompMesh *);
    
    TPZStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);
    
    TPZStructMatrix(const TPZStructMatrix &copy);
    
    virtual ~TPZStructMatrix();
    
    /** @brief Sets number of threads in Assemble process */
    void SetNumThreads(int n){
        this->fAssembleConfig.fNumThreads = n;
    }
    
    int GetNumThreads() const{
        return this->fAssembleConfig.fNumThreads;
    }
    
    virtual TPZMatrix<STATE> * Create();
    
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs);
    
    virtual TPZStructMatrix * Clone();
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs);
    
    /** @brief Assemble the global right hand side */
    virtual void Assemble(TPZFMatrix<STATE> & rhs);
    
public:
    
    /** @brief Determine that the assembly refers to a range of equations */
    void SetEquationRange(long mineq, long maxeq)
    {
        fAssembleConfig.fEquationFilter.Reset();
        fAssembleConfig.fEquationFilter.SetMinMaxEq(mineq, maxeq);
    }
    
    /** @brief Verify if a range has been specified */
    virtual bool HasRange() const
    {
        return fAssembleConfig.fEquationFilter.IsActive();
    }
    
    /** @brief access method for the equation filter */
    TPZEquationFilter &EquationFilter()
    {
        return fAssembleConfig.fEquationFilter;
    }
    
    /** @brief number of equations after applying the filter */
    long NReducedEquations() const
    {
        return fAssembleConfig.fEquationFilter.NActiveEquations();
    }
    
    /** @brief Access method for the mesh pointer */
    TPZCompMesh *Mesh() const
    {
        return fAssembleConfig.fMesh;
    }
    
    /** @brief Set the set of material ids which will be considered when assembling the system */
    void SetMaterialIds(const std::set<int> &materialids);
    
    /** @brief Establish whether the element should be computed */
    bool ShouldCompute(int matid) const
    {
        const unsigned int size = fAssembleConfig.fMaterialIds.size();
        return size == 0 || fAssembleConfig.fMaterialIds.find(matid) != fAssembleConfig.fMaterialIds.end();
    }
    /** @brief Returns the material ids */
    const std::set<int> &MaterialIds()
    {
        return fAssembleConfig.fMaterialIds;
    }
    
protected:
    
    
    friend struct ThreadData;
public:
    struct TPZAssembleConfig
    {
        
        TPZAssembleConfig() : fMesh(0), fCompMesh(0), fEquationFilter(), fMaterialIds(), fNumThreads(0)
        {
            
        }
        
        TPZAssembleConfig(const TPZAssembleConfig &copy) : fMesh(copy.fMesh), fCompMesh(copy.fCompMesh), fEquationFilter(copy.fEquationFilter),
            fMaterialIds(copy.fMaterialIds), fNumThreads(copy.fNumThreads)
        {
        }
        
        TPZAssembleConfig &operator=(const TPZAssembleConfig &copy)
        {
            fMesh = copy.fMesh;
            fCompMesh = copy.fCompMesh;
            fEquationFilter = copy.fEquationFilter;
            fMaterialIds = copy.fMaterialIds;
            fNumThreads = copy.fNumThreads;
            return *this;
        }
        /** @brief Pointer to the computational mesh from which the matrix will be generated */
        TPZCompMesh * fMesh;
        /** @brief Autopointer control of the computational mesh */
        TPZAutoPointer<TPZCompMesh> fCompMesh;
        /** @brief Object which will determine which equations will be assembled */
        TPZEquationFilter fEquationFilter;
        
        /** @brief Set of material ids to be considered. It is a private attribute. */
        /** Use ShouldCompute method to know if element must be assembled or not    */
        std::set<int> fMaterialIds;
        
        /** @brief Number of threads in Assemble process */
        int fNumThreads;
    };
    
protected:
    
    TPZAssembleConfig fAssembleConfig;
    
    TPZAssemble *fStiffAssemble;
    
    TPZAssemble *fRhsAssemble;
    
};

#endif


/**
 * @file
 * @brief Contains the TPZStructMatrix class which responsible for a interface among Matrix and Finite Element classes.
 */

#ifndef TPZ_STRUCT_MATRIX_H
#define TPZ_STRUCT_MATRIX_H

#include "pzstrmatrixcs.h"
#include "pzstrmatrixgc.h"
#include "pzstrmatrixot.h"
#include "pzstrmatrixtbb.h"
#include "pzstrmatrixflowtbb.h"
#include "pzstrmatrixst.h"

/** This is the original and stable version of multi_thread_assemble (producer-consumer) */
//typedef TPZStructMatrixOR TPZStructMatrix;

/** This version has a clean code with openmp parallism */
//typedef TPZStructMatrixST TPZStructMatrix;

/** This version uses locks in the assemble contribuition with tbb (Nathan-Borin) */
//typedef TPZStructMatrixCS TPZStructMatrix;

/** This version uses graph coloring to define the order to process the elements (Devloo-Gilvan) */
//typedef TPZStructMatrixGC TPZStructMatrix;

/** This version uses graph coloring to define the order to process the elements (Devloo-Gilvan) and
 * each color is processed and syncronized */
//typedef TPZStructMatrixOT TPZStructMatrix;

/** This version uses the graph coloring and create a tbb::flow::graph to process in parallel */
//https://trac.macports.org/wiki/MigrationTBB
//typedef TPZStructMatrixTBB TPZStructMatrix;

/** This version uses the graph coloring and create a tbb::flow::graph to process in parallel 
 *  every node of the tbb flow graph computes calc and the assemble
 */
//typedef TPZStructMatrixTBBFlow TPZStructMatrix;

#endif
