//
// Created by Francisco Teixeira Orlandini on 11/17/17.
//

#ifndef PZ_TPZEIGENANALYSIS_H
#define PZ_TPZEIGENANALYSIS_H
#include "pzanalysis.h"     //For TPZAnalysis
#include "tpzautopointer.h" //For TPZAutoPointer
#include "pzmatrix.h"       //For TPZFMatrix

template<typename TVar>
class TPZEigenSolver;

/**
 * @brief Specialization of TPZAnalysis dedicated to eigenvalue problems
 */
class TPZEigenAnalysis : public TPZAnalysis{
public:
    /** @brief Create an empty TPZEigenAnalysis object*/
    TPZEigenAnalysis();
    /** @brief Create an TPZEigenAnalysis object from one mesh pointer */
    TPZEigenAnalysis(TPZCompMesh *mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);
    /** @brief Create an TPZEigenAnalysis object from one mesh auto pointer object */
    TPZEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);

    TPZEigenSolver<STATE> &Solver() const;

    void SetSolver(TPZEigenSolver<STATE> &fSolver);

    void Assemble() override;

    void Solve() override;

    //@TODO: IMPLEMENTAR CLASSID
    // virtual int ClassId() const override;

    TPZFMatrix<SPZAlwaysComplex<STATE>::type> GetEigenvectors() const;
    TPZManVector<SPZAlwaysComplex<STATE>::type> GetEigenvalues() const;
protected:
    template<class TVar>
    void TransferEigenVector(const TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &matrix, TPZFMatrix<TVar> &tpzfMatrix,
                             const unsigned int &i, const bool isAbsVal);
    /**
     * @brief Pointer to the Eigen solver
     */
    TPZEigenSolver<STATE> *fSolver;
    /**
    * @brief Stores the computed eigenvalues
    */
    TPZManVector<SPZAlwaysComplex<STATE>::type,10> fEigenvalues;
    /**
     * @brief Stores the computed eigenvectors
     */
    TPZFMatrix<SPZAlwaysComplex<STATE>::type> fEigenvectors;

/**************************************************************************************
 * These private members and deleted functions are not related to TPZEigenAnalysis.   *
 * Until TPZAnalysis is refactored, they are set as private in order to not be used   *
 **************************************************************************************/
    TPZFMatrix<STATE> &Rhs() = delete;
    void SetStep(int step) = delete;
    int GetStep() = delete;
    void SetTime(REAL time) = delete;
    TPZAlgebraicSystemSolver<STATE> *BuildPreconditioner(EPrecond preconditioner, bool overlap) = delete;
    void AnimateRun(long num_iter, int steps, TPZVec<std::string> &scalnames,
                    TPZVec<std::string> &vecnames, const std::string &plotfile) = delete;
    void CreateListOfCompElsToComputeError(TPZAdmChunkVector<TPZCompEl *> &elvecToComputeError) = delete;
private:
    //using TPZAnalysis::PostProcess;
    using TPZAnalysis::PrePostProcessTable;
    using TPZAnalysis::PostProcessTable;
    using TPZAnalysis::PostProcessError;
    using TPZAnalysis::PostProcessErrorSerial;
    using TPZAnalysis::PostProcessErrorParallel;
};


#endif //PZ_TPZEIGENANALYSIS_H
