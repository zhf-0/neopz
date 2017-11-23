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

    const TPZAutoPointer<TPZStructMatrix> &GetAMatrix() const;

    void SetAMatrix(const TPZAutoPointer<TPZStructMatrix> &fAStructMatrix);

    const TPZAutoPointer<TPZStructMatrix> &GetBMatrix() const;

    void SetBMatrix(const TPZAutoPointer<TPZStructMatrix> &fBStructMatrix);

    TPZEigenSolver<STATE> * GetSolver() const;

    void SetSolver(TPZEigenSolver<STATE> * &fSolver);

    void Solve();
protected:
    /** @brief Structural matrix A (as in Ax = uBx)*/
    TPZAutoPointer<TPZStructMatrix>  fAStructMatrix;
    /** @brief Structural matrix B (as in Ax = uBx)*/
    TPZAutoPointer<TPZStructMatrix>  fBStructMatrix;
    /**
     * @brief Pointer to the Eigen solver
     */
    TPZEigenSolver<STATE> * fSolver;
    /** @brief Whether to solve the eigenvalue problem
     *   is generalised (Ax=uBx) or not (Ax=ux)*/
    bool fIsGeneralised = false;
    /** @brief Desired number of eigenvalues to be computed*/
    int fHowManyEigenValues = 1;
    /**
     * @brief Enum for defining ranges in the spectrum
     */
    enum EDesiredEigen {
        /** Most Negative EigenValues */
         MNE,
        /** Least Negative Eigenvalues */
        LNE,
        /** Least Positive Eigenvalues */
        LPE,
        /** Most Positive Eigenvalues */
        MPE,
        /** Specified Value on the Complex Plane */
        SVCP
    };
    /**
     * @brief Where in the spectrum to search for eigenvalues
     */
    EDesiredEigen fDesiredPartOfSpectrum = MNE;
    /**
     * @brief If fDesiredPartOfSpectrum is SVCP, eigenvalues will be
     * searched for around this value. It is always complex, regardless of
     * what type STATE refers to.
     */
    std::complex<REAL> fSpecifiedValue = 0.;

    /**
     * @brief Stores the computed eigenvalues
     */
    TPZFMatrix<STATE> fEigenvalues;
    /**
     * @brief Stores the computed eigenvectors
     */
    TPZFMatrix<STATE> fEigenvectors;

/************************************************************************************
 * These private members are not really related to TPZEigenAnalysis.                *
 * Until TPZAnalysis is refactored, they are set as private in order to not be used *
 ************************************************************************************/
private:
    using TPZAnalysis::Rhs;
    using TPZAnalysis::Solution;
    using TPZAnalysis::StructMatrix;
    using TPZAnalysis::SetStep;
    using TPZAnalysis::GetStep;
    using TPZAnalysis::SetTime;
    using TPZAnalysis::GetTime;
    using TPZAnalysis::PostProcess;
    using TPZAnalysis::PrePostProcessTable;
    using TPZAnalysis::PostProcessTable;
    using TPZAnalysis::LoadSolution;
    using TPZAnalysis::SetExact;
    using TPZAnalysis::PostProcessError;
    using TPZAnalysis::PostProcessErrorSerial;
    using TPZAnalysis::PostProcessErrorParallel;
    using TPZAnalysis::CreateListOfCompElsToComputeError;

    using TPZAnalysis::SetStructuralMatrix;
};


#endif //PZ_TPZEIGENANALYSIS_H
