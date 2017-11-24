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

    TPZEigenSolver<STATE> * GetSolver() const;

    void SetSolver(TPZEigenSolver<STATE> * &fSolver);

    void Assemble() override;

    void Solve() override;
protected:
    /**
     * @brief Pointer to the Eigen solver
     */
    TPZEigenSolver<STATE> * fSolver;

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
    using TPZAnalysis::BuildPreconditioner;
    using TPZAnalysis::AnimateRun;
    using TPZAnalysis::PostProcessTable;
    using TPZAnalysis::LoadSolution;
    using TPZAnalysis::SetExact;
    using TPZAnalysis::PostProcessError;
    using TPZAnalysis::PostProcessErrorSerial;
    using TPZAnalysis::PostProcessErrorParallel;
    using TPZAnalysis::CreateListOfCompElsToComputeError;
};


#endif //PZ_TPZEIGENANALYSIS_H
