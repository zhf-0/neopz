//
// Created by Francisco Teixeira Orlandini on 11/23/17.
//

#ifndef PZ_TPZEIGENSOLVER_H
#define PZ_TPZEIGENSOLVER_H

#include "TPZSolver.h"

/**
* @ingroup solver
* @brief  Defines a class of solvers for eigenvalue problems. \ref solver "Solver"
*/
template <typename TVar>
class TPZEigenSolver : public TPZSolver<TVar> {
public:
    bool IsGeneralised() const;

    TPZAutoPointer<TPZMatrix<TVar>> MatrixA();

    TPZAutoPointer<TPZMatrix<TVar> > MatrixB();

    void SetMatrixA(TPZAutoPointer<TPZMatrix<TVar>> pointer);

    void SetMatrixB(TPZAutoPointer<TPZMatrix<TVar>> mat);

    void Solve(TPZVec<typename SPZAlwaysComplex<TVar>::type> &eigenValues, TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &eigenVectors);

    virtual int ClassId() const override;

    void SetAsGeneralised(bool fIsGeneralised);
protected:
    /** @brief Whether to solve the eigenvalue problem
         *   is generalised (Ax=uBx) or not (Ax=ux)*/
    bool fIsGeneralised = false;
    /**
     * @brief Whether to calculate the eigenvectors (and not the eigenvalues only)
     */
    bool fMustCalculateEigenVectors = false;
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
    SPZAlwaysComplex<TVar> fSpecifiedValue;

    /**
         * @brief Stores the computed eigenvalues
         */
    TPZFMatrix<TVar> fEigenvalues;

    /**
     * @brief Stores the computed eigenvectors
     */
    TPZFMatrix<TVar> fEigenvectors;

    /** @brief Container classes */
    TPZAutoPointer<TPZMatrix<TVar> > fMatrixA;

    /** @brief Container classes */
    TPZAutoPointer<TPZMatrix<TVar> > fMatrixB;
};


#endif //PZ_TPZEIGENSOLVER_H
