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

    void SetMatrixA(TPZAutoPointer<TPZMatrix<STATE>> pointer);

    void SetMatrixB(TPZAutoPointer<TPZMatrix<STATE>> mat);

    void Solve(TPZVec<TVar> &eigenValues, TPZFMatrix<TVar> &eigenVectors);

    virtual int ClassId() const override;
protected:
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

    /** @brief Container classes */
    TPZAutoPointer<TPZMatrix<TVar> > fMatrixA;

    /** @brief Container classes */
    TPZAutoPointer<TPZMatrix<TVar> > fMatrixB;
};


#endif //PZ_TPZEIGENSOLVER_H
