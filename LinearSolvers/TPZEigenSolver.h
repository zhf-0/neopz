//
// Created by Francisco Teixeira Orlandini on 11/23/17.
//

#ifndef PZ_TPZEIGENSOLVER_H
#define PZ_TPZEIGENSOLVER_H

#include "TPZSolver.h"
#include "TPZEigenHandler.h"


template<typename TVar>
class TPZEigenHandler;
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
* @ingroup solver
* @brief  Defines a class of solvers for eigenvalue problems. \ref solver "Solver"
*/
template <typename TVar>
class TPZEigenSolver : public TPZSolver<TVar> {
public:
  TPZEigenSolver();

  TPZEigenSolver(const TPZEigenSolver &copy);

  void Solve(TPZVec<typename SPZAlwaysComplex<TVar>::type> &eigenValues, TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &eigenVectors);

  void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
             TPZFMatrix<TVar>  *residual = 0) override {
    DebugStop();
  }

  TPZSolver<TVar> *Clone() const override{
    return new TPZEigenSolver<TVar>(*this);
  }

  virtual int ClassId() const override;

  TPZAutoPointer<TPZMatrix<TVar>> MatrixA();

  TPZAutoPointer<TPZMatrix<TVar> > MatrixB();

  void SetMatrixA(TPZAutoPointer<TPZMatrix<TVar>> mat);

  void SetMatrixB(TPZAutoPointer<TPZMatrix<TVar>> mat);

  bool IsGeneralised() const;

  void SetAsGeneralised(bool isGeneralised);

  int GetHowManyEigenvalues() const;

  void SetHowManyEigenValues(int howManyEigenValues);

  TPZEigenHandler<TVar> &GetEigenHandler() const;

  void SetEigenHandler(TPZEigenHandler<TVar> &eig);

  void SetEigenHandler(TPZAutoPointer<TPZEigenHandler<TVar>> &eig);

  bool IsAbsoluteValue();

  void SetAbsoluteValue(bool isAbsoluteValue);

  EDesiredEigen GetDesiredPartOfSpectrum() const;

  void SetDesiredPartOfSpectrum(EDesiredEigen desiredPartOfSpectrum);

  typename SPZAlwaysComplex<TVar>::type GetSpecifiedValue() const;

  void SetSpecifiedValue(typename SPZAlwaysComplex<TVar>::type specifiedValue);
protected:

  void AutoSetEigenHandler();
  /**
   * @brief Whether to display the absolute value(true) or the real part
   * of the eigenvectors over the computational mesh
   */
  bool fShowAbsoluteValue;
  /** @brief Whether to solve the eigenvalue problem
       *   is generalised (Ax=uBx) or not (Ax=ux)*/
  bool fIsGeneralised;
  /**
   * @brief Whether to calculate the eigenvectors (and not the eigenvalues only)
   */
  bool fMustCalculateEigenVectors;
  /** @brief Desired number of eigenvalues to be computed*/
  int fHowManyEigenValues;
  /**
   * @brief Where in the spectrum to search for eigenvalues
   */
  EDesiredEigen fDesiredPartOfSpectrum = MNE;
  /**
   * @brief If fDesiredPartOfSpectrum is SVCP, eigenvalues will be
   * searched for around this value. It is always complex, regardless of
   * what type STATE refers to.
   */
  typename SPZAlwaysComplex<TVar>::type fSpecifiedValue;

  /**
       * @brief Stores the computed eigenvalues
       */
  TPZManVector<typename SPZAlwaysComplex<TVar>::type,10> fEigenvalues;

  /**
   * @brief Stores the computed eigenvectors
   */
  TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> fEigenvectors;

  /** @brief Container classes */
  TPZAutoPointer<TPZMatrix<TVar> > fMatrixA;

  /** @brief Container classes */
  TPZAutoPointer<TPZMatrix<TVar> > fMatrixB;

  TPZAutoPointer<TPZEigenHandler<TVar>> fEig;
};


#endif //PZ_TPZEIGENSOLVER_H
