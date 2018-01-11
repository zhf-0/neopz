//
// Created by Francisco Teixeira Orlandini on 11/23/17.
//

#ifndef PZ_TPZEIGENSOLVER_H
#define PZ_TPZEIGENSOLVER_H

#include "TPZSolver.h"
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

  /**
   * @brief Solves the Ax=w*x eigenvalue problem and calculates the eigenvectors
   * @param A The matrix (input)
   * @param w Stores the eigenvalues
   * @param eigenVectors Stores the correspondent eigenvectors
   * @return it returns 1 if executed correctly
   */
  virtual int SolveEigenProblem(TPZMatrix <TVar> &A, TPZVec<typename SPZAlwaysComplex<TVar>::type> &w,
                                TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &eigenVectors) = 0;

  /**
   * @brief Solves the Ax=w*x eigenvalue problem and does not calculate the eigenvectors
   * @param A The matrix (input)
   * @param w Stores the eigenvalues
   * @return it returns 1 if executed correctly
   */
  virtual int SolveEigenProblem(TPZMatrix <TVar> &A, TPZVec<typename SPZAlwaysComplex<TVar>::type> &w) = 0;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param A The lhs matrix (input)
   * @param B The rhs matrix (input)
   * @param w Stores the eigenvalues
   * @param eigenVectors Stores the correspondent eigenvectors
   * @return it returns 1 if executed correctly
   */
  virtual int SolveGeneralisedEigenProblem(TPZMatrix <TVar> &A, TPZMatrix <TVar> &B,
                                           TPZVec<typename SPZAlwaysComplex<TVar>::type> &w,
                                           TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &eigenVectors) = 0;

  /**
   * @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param A The lhs matrix (input)
   * @param B The rhs matrix (input)
   * @param w Stores the eigenvalues
   * @return it returns 1 if executed correctly
   */
  virtual int SolveGeneralisedEigenProblem(TPZMatrix <TVar> &A, TPZMatrix <TVar> &B,
                                           TPZVec<typename SPZAlwaysComplex<TVar>::type> &w) = 0;

  void Solve(TPZVec<typename SPZAlwaysComplex<TVar>::type> &eigenValues, TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &eigenVectors);

  void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
             TPZFMatrix<TVar>  *residual = 0) override {
    DebugStop();
  }

  int ClassId() const override;//implement me

  TPZAutoPointer<TPZMatrix<TVar>> MatrixA();

  TPZAutoPointer<TPZMatrix<TVar> > MatrixB();

  void SetMatrixA(TPZAutoPointer<TPZMatrix<TVar>> mat);

  void SetMatrixB(TPZAutoPointer<TPZMatrix<TVar>> mat);

  virtual bool IsGeneralised() const;//just set its attribute. must be implemented in child class

  virtual void SetAsGeneralised(bool isGeneralised);//just set its attribute. must be implemented in child class

  virtual int GetHowManyEigenvalues() const;//just set its attribute. must be implemented in child class

  virtual void SetHowManyEigenValues(int howManyEigenValues);//just set its attribute. must be implemented in child class

  virtual bool IsAbsoluteValue();//just set its attribute. must be implemented in child class

  virtual void SetAbsoluteValue(bool isAbsoluteValue);//just set its attribute. must be implemented in child class

  virtual EDesiredEigen GetDesiredPartOfSpectrum() const;//just set its attribute. must be implemented in child class

  virtual void SetDesiredPartOfSpectrum(EDesiredEigen desiredPartOfSpectrum);//just set its attribute. must be implemented in child class

  virtual typename SPZAlwaysComplex<TVar>::type GetSpecifiedValue() const;//just set its attribute. must be implemented in child class

  virtual void SetSpecifiedValue(typename SPZAlwaysComplex<TVar>::type specifiedValue);//just set its attribute. must be implemented in child class
protected:
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
};


#endif //PZ_TPZEIGENSOLVER_H
