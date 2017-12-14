//
// Created by Francisco Teixeira Orlandini on 12/14/17.
//

#ifndef PZ_TPZEIGENHANDLER_H
#define PZ_TPZEIGENHANDLER_H
#include <pzmatrix.h>
#include <pzfmatrix.h>
template<class TVar>
class SPZAlwaysComplex;

template<class TVar>
class TPZEigenHandler {
public:
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
};
#endif //PZ_TPZEIGENHANDLER_H
