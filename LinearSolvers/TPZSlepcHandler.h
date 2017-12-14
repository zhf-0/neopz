//
// Created by Francisco Teixeira Orlandini on 12/13/17.
//

#ifndef PZ_TPZSLEPCHANDLER_H
#define PZ_TPZSLEPCHANDLER_H

#include <pzysmp.h>
#include <TPZEigenHandler.h>

template<class TVar>
class SPZAlwaysComplex;


template<class TVar>
class TPZSlepcHandler : public TPZEigenHandler<TVar> {
  friend class TPZFYsmpMatrix<TVar>;
public:
  int SolveEigenProblem(TPZMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors) override;

  /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param w Stores the eigenvalues
   */
  int SolveEigenProblem(TPZMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w) override;

  /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param w Stores the eigenvalues
   * @param Stores the correspondent eigenvectors
   */
  int SolveGeneralisedEigenProblem(TPZMatrix<TVar> &A, TPZMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors) override;
  /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param w Stores the eigenvalues
   */
  int SolveGeneralisedEigenProblem(TPZMatrix<TVar> &A, TPZMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w) override;

  /*******************
  *    TPZFYSMPMATRIX    *
  *******************/
  int SolveEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors);

  /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param w Stores the eigenvalues
   */
  int SolveEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w);

  /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
   * @param w Stores the eigenvalues
   * @param Stores the correspondent eigenvectors
   */
  int SolveGeneralisedEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZFYsmpMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors);
  /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
   * @param w Stores the eigenvalues
   */
  int SolveGeneralisedEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZFYsmpMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w);

private:
  bool fAmIInitialised = false;
};


#endif //PZ_TPZSLEPCHANDLER_H
