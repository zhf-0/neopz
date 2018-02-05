//
// Created by Francisco Teixeira Orlandini on 12/13/17.
//

#ifndef PZ_TPZSLEPCHANDLER_H
#define PZ_TPZSLEPCHANDLER_H

#ifdef USING_SLEPC
#include "TPZPetscWrapper.h"
#include <pzysmp.h>
#include <TPZEigenSolver.h>
#include <slepceps.h>

template<class TVar>
class SPZAlwaysComplex;

class TPZSlepcSTHandler;

template<class TVar>
class TPZSlepcEPSHandler : public TPZEigenSolver<TVar> , TPZPetscWrapper {
  friend class TPZFYsmpMatrix<TVar>;
public:
  TPZSlepcEPSHandler();
  ~TPZSlepcEPSHandler();
  TPZSolver<TVar> *Clone() const override{
    //@TODO: Implement me!
    return (TPZSolver<TVar> *)this;
  }

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

  void SetST(const TPZSlepcSTHandler &st);

  void GetST(TPZSlepcSTHandler *st);

  void SetProblemType(const EPSProblemType epsProblem);

  EPSProblemType GetProblemType() const;

  void SetTolerances(const PetscReal &tol, const PetscInt &max_its);

  void GetTolerances(PetscReal *tol, PetscInt *max_its) const;

  void SetEPSDimensions(const PetscInt &nev, const PetscInt &ncv, const PetscInt &mpd);

  void GetEPSDimensions(PetscInt *nev, PetscInt *ncv, PetscInt *mpd) const;

  void SetConvergenceTest(const EPSConv&test);

  EPSConv GetConvergenceTest() const;

  void SetTrueResidual(const bool &opt);

  void SetType(const EPSType &type);

  EPSType GetType() const;

  void SetVerbose(bool fVerbose);

  void SetTargetEigenvalue(const PetscScalar &target);

  PetscScalar GetTargetEigenvalue() const;

  void SetWhichEigenpairs(const EPSWhich eps_which);

  EPSWhich GetWhichEigenpairs() const;

  void SetKrylovOptions(const bool &pLocking, const PetscReal &restart);

  void GetKrylovOptions(bool *pLocking , PetscReal *restart) const;

private:

  bool fVerbose = true;

  EPS fEps;
};

#endif //USING_SLEPC
#endif //PZ_TPZSLEPCHANDLER_H
