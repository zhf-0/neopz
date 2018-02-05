//
// Created by Francisco Teixeira Orlandini on 2/2/18.
// This class is poorly designed at the moment. Ideally, the solver and
// the preconditioner should both have their own classes. The SetType
// method could be better (currently it is designed thinking of using ShiftInvert).
//

#ifndef PZ_TPZSLEPCSTHANDLER_H
#define PZ_TPZSLEPCSTHANDLER_H

#include <pzreal.h>
#ifdef USING_SLEPC

#include <slepceps.h>
#include "TPZPetscWrapper.h"
#include "TPZSlepcEPSHandler.h"

class TPZSlepcSTHandler : public TPZPetscWrapper{
  template<typename TVar>
  friend class TPZSlepcEPSHandler;
public:
  TPZSlepcSTHandler();
  ~TPZSlepcSTHandler();
  void SetType(const STType &st, const STATE &shift = 0.);
  void SetSolver(const KSPType &solver);
  void SetSolverTol(const PetscReal &rtol, const PetscReal &atol, const PetscReal &dtol, const PetscInt &max_its);
  void GetSolverTol(PetscReal *rtol, PetscReal *atol, PetscReal *dtol, PetscInt *max_its);
  void SetPrecond(const PCType &precond);
private:
  ST fSt;
  KSP fKsp;
  PC fPc;
  bool fHasSolver;
  bool fHasPrecond;
};

#endif //USING_SLEPC
#endif //PZ_TPZSLEPCSTHANDLER_H
