//
// Created by Francisco Teixeira Orlandini on 2/2/18.
//

#include <pzerror.h>
#include "TPZSlepcSTHandler.h"
#ifdef USING_SLEPC

TPZSlepcSTHandler::TPZSlepcSTHandler() : fSt(0) , fKsp(0) , fPc(0) ,
                                         fHasSolver(false), fHasPrecond(false) {
  STCreate(PETSC_COMM_WORLD,&fSt);
}

TPZSlepcSTHandler::~TPZSlepcSTHandler() {
  STDestroy(&fSt);
  KSPDestroy(&fKsp);
  PCDestroy(&fPc);
}

void TPZSlepcSTHandler::SetType(const STType &st, const STATE &shift) {
  PetscErrorCode ierr = 1;
  ierr = STSetType(fSt, st);
  if(ierr != 0){
    PZError<<"Invalid spectral transformation\n";
    DebugStop();
  }
  if(strcmp(st,STSINVERT) == 0){
    STSetShift(fSt,shift);
  }
}

void TPZSlepcSTHandler::SetSolver(const KSPType &solver) {
  PetscErrorCode ierr = 1;
  if(!fHasSolver){
    KSPCreate(PETSC_COMM_WORLD,&fKsp);
    fHasSolver = true;
  }
  ierr = KSPSetType(fKsp, solver);
  if(ierr != 0){
    PZError<<"Invalid ST solver\n";
    DebugStop();
  }
  STSetKSP(fSt,fKsp);
}

void TPZSlepcSTHandler::SetSolverTol(const PetscReal &rtol, const PetscReal &atol, const PetscReal &dtol, const PetscInt &max_its) {
  if(!fHasSolver){
    KSPCreate(PETSC_COMM_WORLD,&fKsp);
    fHasSolver = true;
  }
  KSPSetTolerances(fKsp , rtol , atol , dtol , max_its);
  STSetKSP(fSt,fKsp);
}

void TPZSlepcSTHandler::GetSolverTol(PetscReal *rtol, PetscReal *atol, PetscReal *dtol, PetscInt *max_its) {
  KSP ksp;
  PetscErrorCode ierr = STGetKSP(fSt, &ksp);
  if(ierr) DebugStop();
  ierr = KSPGetTolerances(ksp, rtol , atol , dtol , max_its);
  if(ierr) DebugStop();
}

void TPZSlepcSTHandler::SetPrecond(const PCType &precond) {
  PetscErrorCode ierr = 1;
  if(!fHasPrecond){
    PCCreate(PETSC_COMM_WORLD,&fPc);
    fHasPrecond = true;
  }

  ierr = PCSetType(fPc, precond);
  if(ierr != 0) DebugStop();
  if(!fHasSolver){
    KSPCreate(PETSC_COMM_WORLD,&fKsp);
    fHasSolver = true;
  }
  KSPSetPC(fKsp,fPc);
}
#endif //USING_SLEPC
