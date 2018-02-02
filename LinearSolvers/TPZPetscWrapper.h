//
// Created by Francisco Teixeira Orlandini on 2/2/18.
// Every SLEPc/PETSc handler/wrapper in NeoPZ **MUST** inherit
// from this class
//

#ifndef PZ_TPZPETSCWRAPPER_H
#define PZ_TPZPETSCWRAPPER_H

#ifdef USING_PETSC
class TPZPetscWrapper{
protected:
  TPZPetscWrapper();
};

#endif//USING_PETSC
#endif //PZ_TPZPETSCWRAPPER_H
