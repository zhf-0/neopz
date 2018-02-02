//
// Created by Francisco Teixeira Orlandini on 2/2/18.
//

#ifdef USING_PETSC
#include <petscsys.h>
#ifdef USING_SLEPC
#include <slepcsys.h>
#endif//USING_SLEPC
#include "TPZPetscWrapper.h"

TPZPetscWrapper::TPZPetscWrapper() {
  static bool firstTimeHere = true;
  if(firstTimeHere){
    #ifdef USING_SLEPC
    SlepcInitialize((int *)0, (char ***)0, (const char*)0,(const char*)0 );
    #else
    PetscInitialize((int *)0, (char ***)0, (const char*)0,(const char*)0 );
    #endif
  }
}
#endif//USING_PETSC
