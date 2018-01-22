
#include <pzsysmp.h>
#include <pzysmp.h>
#include "petscmat.h"
#include "slepceps.h"

int TestintTPZSYsmpMatrix();

int main(int argc, char **argv){

    #if !defined(PETSC_USE_COMPLEX)
    std::cout<<"ops"<<std::endl;
    #else
    std::cout<<"hehe"<<std::endl;
    #endif

  if(!TestintTPZSYsmpMatrix()){
    DebugStop();
  }
  return 0;
}

int TestintTPZSYsmpMatrix() {
  static char help[] = "Testing MatCreateSeqSBAIJWithArrays().\n\n";
  const int nRows = 3;
  const int nCols = nRows;

  TPZFYsmpMatrix<STATE> matrix(nRows,nCols);
  TPZVec<long> ia ({0,1,3,4});
  TPZVec<long> ja({0,1,2,2});

  TPZVec<STATE> aVec({2.,3.,4.,9.});
  matrix.SetData(ia,ja,aVec);
  matrix.Print("mat");
  const int blockSize = 1;
  Mat A;
  //@TODO:Think on int/long issue
  //@TODO:Think on PetsScalar issue
  TPZVec<int> iaP ({0,1,3,4});
  TPZVec<int> jaP({0,1,2,2});

  PetscErrorCode err;
  err=SlepcInitialize((int *)0,(char ***)0,(char *)0,help);

  if(err){
    DebugStop();
  }
  err=MatCreateSeqSBAIJWithArrays(MPI_COMM_WORLD,blockSize,nRows,nCols,(PetscInt *)iaP.begin(),(PetscInt *)jaP.begin(),(PetscScalar *)aVec.begin(),&A);
  if(err){
    DebugStop();
  }

  err = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATHEMATICA);
  if(err){
    DebugStop();
  }
  err=MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  if(err){
    DebugStop();
  }
  err=PetscViewerFlush(PETSC_VIEWER_STDOUT_WORLD);
  if(err){
    DebugStop();
  }
  EPS eps;
  EPSCreate(PETSC_COMM_WORLD,&eps);

  //Set operators. In this case, it is a standard eigenvalue problem
  EPSSetOperators(eps,A,NULL);
  EPSSetProblemType(eps,EPS_HEP);

  //Set solver parameters at runtime
  EPSSetFromOptions(eps);

  EPSSolve(eps);
  EPSType type;
  PetscInt nev;
  PetscBool terse;
  EPSGetType(eps,&type);
  PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);
  EPSGetDimensions(eps,&nev,NULL,NULL);
  PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  PetscOptionsHasName(NULL,NULL,"-terse",&terse);
  if (terse) {
    EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);
  } else {
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
    EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);
    EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);
    PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
  }
  EPSDestroy(&eps);
  MatDestroy(&A);
  SlepcFinalize();
  return 1;
}
