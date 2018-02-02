//
// Created by Francisco Teixeira Orlandini on 12/13/17.
//
#include "TPZSlepcHandler.h"

#ifndef USING_SLEPC

void SlepcErrorMessage(){
  std::cout<<"SLEPc is not available. NeoPZ must be reconfigured."<<std::endl;
  DebugStop();
}
template<class TVar>
void TPZSlepcHandler::SetAsGeneralised(bool isGeneralised){
  SlepcErrorMessage();
}

template<class TVar>
void TPZSlepcHandler::SetHowManyEigenValues(int howManyEigenValues){
  SlepcErrorMessage();
}

template<class TVar>
void TPZSlepcHandler::SetAbsoluteValue(bool isAbsoluteValue){
  SlepcErrorMessage();
}

template<class TVar>
void TPZSlepcHandler::SetDesiredPartOfSpectrum(EDesiredEigen desiredPartOfSpectrum){
  SlepcErrorMessage();
}

template<class TVar>
void TPZSlepcHandler::SetSpecifiedValue(typename SPZAlwaysComplex<TVar>::type specifiedValue){
  SlepcErrorMessage();
}

template<class TVar>
int TPZSlepcHandler<TVar>::SolveEigenProblem(TPZMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors){
  SlepcErrorMessage();
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveEigenProblem(TPZMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  SlepcErrorMessage();
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveGeneralisedEigenProblem(TPZMatrix<TVar> &A, TPZMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
  SlepcErrorMessage();
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveGeneralisedEigenProblem(TPZMatrix<TVar> &A, TPZMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  SlepcErrorMessage();
  return 0;
}

#else
#include <slepceps.h>
#include <petsctime.h>
#include <petscksp.h>

template<class TVar>
void TPZSlepcHandler<TVar>::SetAsGeneralised(bool isGeneralised){
  TPZEigenSolver<TVar>::SetAsGeneralised(isGeneralised);
}

template<class TVar>
void TPZSlepcHandler<TVar>::SetHowManyEigenValues(int howManyEigenValues){
  TPZEigenSolver<TVar>::SetHowManyEigenValues(howManyEigenValues);
}

template<class TVar>
void TPZSlepcHandler<TVar>::SetAbsoluteValue(bool isAbsoluteValue){
  TPZEigenSolver<TVar>::SetAbsoluteValue(isAbsoluteValue);
}

template<class TVar>
void TPZSlepcHandler<TVar>::SetDesiredPartOfSpectrum(EDesiredEigen desiredPartOfSpectrum){
  TPZEigenSolver<TVar>::SetDesiredPartOfSpectrum(desiredPartOfSpectrum);
}

template<class TVar>
void TPZSlepcHandler<TVar>::SetSpecifiedValue(typename SPZAlwaysComplex<TVar>::type specifiedValue){
  TPZEigenSolver<TVar>::SetSpecifiedValue(specifiedValue);
}


/*******************
*    GENERAL       *
*******************/
template<class TVar>
int TPZSlepcHandler<TVar>::SolveEigenProblem(TPZMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors){

  TPZFYsmpMatrix<TVar> *Afsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&A);
  if(Afsparse){
    return SolveEigenProblem((TPZFYsmpMatrix<TVar>&)A,w,eigenVectors);
  }
  else{
    std::cout<<"TPZSlepcHandler does not support this matrix format"<<std::endl;
    DebugStop();
  }
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveEigenProblem(TPZMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  TPZFYsmpMatrix<TVar> *Afsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&A);
  if(Afsparse){
    return SolveEigenProblem((TPZFYsmpMatrix<TVar>&)A,w);
  }
  else{
    std::cout<<"TPZSlepcHandler does not support this matrix format"<<std::endl;
    DebugStop();
  }
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveGeneralisedEigenProblem(TPZMatrix<TVar> &A, TPZMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
  TPZFYsmpMatrix<TVar> *Afsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&A);
  TPZFYsmpMatrix<TVar> *Bfsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&B);
  if(Afsparse && Bfsparse){
    return SolveGeneralisedEigenProblem((TPZFYsmpMatrix<TVar>&)A,(TPZFYsmpMatrix<TVar>&)B,w,eigenVectors);
  }
  else{
    std::cout<<"TPZSlepcHandler does not support this matrix format"<<std::endl;
    DebugStop();
  }
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveGeneralisedEigenProblem(TPZMatrix<TVar> &A, TPZMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  TPZFYsmpMatrix<TVar> *Afsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&A);
  TPZFYsmpMatrix<TVar> *Bfsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&B);
  if(Afsparse && Bfsparse){
    return SolveGeneralisedEigenProblem((TPZFYsmpMatrix<TVar>&)A,(TPZFYsmpMatrix<TVar>&)B,w);
  }
  else{
    std::cout<<"TPZSlepcHandler does not support this matrix format"<<std::endl;
    DebugStop();
  }
  return 0;
}


/*******************
*    TPZFYSMPMATRIX    *
*******************/
template<class TVar>
int TPZSlepcHandler<TVar>::SolveGeneralisedEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZFYsmpMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
  /*****************************
   *  INITIALIZE STRUCTURES
   *****************************/
  SlepcInitialize((int *)0, (char ***)0, (const char*)0,(const char*)0 );
  /*****************************
   *  CREATE MATRICES
   *****************************/
  PetscReal error;
  const int blockSize = 1;
  const int nRows = A.Rows();
  const int nCols = A.Cols();

//	for (int i = 0; i < nRows; ++i) {
//		for (int j = 0; j < nRows; ++j) {
//
//			if(std::imag(A.GetVal(i,j))!=0){
//				std::cout<<"A("<<i<<","<<j<<") = "<<A.GetVal(i,j)<<std::endl;
//			}
//			if(std::imag(B.GetVal(i,j))!=0){
//				std::cout<<"A("<<i<<","<<j<<") = "<<A.GetVal(i,j)<<std::endl;
//			}
//		}
//	}
  Mat fAmat, fBmat;
  PetscErrorCode ierr;
  PetscInt *fIaP, *jaP;
  std::cout<<"Creating PETSc Amat...";
  ierr = PetscMalloc1(A.fIA.size(),&fIaP);CHKERRQ(ierr);
  ierr = PetscMalloc1(A.fJA.size(),&jaP);CHKERRQ(ierr);

  for (int j = 0; j < A.fIA.size(); ++j) {
    fIaP[j]=A.fIA[j];
  }
  for (int j = 0; j < A.fJA.size(); ++j) {
    jaP[j]=A.fJA[j];
  }
  ierr = MatCreateSeqBAIJWithArrays(MPI_COMM_WORLD,blockSize,nRows,nCols,fIaP,jaP,(PetscScalar *)A.fA.begin(),&fAmat);CHKERRQ(ierr);

  std::cout<<"Created!"<<std::endl;
  std::cout<<"Creating PETSc fBmat...";
  PetscInt *ibP, *jbP;
  ierr = PetscMalloc1(B.fIA.size(),&ibP);CHKERRQ(ierr);
  ierr = PetscMalloc1(B.fJA.size(),&jbP);CHKERRQ(ierr);
  for (int j = 0; j < B.fIA.size(); ++j) {
    ibP[j]=B.fIA[j];
  }
  for (int j = 0; j < B.fJA.size(); ++j) {
    jbP[j]=B.fJA[j];
  }

  ierr=MatCreateSeqBAIJWithArrays(MPI_COMM_WORLD,blockSize,nRows,nCols,ibP,jbP,(PetscScalar *)B.fA.begin(),&fBmat);CHKERRQ(ierr);

  std::cout<<"Created!"<<std::endl;

  /*****************************
   *  DEFINE PROBLEM TYPE
   *****************************/
  PC pc;
  KSP ksp;
  ST st;
  EPS eps;

  //EPSSetType(eps,EPSPOWER);
  //EPSSetType(eps,EPSLANCZOS);
  //EPSSetType(eps,EPSARPACK);
  //EPSSetType(eps,EPSGD);
  //EPSSetType(eps,EPSJD);
  //EPSSetType(eps,EPSCISS);

  //PETSc: Preconditioner class:
  PCCreate( PETSC_COMM_WORLD, &pc );
  PCSetType(pc, PCLU);
  PCSetFromOptions(pc);


//PETSc: Solver class (uses Preconditioner class):
  KSPCreate( PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPPREONLY);
  KSPSetPC(ksp, pc); // set PC created above
  //KSPSetConvergenceTest();
  KSPSetFromOptions(ksp);


//SLEPc: SpectralTransformation class (uses Solver class)
  STCreate( PETSC_COMM_WORLD, &st );
  STSetType(st,STSINVERT);
  const PetscScalar shift = -600000.;
  STSetShift(st, shift); // some specific settings to transformation
  STSetKSP(st, ksp); // set KSP created above


//SLEPc: Eigensolver class (uses SpectralTransformation class)

  PetscReal eps_tol;
  PetscInt eps_max_its;
  EPSCreate( PETSC_COMM_WORLD, &eps );
  EPSGetTolerances(eps, &eps_tol, &eps_max_its);
  eps_tol = 1e-10;
  eps_max_its = 100;
  EPSSetTolerances(eps,eps_tol,eps_max_its);
  EPSSetConvergenceTest(eps,EPS_CONV_REL);
  EPSSetWhichEigenpairs(eps, EPS_TARGET_REAL);
  EPSSetTarget(eps,shift);
  EPSSetST(eps,st); // set spectral transformation created above
  //EPSSetTrueResidual(eps,PETSC_TRUE);
  EPSSetOperators(eps, fAmat, fBmat);
  EPSSetProblemType(eps, EPS_GNHEP);
  EPSSetType(eps,EPSKRYLOVSCHUR);
  EPSKrylovSchurSetLocking(eps,PETSC_FALSE);
  EPSKrylovSchurSetRestart(eps,0.5);
  const PetscInt nev = 5;
  const PetscInt ncv = 30;
  //const PetscInt mpd = ncv - nev;
  EPSSetDimensions(eps, nev, ncv, PETSC_DECIDE);
  EPSSetFromOptions(eps);
  /*****************************
   *  SOLVE
   *****************************/
  PetscLogDouble t1,t2;
  ierr = PetscTime(&t1);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscTime(&t2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Elapsed Time in EPSSolve: %f\n",t2-t1);CHKERRQ(ierr);
  /*
     Optional: Get some information from the solver and display it
  */
  PetscInt its,lits,maxit;
  PetscReal tol;
  EPSType type;

  ierr = EPSGetIterationNumber(eps, &its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPGetTotalIterations(ksp, &lits);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %D\n",lits);CHKERRQ(ierr);
  ierr = EPSGetType(eps, &type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  //ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Show detailed info unless -terse option is given by user
   */
  const bool verbose = true;
  if (verbose) {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  } else {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  }
  

  return 1;
}

template<class TVar>
int TPZSlepcHandler<TVar>::SolveGeneralisedEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZFYsmpMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  DebugStop();
  return 0;
}

template<class TVar>
int TPZSlepcHandler<TVar>::SolveEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
  DebugStop();
  return 0;
}

template<class TVar>
int TPZSlepcHandler<TVar>::SolveEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  DebugStop();
  return 0;
}

#endif

template class TPZSlepcHandler< std::complex<float> >;
template class TPZSlepcHandler< std::complex<double> >;
template class TPZSlepcHandler< std::complex<long double> >;

template class TPZSlepcHandler<float >;
template class TPZSlepcHandler<double >;
template class TPZSlepcHandler<long double>;
