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
  EPS eps;
  EPSCreate( PETSC_COMM_WORLD, &eps);
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
  bool fIsParallelStorage;
  if (fIsParallelStorage)
  {
    ierr = MatCreate(PETSC_COMM_WORLD,&fAmat);CHKERRQ(ierr);
    ierr = MatSetSizes(fAmat,PETSC_DECIDE,PETSC_DECIDE,nRows,nCols);CHKERRQ(ierr);
    ierr = MatSetType(fAmat,MATMPIAIJ);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocationCSR(fAmat,fIaP,jaP,(PetscScalar *)A.fA.begin());CHKERRQ(ierr);
    ierr = PetscFree(fIaP);CHKERRQ(ierr);
    ierr = PetscFree(jaP);CHKERRQ(ierr);
  }
  else
  {
    error=MatCreateSeqBAIJWithArrays(MPI_COMM_WORLD,blockSize,nRows,nCols,fIaP,jaP,(PetscScalar *)A.fA.begin(),&fAmat);
  }
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
  if(fIsParallelStorage)
  {
    ierr = MatCreate(PETSC_COMM_WORLD,&fBmat);CHKERRQ(ierr);
    ierr = MatSetSizes(fBmat,PETSC_DECIDE,PETSC_DECIDE,nRows,nCols);CHKERRQ(ierr);
    ierr = MatSetType(fBmat,MATMPIAIJ);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocationCSR(fBmat,ibP,jbP,(PetscScalar *)B.fA.begin());CHKERRQ(ierr);
    ierr = PetscFree(ibP);CHKERRQ(ierr);
    ierr = PetscFree(jbP);CHKERRQ(ierr);
  }
  else
  {
    ierr=MatCreateSeqBAIJWithArrays(MPI_COMM_WORLD,blockSize,nRows,nCols,ibP,jbP,(PetscScalar *)B.fA.begin(),&fBmat);CHKERRQ(ierr);
  }
  std::cout<<"Created!"<<std::endl;

  /*****************************
   *  DEFINE PROBLEM TYPE
   *****************************/
  EPSSetOperators(eps, fAmat, fBmat);
  EPSSetProblemType(eps, EPS_GNHEP);
  //EPSSetType(eps,EPSPOWER);
  //EPSSetType(eps,EPSLANCZOS);
  //EPSSetType(eps,EPSARPACK);
  //EPSSetType(eps,EPSGD);
  //EPSSetType(eps,EPSJD);
  //EPSSetType(eps,EPSCISS);

  /*****************************
   *  SET SOLVER OPTIONS
   *****************************/
  EPSSetType(eps,EPSKRYLOVSCHUR);
  RG rg;
  EPSGetRG(eps,&rg);
  {
    PetscScalar center = -150000;
    PetscReal radius=125000, vscale=1;
    RGSetType(rg,RGELLIPSE);
    RGEllipseSetParameters(rg,center,radius,vscale);
  }
//	{
//		RGSetType(rg,RGINTERVAL);
//		PetscReal ra=-300000,rb=-5000,ia=-1,ib=1;
//		RGIntervalSetEndpoints(rg,ra,rb,ia,ib);
//	}
  const PetscScalar target = -600000;
  EPSSetWhichEigenpairs(eps, 	EPS_TARGET_REAL);
  EPSSetTarget(eps, target);
  PetscInt nev = 8;
  PetscInt ncv = nRows;
  PetscInt mpd = nRows - nev;
  EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE);
//	ierr = EPSSetTolerances(eps, 1e+2, 6); CHKERRQ(ierr);


  ST st;
  EPSGetST(eps,&st);
  ierr = STSetType(st, STSINVERT); CHKERRQ(ierr);

  KSP ksp;
  STGetKSP(st,&ksp);
  ierr = KSPSetType(ksp, KSPPREONLY); CHKERRQ(ierr); // because we're using a direct solver (PCLU)
  KSPSetErrorIfNotConverged(ksp,PETSC_FALSE);
  PetscReal rtol, abstol, dtol;
  PetscInt maxits;
  ierr = KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxits);CHKERRQ(ierr);
  rtol=1e-2;//it was 1e-8
  //abstol=1e-6;//it was 1e-50
  ierr = KSPSetTolerances(ksp, rtol, abstol, dtol, maxits);CHKERRQ(ierr);

  PC pc;
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
  //#ifdef USING_MKL
  //ierr = PCFactorSetMatSolverPackage(pc, MATSOLVERMKL_PARDISO); CHKERRQ(ierr);
  //#else
  ierr = PCFactorSetMatSolverPackage(pc, MATSOLVERMUMPS); CHKERRQ(ierr);
  //#endif

  EPSView(eps,PETSC_VIEWER_STDOUT_WORLD);

  /*****************************
   *  SOLVE
   *****************************/
  PetscLogDouble t1,t2;
  ierr = PetscTime(&t1);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscTime(&t2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Elapsed Time in EPSSolve: %f\n",t2-t1);CHKERRQ(ierr);
  PetscInt nconv;
  EPSGetConverged(eps, &nconv);
  PetscScalar kr, ki;
  Vec xr, xi;
  MatCreateVecs(fAmat,NULL,&xr);
  MatCreateVecs(fAmat,NULL,&xi);
  for (int i = 0; i < nconv; ++i) {
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error);
    std::cout<<kr<<std::endl;
    std::cout<<error<<std::endl;
  }

  EPSDestroy(&eps);
  ierr=MatDestroy(&fAmat);CHKERRQ(ierr);
  ierr=MatDestroy(&fBmat);CHKERRQ(ierr);
  if(fIsParallelStorage){
    ierr = PetscFree(ibP);CHKERRQ(ierr);
    ierr = PetscFree(jbP);CHKERRQ(ierr);
  }
  VecDestroy(&xr);
  VecDestroy(&xi);

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
