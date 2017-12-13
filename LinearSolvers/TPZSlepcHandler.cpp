//
// Created by Francisco Teixeira Orlandini on 12/13/17.
//
#include "TPZSlepcHandler.h"

/*******************
*    TPZFYSMPMATRIX    *
*******************/
template<class TVar>
int TPZSlepcHandler<TVar>::SolveEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors){
  DebugStop();
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  DebugStop();
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveGeneralisedEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZFYsmpMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
  DebugStop();
  return 0;
}
template<class TVar>
int TPZSlepcHandler<TVar>::SolveGeneralisedEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZFYsmpMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  DebugStop();
  return 0;
}

#ifdef USING_SLEPC

#include <slepceps.h>
#include <ldap.h>
#include <petsctime.h>
#include "slepcrg.h"
/*******************
*    TPZFYSMPMATRIX    *
*******************/
template<>
int TPZSlepcHandler<STATE>::SolveGeneralisedEigenProblem(TPZFYsmpMatrix<STATE> &A, TPZFYsmpMatrix< STATE > &B , TPZVec < typename SPZAlwaysComplex<STATE>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<STATE>::type > &eigenVectors){
  /*****************************
   *  INITIALIZE STRUCTURES
   *****************************/

  EPS eps;

  bool parallelStructures = true;
  if(parallelStructures){
    char *var;
    setenv("OMP_NUM_THREADS", "8", 1);
    var = getenv("OMP_NUM_THREADS");
    setenv("MPC_NUM_THREADS", "8", 1);
    std::cout<<var<<std::endl;
  }
  SlepcInitialize((int *)0, (char ***)0, (const char*)0,(const char*)0 );

  EPSCreate( PETSC_COMM_WORLD, &eps);
  /*****************************
   *  CREATE MATRICES
   *****************************/
  Mat Amat, Bmat;
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

  PetscErrorCode ierr;
  PetscInt *iaP, *jaP;
  std::cout<<"Creating PETSc Amat...";
  ierr = PetscMalloc1(A.fIA.size(),&iaP);CHKERRQ(ierr);
  ierr = PetscMalloc1(A.fJA.size(),&jaP);CHKERRQ(ierr);

  for (int j = 0; j < A.fIA.size(); ++j) {
    iaP[j]=A.fIA[j];
  }
  for (int j = 0; j < A.fJA.size(); ++j) {
    jaP[j]=A.fJA[j];
  }
  if (parallelStructures)
  {
    ierr = MatCreate(PETSC_COMM_WORLD,&Amat);CHKERRQ(ierr);
    ierr = MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,nRows,nCols);CHKERRQ(ierr);
    ierr = MatSetType(Amat,MATMPIAIJ);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocationCSR(Amat,iaP,jaP,(PetscScalar *)A.fA.begin());CHKERRQ(ierr);
    ierr = PetscFree(iaP);CHKERRQ(ierr);
    ierr = PetscFree(jaP);CHKERRQ(ierr);
  }
  else
  {
    error=MatCreateSeqBAIJWithArrays(MPI_COMM_WORLD,blockSize,nRows,nCols,iaP,jaP,(PetscScalar *)A.fA.begin(),&Amat);
  }
  std::cout<<"Created!"<<std::endl;
  std::cout<<"Creating PETSc Bmat...";
  PetscInt *ibP, *jbP;
  ierr = PetscMalloc1(B.fIA.size(),&ibP);CHKERRQ(ierr);
  ierr = PetscMalloc1(B.fJA.size(),&jbP);CHKERRQ(ierr);
  for (int j = 0; j < B.fIA.size(); ++j) {
    ibP[j]=B.fIA[j];
  }
  for (int j = 0; j < B.fJA.size(); ++j) {
    jbP[j]=B.fJA[j];
  }
  if(parallelStructures)
  {
    ierr = MatCreate(PETSC_COMM_WORLD,&Bmat);CHKERRQ(ierr);
    ierr = MatSetSizes(Bmat,PETSC_DECIDE,PETSC_DECIDE,nRows,nCols);CHKERRQ(ierr);
    ierr = MatSetType(Bmat,MATMPIAIJ);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocationCSR(Bmat,ibP,jbP,(PetscScalar *)B.fA.begin());CHKERRQ(ierr);
    ierr = PetscFree(ibP);CHKERRQ(ierr);
    ierr = PetscFree(jbP);CHKERRQ(ierr);
  }
  else
  {
    ierr=MatCreateSeqBAIJWithArrays(MPI_COMM_WORLD,blockSize,nRows,nCols,ibP,jbP,(PetscScalar *)B.fA.begin(),&Bmat);CHKERRQ(ierr);
  }
  std::cout<<"Created!"<<std::endl;

  /*****************************
   *  DEFINE PROBLEM TYPE
   *****************************/
  EPSSetOperators(eps, Amat, Bmat);
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
  ierr = PetscPrintf(PETSC_COMM_WORLD," Elapsed Time in EPSSolve: %f\n",t2-t1);
  PetscInt nconv;
  EPSGetConverged(eps, &nconv);
  PetscScalar kr, ki;
  Vec xr, xi;
  MatCreateVecs(Amat,NULL,&xr);
  MatCreateVecs(Amat,NULL,&xi);
  for (int i = 0; i < nconv; ++i) {
    EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
    EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error);
    std::cout<<kr<<std::endl;
    std::cout<<error<<std::endl;
  }
  EPSDestroy(&eps);
  MatDestroy(&Amat);
  MatDestroy(&Bmat);
  VecDestroy(&xr);
  VecDestroy(&xi);
  if(!parallelStructures){
    ierr = PetscFree(ibP);CHKERRQ(ierr);
    ierr = PetscFree(jbP);CHKERRQ(ierr);
  }
  return 1;
}
template<>
int TPZSlepcHandler<STATE>::SolveGeneralisedEigenProblem(TPZFYsmpMatrix<STATE> &A, TPZFYsmpMatrix< STATE > &B , TPZVec < typename SPZAlwaysComplex<STATE>::type > &w){
  DebugStop();
  return 0;
}



#endif

template class TPZSlepcHandler< std::complex<float> >;
template class TPZSlepcHandler< std::complex<double> >;
template class TPZSlepcHandler< std::complex<long double> >;

template class TPZSlepcHandler<long >;

template class TPZSlepcHandler<float >;
template class TPZSlepcHandler<double >;
template class TPZSlepcHandler<long double>;

template class TPZSlepcHandler<int >;
template class TPZSlepcHandler<TPZFlopCounter>;
