//
// Created by Francisco Teixeira Orlandini on 12/13/17.
//
#include "TPZSlepcEPSHandler.h"
#ifdef USING_SLEPC

#include "TPZSlepcSTHandler.h"
#include <petsctime.h>
#include <petscksp.h>

/*******************
*    GENERAL       *
*******************/
template<class TVar>
TPZSlepcEPSHandler<TVar>::TPZSlepcEPSHandler() : fVerbose(true){
  EPSCreate( PETSC_COMM_WORLD, &fEps );
}
template<class TVar>
TPZSlepcEPSHandler<TVar>::~TPZSlepcEPSHandler() {
  EPSDestroy(&fEps);
}

template<class TVar>
int TPZSlepcEPSHandler<TVar>::SolveEigenProblem(TPZMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors){

  TPZFYsmpMatrix<TVar> *Afsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&A);
  if(Afsparse){
    return SolveEigenProblem((TPZFYsmpMatrix<TVar>&)A,w,eigenVectors);
  }
  else{
    std::cout<<"TPZSlepcEPSHandler does not support this matrix format"<<std::endl;
    DebugStop();
  }
  return 0;
}
template<class TVar>
int TPZSlepcEPSHandler<TVar>::SolveEigenProblem(TPZMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  TPZFYsmpMatrix<TVar> *Afsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&A);
  if(Afsparse){
    return SolveEigenProblem((TPZFYsmpMatrix<TVar>&)A,w);
  }
  else{
    std::cout<<"TPZSlepcEPSHandler does not support this matrix format"<<std::endl;
    DebugStop();
  }
  return 0;
}
template<class TVar>
int TPZSlepcEPSHandler<TVar>::SolveGeneralisedEigenProblem(TPZMatrix<TVar> &A, TPZMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
  TPZFYsmpMatrix<TVar> *Afsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&A);
  TPZFYsmpMatrix<TVar> *Bfsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&B);
  if(Afsparse && Bfsparse){
    return SolveGeneralisedEigenProblem((TPZFYsmpMatrix<TVar>&)A,(TPZFYsmpMatrix<TVar>&)B,w,eigenVectors);
  }
  else{
    std::cout<<"TPZSlepcEPSHandler does not support this matrix format"<<std::endl;
    DebugStop();
  }
  return 0;
}
template<class TVar>
int TPZSlepcEPSHandler<TVar>::SolveGeneralisedEigenProblem(TPZMatrix<TVar> &A, TPZMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  TPZFYsmpMatrix<TVar> *Afsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&A);
  TPZFYsmpMatrix<TVar> *Bfsparse = dynamic_cast<TPZFYsmpMatrix<TVar>* >(&B);
  if(Afsparse && Bfsparse){
    return SolveGeneralisedEigenProblem((TPZFYsmpMatrix<TVar>&)A,(TPZFYsmpMatrix<TVar>&)B,w);
  }
  else{
    std::cout<<"TPZSlepcEPSHandler does not support this matrix format"<<std::endl;
    DebugStop();
  }
  return 0;
}


/*******************
*    TPZFYSMPMATRIX    *
*******************/
template<class TVar>
int TPZSlepcEPSHandler<TVar>::SolveGeneralisedEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZFYsmpMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
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

  EPSSetOperators(fEps, fAmat, fBmat);

  if(fVerbose){
    EPSView(fEps,PETSC_VIEWER_STDOUT_WORLD);
    PetscViewerAndFormat *vf;
    PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
    EPSMonitorSet(fEps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,
                                           PetscReal*,PetscInt,void*))EPSMonitorFirst,vf,
                  (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
  }
  /*****************************
   *  SOLVE
   *****************************/
  PetscLogDouble t1,t2;
  ierr = PetscTime(&t1);CHKERRQ(ierr);
  ierr = EPSSolve(fEps);CHKERRQ(ierr);
  ierr = PetscTime(&t2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Elapsed Time in EPSSolve: %f\n",t2-t1);CHKERRQ(ierr);
  /*
     Optional: Get some information from the solver and display it
  */
  PetscInt its,lits,maxit,nev;
  PetscReal tol;
  EPSType type;
  KSP ksp;
  ST st;
  ierr = EPSGetIterationNumber(fEps, &its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetST(fEps,&st);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPGetTotalIterations(ksp, &lits);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %D\n",lits);CHKERRQ(ierr);
  ierr = EPSGetType(fEps, &type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(fEps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(fEps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Show detailed info unless -terse option is given by user
   */
  if (fVerbose) {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSReasonView(fEps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = EPSErrorView(fEps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  } else {
    ierr = EPSErrorView(fEps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  }
  

  return 1;
}

template<class TVar>
int TPZSlepcEPSHandler<TVar>::SolveGeneralisedEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZFYsmpMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  DebugStop();
  return 0;
}

template<class TVar>
int TPZSlepcEPSHandler<TVar>::SolveEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
  DebugStop();
  return 0;
}

template<class TVar>
int TPZSlepcEPSHandler<TVar>::SolveEigenProblem(TPZFYsmpMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
  DebugStop();
  return 0;
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetST(const TPZSlepcSTHandler &st){
  PetscErrorCode ierr = EPSSetST(fEps,st.fSt);
  if(ierr) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::GetST(TPZSlepcSTHandler *st){
  PetscErrorCode ierr = EPSGetST(fEps,&(st->fSt));
  if(ierr) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetProblemType(const EPSProblemType eps_problem)
{
  const PetscErrorCode ierr = EPSSetProblemType (fEps, eps_problem);
  if(ierr!=0){
    PZError<<"Invalid problem type.\n";
    PZError<<"Valid problems are:\n";
    PZError<<"EPS_HEP, EPS_GHEP, EPS_NHEP, EPS_GNHEP, EPS_PGNHEP or EPS_GHIEP\n";
    DebugStop();
  }
  this->SetAsGeneralised(!(eps_problem == EPS_HEP || eps_problem == EPS_NHEP));
}

template<class TVar>
EPSProblemType TPZSlepcEPSHandler<TVar>::GetProblemType() const
{
  EPSProblemType eps_problem;
  const PetscErrorCode ierr = EPSGetProblemType (fEps, &eps_problem);
  if(ierr!=0) DebugStop();
  return eps_problem;
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetTargetEigenvalue (const PetscScalar &target){
  // set target eigenvalues to solve for
  // in all transformation except STSHIFT there is a direct connection between
  // the target and the shift, read more on p41 of SLEPc manual.
  const PetscErrorCode ierr = EPSSetTarget (fEps, target );
  if(ierr!=0) DebugStop();
}

template<class TVar>
PetscScalar TPZSlepcEPSHandler<TVar>::GetTargetEigenvalue () const {
  // set target eigenvalues to solve for
  // in all transformation except STSHIFT there is a direct connection between
  // the target and the shift, read more on p41 of SLEPc manual.
  PetscScalar target;
  const PetscErrorCode ierr = EPSGetTarget (fEps, &target );
  if(ierr!=0) DebugStop();
  return target;
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetWhichEigenpairs (const EPSWhich which){
  // set which portion of the eigenspectrum to solve for
  const PetscErrorCode ierr = EPSSetWhichEigenpairs (fEps, which);
  if(ierr!=0) DebugStop();
}

template<class TVar>
EPSWhich TPZSlepcEPSHandler<TVar>::GetWhichEigenpairs () const {
  // set which portion of the eigenspectrum to solve for
  EPSWhich which;
  const PetscErrorCode ierr = EPSGetWhichEigenpairs (fEps, &which);
  if(ierr!=0) DebugStop();
  return which;
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetKrylovOptions (const bool &pLocking, const PetscReal &restart){
  EPSType currentType;
  EPSGetType(fEps, &currentType);
  if(strcmp(currentType,EPSKRYLOVSCHUR)){
    PZError<<"EPSType is not EPSKRYLOVSCHUR\n";
    DebugStop();
  }

  PetscBool locking = pLocking ? PETSC_TRUE : PETSC_FALSE;

  PetscErrorCode ierr = EPSKrylovSchurSetLocking (fEps, locking);
  if(ierr!=0) DebugStop();
  ierr = EPSKrylovSchurSetRestart(fEps, restart);
  if(ierr!=0) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::GetKrylovOptions (bool *pLocking, PetscReal *restart) const{
  EPSType currentType;
  EPSGetType(fEps, &currentType);
  if(strcmp(currentType,EPSKRYLOVSCHUR)){
    PZError<<"EPSType is not EPSKRYLOVSCHUR\n";
    DebugStop();
  }

  PetscBool locking;


  PetscErrorCode ierr = EPSKrylovSchurGetLocking(fEps, &locking);
  if(ierr!=0) DebugStop();
  *pLocking = locking;
  ierr = EPSKrylovSchurGetRestart(fEps, restart);
  if(ierr!=0) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetTolerances(const PetscReal &tol, const PetscInt &max_its) {
  const PetscErrorCode ierr = EPSSetTolerances(fEps,tol,max_its);
  if(ierr != 0) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::GetTolerances(PetscReal *tol, PetscInt *max_its) const {
  const PetscErrorCode ierr = EPSGetTolerances(fEps,tol,max_its);
  if(ierr != 0) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetConvergenceTest(const EPSConv &test) {
  const PetscErrorCode ierr = EPSSetConvergenceTest(fEps,test);
  if(ierr != 0) DebugStop();
}

template<class TVar>
EPSConv TPZSlepcEPSHandler<TVar>::GetConvergenceTest() const {
  EPSConv test;
  const PetscErrorCode ierr = EPSGetConvergenceTest(fEps,&test);
  if(ierr != 0) DebugStop();
  return test;
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetTrueResidual(const bool &pOpt) {
  ST st;
  EPSGetST(fEps, &st);
  STType type;
  STGetType(st, &type);

  if(strcmp(type,STSINVERT)){
    PZError<<__PRETTY_FUNCTION__<<"is only available if STTYpe is STSINVERT\n";
    DebugStop();
  }
  const PetscBool opt= pOpt ? PETSC_TRUE : PETSC_FALSE;
  const PetscErrorCode ierr = EPSSetTrueResidual(fEps,opt);
  if(ierr != 0) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetType(const EPSType &type) {
  const PetscErrorCode ierr = EPSSetType(fEps,type);
  if(ierr != 0) DebugStop();
}

template<class TVar>
EPSType TPZSlepcEPSHandler<TVar>::GetType() const {
  EPSType type;
  const PetscErrorCode ierr = EPSGetType(fEps,&type);
  if(ierr != 0) DebugStop();
  return type;
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetEPSDimensions(const PetscInt &nev, const PetscInt &ncv, const PetscInt &mpd) {
  const PetscErrorCode ierr = EPSSetDimensions(fEps,nev,ncv,mpd);
  if(ierr != 0) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::GetEPSDimensions(PetscInt *nev, PetscInt *ncv, PetscInt *mpd) const {
  const PetscErrorCode ierr = EPSGetDimensions(fEps,nev,ncv,mpd);
  if(ierr != 0) DebugStop();
}

template<class TVar>
void TPZSlepcEPSHandler<TVar>::SetVerbose(bool fVerbose) {
  TPZSlepcEPSHandler::fVerbose = fVerbose;
}

template class TPZSlepcEPSHandler< std::complex<float> >;

template class TPZSlepcEPSHandler< std::complex<double> >;
template class TPZSlepcEPSHandler< std::complex<long double> >;
template class TPZSlepcEPSHandler<float >;

template class TPZSlepcEPSHandler<double >;
template class TPZSlepcEPSHandler<long double>;

#endif
