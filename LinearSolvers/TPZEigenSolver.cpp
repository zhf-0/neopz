//
//
// Created by Francisco Teixeira Orlandini on 11/23/17.
#include "TPZEigenSolver.h"
#include "TPZSlepcHandler.h"
#include "TPZLapackWrapper.h"
#include "pzsbndmat.h"
#include "pzysmp.h"

template<typename TVar>
TPZEigenSolver<TVar>::TPZEigenSolver() : fIsGeneralised(false), fMustCalculateEigenVectors(true),
                                         fHowManyEigenValues(1),fDesiredPartOfSpectrum(MNE),fSpecifiedValue(0.),
                                         fEig(NULL){
}

template<typename TVar>
TPZEigenSolver<TVar>::TPZEigenSolver(const TPZEigenSolver &copy) {
  fIsGeneralised = copy.fIsGeneralised;
  fMustCalculateEigenVectors = copy.fMustCalculateEigenVectors;
  fHowManyEigenValues = copy.fHowManyEigenValues;
  fDesiredPartOfSpectrum = copy.fDesiredPartOfSpectrum;
  fSpecifiedValue = copy.fSpecifiedValue;
  fEigenvalues = copy.fEigenvalues;
  fEigenvectors = copy.fEigenvectors;
  fMatrixA = copy.fMatrixA;
  fMatrixB = copy.fMatrixB;
  fEig = copy.fEig;
}

template <class TVar>
void TPZEigenSolver<TVar>::Solve(TPZVec<typename SPZAlwaysComplex<TVar>::type> &eigenValues, TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &eigenVectors){
  if(!this->MatrixA() || (!this->MatrixB() && this->IsGeneralised() )) {
    std::cout << "TPZEigenSolver::Solve called without a matrix pointer"<<std::endl;
    if(this->MatrixA()){
      std::cout<<"Missing B matrix."<<std::endl;
    }
    DebugStop();
  }

  TPZAutoPointer<TPZMatrix<TVar> > matA = this->MatrixA();
  if(eigenValues.size() != matA->Rows()) {
    eigenValues.Resize(matA->Rows());
  }
  if(this->IsGeneralised()){
    TPZAutoPointer<TPZMatrix<TVar> > matB = this->MatrixB();
    if(!matA->SolveGeneralisedEigenProblem(matB,eigenValues,eigenVectors,fEig.operator->())) DebugStop();
  }else{
    if(!matA->SolveEigenProblem(eigenValues,eigenVectors,fEig.operator->())) DebugStop();
  }
}
/////////////////////FALTA REVISAO//////////////////////

template<typename TVar>
TPZAutoPointer<TPZMatrix<TVar>> TPZEigenSolver<TVar>::MatrixA() {
  return fMatrixA;
}

template<typename TVar>
TPZAutoPointer<TPZMatrix<TVar>> TPZEigenSolver<TVar>::MatrixB() {
  return fMatrixB;
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetMatrixA(TPZAutoPointer<TPZMatrix<TVar>> mat) {
  fMatrixA = mat;
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetMatrixB(TPZAutoPointer<TPZMatrix<TVar>> mat) {
  if(!fIsGeneralised){
    DebugStop();//Why are you setting rhs Matrix if is not a generalised eigenvalue problem?
  }
  fMatrixB = mat;
}

template<typename TVar>
int TPZEigenSolver<TVar>::ClassId() const{
  return 666;//@TODO: Implementar corretamente!
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetAsGeneralised(bool isGeneralised) {
  fIsGeneralised = isGeneralised;
}

template<typename TVar>
bool TPZEigenSolver<TVar>::IsGeneralised() const {
  return fIsGeneralised;
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetAbsoluteValue(bool isAbsoluteValue){
  fShowAbsoluteValue = isAbsoluteValue;
}

template<typename TVar>
bool TPZEigenSolver<TVar>::IsAbsoluteValue(){
  return fShowAbsoluteValue;
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetHowManyEigenValues(int howManyEigenValues) {
  fHowManyEigenValues = howManyEigenValues;
}

template<typename TVar>
int TPZEigenSolver<TVar>::GetHowManyEigenvalues() const {
  return fHowManyEigenValues;
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetDesiredPartOfSpectrum(EDesiredEigen desiredPartOfSpectrum) {
  fDesiredPartOfSpectrum = desiredPartOfSpectrum;
}

template<typename TVar>
EDesiredEigen TPZEigenSolver<TVar>::GetDesiredPartOfSpectrum() const {
  return fDesiredPartOfSpectrum;
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetSpecifiedValue(typename SPZAlwaysComplex<TVar>::type specifiedValue) {
  fSpecifiedValue = specifiedValue;
}

template<typename TVar>
typename SPZAlwaysComplex<TVar>::type TPZEigenSolver<TVar>::GetSpecifiedValue() const {
  return fSpecifiedValue;
}

template class TPZEigenSolver< std::complex<float> >;
template class TPZEigenSolver< std::complex<double> >;
template class TPZEigenSolver< std::complex<long double> >;

template class TPZEigenSolver<float >;
template class TPZEigenSolver<double >;
template class TPZEigenSolver<long double>;
