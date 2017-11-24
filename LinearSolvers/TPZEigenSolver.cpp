//
// Created by Francisco Teixeira Orlandini on 11/23/17.
//

#include "TPZEigenSolver.h"
template<typename TVar>
bool TPZEigenSolver<TVar>::IsGeneralised() const {
    return fIsGeneralised;
}
template<typename TVar>
TPZAutoPointer<TPZMatrix<TVar>> TPZEigenSolver<TVar>::MatrixA() {
    return fMatrixA;
}

template<typename TVar>
TPZAutoPointer<TPZMatrix<TVar>> TPZEigenSolver<TVar>::MatrixB() {
    return fMatrixB;
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetMatrixA(TPZAutoPointer<TPZMatrix<STATE>> mat) {
    fMatrixA = mat;
}

template<typename TVar>
void TPZEigenSolver<TVar>::SetMatrixB(TPZAutoPointer<TPZMatrix<STATE>> mat) {
    if(!fIsGeneralised){
        DebugStop();//Why are you setting rhs Matrix if is not a generalised eigenvalue problem?
    }
    fMatrixB = mat;
}
