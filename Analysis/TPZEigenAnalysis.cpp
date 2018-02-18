//
//
// Created by Francisco Teixeira Orlandini on 11/17/17.

#include "TPZEigenAnalysis.h"
#include "TPZEigenSolver.h"
#include "pzmaterial.h"

TPZEigenAnalysis::TPZEigenAnalysis() : TPZAnalysis() , fSolver(0){

}

TPZEigenAnalysis::TPZEigenAnalysis(TPZCompMesh *mesh, bool mustOptimizeBandwidth, std::ostream &out)
        :TPZAnalysis(mesh, mustOptimizeBandwidth,out), fSolver(0) {

}

TPZEigenAnalysis::TPZEigenAnalysis(TPZAutoPointer<TPZCompMesh> mesh, bool mustOptimizeBandwidth, std::ostream &out)
        : TPZAnalysis(mesh, mustOptimizeBandwidth, out) ,  fSolver(0){

}

TPZEigenSolver<STATE> &TPZEigenAnalysis::Solver() const {
    return *fSolver;
}

void TPZEigenAnalysis::SetSolver(TPZEigenSolver<STATE> &solver) {
    if(fSolver) delete fSolver;
    fSolver = (TPZEigenSolver<STATE> *) solver.Clone();
}

void TPZEigenAnalysis::Assemble()
{
    if(!fCompMesh || !fStructMatrix || !fSolver)
    {
        std::stringstream sout;
        sout << "TPZAnalysis::Assemble lacking definition for Assemble\n fCompMesh "<< (void *) fCompMesh
             << " fStructMatrix " << (void *) fStructMatrix.operator->()
             << " fSolver " << (void *) fSolver;
        return;
    }
    int numloadcases = ComputeNumberofLoadCases();
    long sz = fCompMesh->NEquations();
    //@TODO: Get rid of this fRhs
    fRhs.Redim(sz,numloadcases);

    if(fSolver->IsGeneralised()){
        std::map<int,TPZMaterial *> &materialVec = fCompMesh->MaterialVec();
        for (auto &&item : materialVec) {
            TPZMaterial * mat = (item.second);
            mat->SetMatrixA();
        }
    }
    if(fSolver->MatrixA() && fSolver->MatrixA()->Rows()==sz){
        fSolver->MatrixA()->Zero();
        fStructMatrix->Assemble(*(fSolver->MatrixA().operator ->()),fRhs,fGuiInterface);
    }
    else{
        TPZAutoPointer<TPZMatrix<STATE> >mat(fStructMatrix->CreateAssemble(fRhs,fGuiInterface));
        fSolver->SetMatrixA(mat);
    }
    //fSolver->UpdateFrom(fSolver->MatrixA());
    if(fSolver->IsGeneralised()) {
        std::map<int, TPZMaterial *> &materialVec = fCompMesh->MaterialVec();
        for (auto &&item : materialVec) {
            TPZMaterial *mat = (item.second);
            mat->SetMatrixB();
        }
        if (fSolver->MatrixB() && fSolver->MatrixB()->Rows() == sz) {
            fSolver->MatrixB()->Zero();
            fStructMatrix->Assemble(*(fSolver->MatrixB().operator->()), fRhs, fGuiInterface);
        } else {
            TPZAutoPointer<TPZMatrix<STATE> >mat(fStructMatrix->CreateAssemble(fRhs, fGuiInterface));
            fSolver->SetMatrixB(mat);
        }
    }
}

template <class TVar>
void TPZEigenAnalysis::TransferEigenVector(const TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &originalSol,
                                           TPZFMatrix<TVar> &newSol, const unsigned int &i, const bool isAbsVal) {
    originalSol.GetSub(0,i,originalSol.Rows(),1,newSol);
}

template<>
void TPZEigenAnalysis::TransferEigenVector<float>(const TPZFMatrix<typename SPZAlwaysComplex<float>::type> &originalSol,
                                                  TPZFMatrix<float> &newSol, const unsigned int &i, const bool isAbsVal) {
    for (int j = 0; j < originalSol.Rows(); ++j) {
        newSol(j,0) = isAbsVal ? std::abs(originalSol.GetVal(j,i)) : std::real(originalSol.GetVal(j,i));
    }
}

template<>
void TPZEigenAnalysis::TransferEigenVector<double>(const TPZFMatrix<typename SPZAlwaysComplex<double>::type> &originalSol,
                                                   TPZFMatrix<double> &newSol, const unsigned int &i, const bool isAbsVal) {
    for (int j = 0; j < originalSol.Rows(); ++j) {
        newSol(j,0) = isAbsVal ? std::abs(originalSol.GetVal(j,i)) : std::real(originalSol.GetVal(j,i));
    }
}

void TPZEigenAnalysis::Solve() {
    long numeq = fCompMesh->NEquations();
    long nReducedEq = fStructMatrix->NReducedEquations();
    if (nReducedEq == numeq)
    {
        fEigenvalues.Resize(numeq);
        fEigenvectors.Redim(numeq,1);
        fSolver->Solve(fEigenvalues, fEigenvectors);
        fSolution.Resize(numeq, 1);
        for (int i = 0; i < fEigenvalues.size(); ++i) {
            TransferEigenVector(fEigenvectors,fSolution,i, fSolver->IsAbsoluteValue());
        }
    }
    else
    {
        fEigenvalues.Resize(nReducedEq);
        fEigenvectors.Redim(nReducedEq,1);
        //@TODO:THINK ABOUT EQUATION FILTER
        //fStructMatrix->EquationFilter().Gather(delu,fSolution);
        fSolver->Solve(fEigenvalues, fEigenvectors);
        TPZFMatrix<STATE> sol(nReducedEq,1);
        for (int i = 0; i < fEigenvalues.size(); ++i) {
            TransferEigenVector(fEigenvectors,fSolution,i, fSolver->IsAbsoluteValue());
        }
        fStructMatrix->EquationFilter().Scatter(sol,fSolution);
    }
    fCompMesh->LoadSolution(fSolution);
    fCompMesh->TransferMultiphysicsSolution();

}

TPZFMatrix<typename SPZAlwaysComplex<STATE>::type>
TPZEigenAnalysis::GetEigenvectors() const {
    return fEigenvectors;
}

TPZManVector<typename SPZAlwaysComplex<STATE>::type>
TPZEigenAnalysis::GetEigenvalues() const {
    return fEigenvalues;
}
