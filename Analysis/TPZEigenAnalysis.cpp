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
        sout << "TPZAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh
             << " fStructMatrix " << (void *) fStructMatrix.operator->()
             << " fSolver " << (void *) fSolver;
#ifndef WINDOWS
        sout << " at file " << __FILE__ << " line " << __LINE__ ;
#else
        sout << " TPZAnalysis::Assemble() " ;
#endif
#ifdef LOG4CXX
        LOGPZ_ERROR(logger,sout.str().c_str());
#else
        std::cout << sout.str().c_str() << std::endl;
#endif
        return;
    }
    int numloadcases = ComputeNumberofLoadCases();
    long sz = fCompMesh->NEquations();
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

void TPZEigenAnalysis::Solve() {
    long numeq = fCompMesh->NEquations();
    long nReducedEq = fStructMatrix->NReducedEquations();
    if (nReducedEq == numeq)
    {
//        @TODO: Fix this call
        fEigenvalues.Resize(numeq);
        fEigenvectors.Redim(numeq,1);
        fSolver->Solve(fEigenvalues, fEigenvectors);
        fSolution.Resize(numeq, 1);
        TransferEigenVector(fEigenvectors,fSolution,0, fSolver->IsAbsoluteValue());
    }
    else
    {
        fEigenvalues.Resize(nReducedEq);
        fEigenvectors.Redim(nReducedEq,1);
        //@TODO:THINK ABOUT EQUATION FILTER
        //fStructMatrix->EquationFilter().Gather(delu,fSolution);
        fSolver->Solve(fEigenvalues, fEigenvectors);
        //@TODO:Must now select the eigenpairs based on fSolver->fHowManyEigenvalues, fSolver->fDesiredPartOfTheSpectrum...
        TPZFMatrix<STATE> sol(nReducedEq,1);
        TransferEigenVector(fEigenvectors,sol,0, fSolver->IsAbsoluteValue());
        //For now, the first solution is loaded.
        fStructMatrix->EquationFilter().Scatter(sol,fSolution);
    }
#ifdef LOG4CXX
    std::stringstream sout;
    TPZStepSolver<STATE> *step = dynamic_cast<TPZStepSolver<STATE> *> (fSolver);
    if(!step) DebugStop();
    long nsing = step->Singular().size();
	if(nsing && logger->isWarnEnabled()) {
		sout << "Number of singular equations " << nsing;
		std::list<long>::iterator it = step->Singular().begin();
		if(nsing) sout << "\nSingular modes ";
		while(it != step->Singular().end())
		{
			sout << *it << " ";
			it++;
		}
		if(nsing) sout << std::endl;
		LOGPZ_WARN(logger,sout.str())
	}
#endif
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Solution norm " << Norm(fSolution) << std::endl;
		fSolution.Print("delu",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
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

template <class TVar>
void TPZEigenAnalysis::TransferEigenVector(const TPZFMatrix<typename SPZAlwaysComplex<TVar>::type> &originalSol,
                                           TPZFMatrix<TVar> &newSol, const unsigned int &i, const bool isAbsVal) {
    originalSol.GetSub(i,0,originalSol.Cols(),1,newSol);

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
