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

TPZEigenSolver<STATE> * TPZEigenAnalysis::GetSolver() const {
    return fSolver;
}

void TPZEigenAnalysis::SetSolver(TPZEigenSolver<STATE> * &solver) {
    fSolver = solver;
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
    if(fRhs.Rows() != numeq )
    {
        DebugStop();
    }
    long nReducedEq = fStructMatrix->NReducedEquations();
    if (nReducedEq == numeq)
    {
//        TPZFMatrix<STATE> residual(fRhs);
//        TPZFMatrix<STATE> delu(numeq,1,0.);
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            TPZFMatrix<STATE> res2(fRhs);
            fSolver->Matrix()->Residual(fSolution,fRhs,res2);
            std::stringstream sout;
            sout << "Residual norm " << Norm(res2) << std::endl;
    //		res2.Print("Residual",sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
//        @TODO: Fix this call
//        fSolver->Solve(residual, delu);
//        fSolution = delu;
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            if(!fSolver->Matrix()->IsDecomposed())
            {
                TPZFMatrix<STATE> res2(fRhs);
                fSolver->Matrix()->Residual(delu,fRhs,res2);
                std::stringstream sout;
                sout << "Residual norm " << Norm(res2) << std::endl;
                //            res2.Print("Residual",sout);
                LOGPZ_DEBUG(logger,sout.str())
            }
        }
#endif

    }
    else
    {
//        TPZFMatrix<STATE> residual(nReducedEq,1,0.);
//        TPZFMatrix<STATE> delu(nReducedEq,1,0.);
//        fStructMatrix->EquationFilter().Gather(fRhs,residual);
//        @TODO: Fix this call
//        fSolver->Solve(residual, delu);
//        fSolution.Redim(numeq,1);
//        fStructMatrix->EquationFilter().Scatter(delu,fSolution);
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
