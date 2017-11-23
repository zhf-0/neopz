//
// Created by Francisco Teixeira Orlandini on 11/17/17.
//

#include "TPZEigenAnalysis.h"

const TPZAutoPointer<TPZStructMatrix> &TPZEigenAnalysis::GetAMatrix() const {
    return fAStructMatrix;
}

void TPZEigenAnalysis::SetAMatrix(const TPZAutoPointer<TPZStructMatrix> &fAStructMatrix) {
    TPZEigenAnalysis::fAStructMatrix = fAStructMatrix;
}

const TPZAutoPointer<TPZStructMatrix> &TPZEigenAnalysis::GetBMatrix() const {
    return fBStructMatrix;
}

void TPZEigenAnalysis::SetBMatrix(const TPZAutoPointer<TPZStructMatrix> &fBStructMatrix) {
    TPZEigenAnalysis::fBStructMatrix = fBStructMatrix;
}

TPZEigenSolver<STATE> * TPZEigenAnalysis::GetSolver() const {
    return fSolver;
}

void TPZEigenAnalysis::SetSolver(TPZEigenSolver<STATE> * &solver) {
    fSolver = solver;
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
