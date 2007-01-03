/* Generated by Together */

#include "pzmgsolver.h"
#include "pztransfer.h"
using namespace std;


TPZMGSolver::TPZMGSolver(TPZAutoPointer<TPZTransfer> trf,const  TPZMatrixSolver &sol, int nvar, TPZAutoPointer<TPZMatrix> refmat) : TPZMatrixSolver(refmat),fStep(trf) {
  fCoarse = (TPZMatrixSolver *) sol.Clone();
  //  fTransfer = new TPZMatrixSolver::TPZContainer(trf);
  fNVar = nvar;
}

TPZMGSolver::TPZMGSolver(TPZAutoPointer<TPZTransfer> trf,const  TPZMatrixSolver &sol, int nvar) : TPZMatrixSolver(),fStep(trf) {
  fCoarse = (TPZMatrixSolver *) sol.Clone();
  //  fTransfer = new TPZMatrixSolver::TPZContainer(trf);
  fNVar = nvar;
}

void TPZMGSolver::Solve(const TPZFMatrix &F, TPZFMatrix &result, TPZFMatrix *residual){
  if(!Matrix() || !TransferMatrix()) {
    cout << "TPZMGSolver::Solve called without a matrix pointer\n";
    exit(-1);
  }
  TPZAutoPointer<TPZMatrix> mat = Matrix();
  if(result.Rows() != mat->Rows() || result.Cols() != F.Cols()) {
    result.Redim(mat->Rows(),F.Cols());
  }

  TPZFMatrix FCoarse,UCoarse;
  TPZAutoPointer<TPZTransfer> tr = TransferMatrix();
  tr->TransferResidual(F,FCoarse);
  fCoarse->Solve(FCoarse,UCoarse);
  tr->TransferSolution(UCoarse,result);
  if(residual) Matrix()->Residual(F,result,*residual);
}

TPZMGSolver::TPZMGSolver(const TPZMGSolver & copy): TPZMatrixSolver(copy), fStep(copy.fStep) {
    fCoarse = (TPZMatrixSolver *) copy.fCoarse->Clone();
    fNVar = copy.fNVar;
    //fTransfer = copy.fTransfer;
    //    fTransfer->IncreaseRefCount();
}

TPZSolver * TPZMGSolver::Clone() const {
    return new TPZMGSolver(*this);
}

TPZMGSolver::~TPZMGSolver(){
    delete fCoarse;
    //    fTransfer->DecreaseRefCount();

}
void TPZMGSolver::ResetTransferMatrix(){
  //  fTransfer->SetMatrix(0);
  TPZAutoPointer<TPZTransfer> reset;
  fStep = reset;
}
void TPZMGSolver::SetTransferMatrix(TPZAutoPointer<TPZTransfer> Refmat){
  fStep = Refmat;
  //    if(fTransfer->Matrix() == Refmat || !fTransfer->Matrix()) {
  //        fTransfer->SetMatrix(Refmat);
  //    } else {
  //        fTransfer->DecreaseRefCount();
  //        fTransfer = new TPZContainer(Refmat);
  //    }
}
