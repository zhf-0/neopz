/**
 * @file
 * @brief Contains the implementation of the TPZSBandStructMatrix methods. 
 */

#include "pzsbstrmatrix.h"

#include "pzsbndmat.h"
#include "pzcmesh.h"

TPZStructMatrix * TPZSBandStructMatrix::Clone(){
    return new TPZSBandStructMatrix(*this);
}

TPZMatrix<STATE> * TPZSBandStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	TPZMatrix<STATE> *mat = Create();
	rhs.Redim(mat->Rows(),1);
	Assemble(*mat,rhs,guiInterface);
    return mat;
}

TPZMatrix<STATE> * TPZSBandStructMatrix::Create(){
    if (fEquationFilter.IsActive()) {
        DebugStop();
    }
	long neq = fEquationFilter.NActiveEquations();
	
	long band = fMesh->BandWidth();
	return new TPZSBMatrix<STATE>(neq,band);
}

TPZSBandStructMatrix::TPZSBandStructMatrix(TPZCompMesh *mesh) : TPZRegisterClassId(&TPZSBandStructMatrix::ClassId),TPZStructMatrix(mesh)
{
}

int TPZSBandStructMatrix::ClassId() {
	//CLASSIDFRANreturn TPZStructMatrix::ClassId()^Hash("TPZBandStructMatrix");
	return 666;
}

TPZSBandStructMatrix::TPZSBandStructMatrix() : TPZRegisterClassId(&TPZSBandStructMatrix::ClassId),TPZStructMatrix()
{
}
