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
	rhs.Redim(fMesh->NEquations(),1);//size of the original system (no filter)
	Assemble(*mat,rhs,guiInterface);
    return mat;
}

TPZMatrix<STATE> * TPZSBandStructMatrix::Create(){
	long neq = fEquationFilter.NActiveEquations();
	
	long band = fMesh->BandWidth();
	return new TPZSBMatrix<STATE>(neq,band);
}

TPZSBandStructMatrix::TPZSBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}
