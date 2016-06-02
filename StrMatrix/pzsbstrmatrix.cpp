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

TPZMatrix<STATE> * TPZSBandStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs){
	TPZMatrix<STATE> *mat = Create();
	rhs.Redim(mat->Rows(),1);
	Assemble(*mat,rhs);
    return mat;
}

TPZMatrix<STATE> * TPZSBandStructMatrix::Create(){
    if (fAssembleConfig.fEquationFilter.IsActive()) {
        DebugStop();
    }
	long neq = fAssembleConfig.fEquationFilter.NActiveEquations();
	
	long band = fAssembleConfig.fMesh->BandWidth();
	return new TPZSBMatrix<STATE>(neq,band);
}

TPZSBandStructMatrix::TPZSBandStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}
