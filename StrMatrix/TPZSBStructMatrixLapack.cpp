/**
 * @file
 * @brief Contains the implementation of the TPZSBStructMatrixLapack methods.
 */

#include "TPZSBStructMatrixLapack.h"
#include "TPZSBMatrixLapack.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzvec.h"

TPZStructMatrix * TPZSBStructMatrixLapack::Clone(){
    return new TPZSBStructMatrixLapack(*this);
}

TPZSBStructMatrixLapack::TPZSBStructMatrixLapack(const TPZSBStructMatrixLapack &cp)
:TPZStructMatrix(cp) {
}

TPZSBStructMatrixLapack::TPZSBStructMatrixLapack(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{
}

TPZSBStructMatrixLapack::TPZSBStructMatrixLapack(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrix(cmesh)
{
}

TPZMatrix<STATE> * TPZSBStructMatrixLapack::Create(){
  long band = fMesh->BandWidth();
  long neq = fEquationFilter.NActiveEquations();
  return this->ReallyCreate(neq,band);
}


TPZMatrix<STATE> * TPZSBStructMatrixLapack::ReallyCreate(long neq, long band){
  return new TPZSBMatrixLapack<STATE>(neq,band);
}




TPZSBStructMatrixLapack::~TPZSBStructMatrixLapack(){
}




