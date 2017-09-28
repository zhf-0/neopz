
#include "TPZSkylineNSymStructMatrix.h"
#include "pzskylnsymmat.h"

TPZSkylineNSymStructMatrix::TPZSkylineNSymStructMatrix(TPZCompMesh *cmesh)
                           : TPZRegisterClassId(&TPZSkylineNSymStructMatrix::ClassId),TPZSkylineStructMatrix(cmesh)
{
  ///nothing here
}

TPZSkylineNSymStructMatrix::TPZSkylineNSymStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh)
                            : TPZRegisterClassId(&TPZSkylineNSymStructMatrix::ClassId),TPZSkylineStructMatrix(cmesh)
{
  ///nothing here
}


TPZSkylineNSymStructMatrix::TPZSkylineNSymStructMatrix(const TPZSkylineStructMatrix &cp):
        TPZRegisterClassId(&TPZSkylineNSymStructMatrix::ClassId),TPZSkylineStructMatrix(cp)
{
  ///nothing here
}

TPZSkylineNSymStructMatrix::~TPZSkylineNSymStructMatrix()
{
  ///nothing here
}

TPZMatrix<STATE> * TPZSkylineNSymStructMatrix::ReallyCreate(long neq, const TPZVec<long> &skyline)
{
  return new TPZSkylNSymMatrix<STATE>(neq,skyline);
}

int TPZSkylineNSymStructMatrix::ClassId() {
  //CLASSIDFRANreturn TPZSkylineStructMatrix::ClassId()^Hash("TPZSkylineNSymStructMatrix");
  return 666;
}

