
#include <pzsespmat.h>
#include <pzsysmp.h>
#include "petscmat.h"
static char help[] = "Testing MatCreateSeqSBAIJWithArrays().\n\n";
int main(int argc, char **argv){
  const int nRows = 3;
  const int nCols = nRows;
  TPZSYsmpMatrix<STATE> matrix(nRows,nCols);
  TPZVec<long> ia ({0,2,3,4});
  TPZVec<long> ja({0,2,1,2});

  TPZVec<STATE> aVec({1.,3.,4.,5.});
  matrix.SetData(ia,ja,aVec);
  matrix.Print("mat");
  const int blockSize = 1;
  Mat A;
  //@TODO:Think on int/long issue
  //@TODO:Think on PetsScalar issue
  TPZVec<int> iaP ({0,2,3,4});
  TPZVec<int> jaP({0,2,1,2});
  TPZVec<std::complex<double>> aVecP({1.,3.,4.,5.});

  PetscInitialize(&argc,&argv,(char *)0,help);
  MatCreateSeqSBAIJWithArrays(MPI_COMM_WORLD,blockSize,nRows,nCols,(PetscInt *)iaP.begin(),(PetscInt *)jaP.begin(),(PetscScalar *)aVecP.begin(),&A);
  MatView(A, PETSC_VIEWER_STDOUT_WORLD);
  PetscViewerFlush(PETSC_VIEWER_STDOUT_WORLD); // Diag: 1, 2, 1
  return 0;
}
