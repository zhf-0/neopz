/**
 * @file
 * @brief Contains the TPZSBStructMatrixLapack class which implements Symmetric/Hermitian Band Structural Matrices.
 */

#ifndef TPZSBSTRUCTMATRIXLAPACK_H
#define TPZSBSTRUCTMATRIXLAPACK_H

#include "pzstrmatrix.h"

class TPZCompMesh;
template<class TVar>
class TPZFMatrix;
template<class TVar>
class TPZMatrix;

/**
 * @brief Implements Symmetric/Hermitian Band Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSBStructMatrixLapack : public TPZStructMatrix {
protected:
    
    /** @brief the equations which should actually be assembled */
    TPZVec<long> fActiveEquations;
    
    /** @brief Equation destination */
    TPZVec<long> fEquationDestination;
    
    /** Returns the skyline matrix object */
    virtual TPZMatrix<STATE> * ReallyCreate(long neq, long band);
    
public:    
	
	TPZSBStructMatrixLapack(TPZCompMesh *);
    
  TPZSBStructMatrixLapack(TPZAutoPointer<TPZCompMesh> cmesh);
	
	TPZSBStructMatrixLapack(const TPZSBStructMatrixLapack &cp);
    
  ~TPZSBStructMatrixLapack();
	
  virtual TPZMatrix<STATE> * Create();
	
  virtual TPZStructMatrix * Clone();
    
public:
	
};

#endif //TPZSBSTRUCTMATRIXLAPACK_H
