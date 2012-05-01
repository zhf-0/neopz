/**
 * @file
 * @brief Contains TPZStepSolver class which defines step solvers class.
 */

#ifndef TPZSTEPSOLVER_H
#define TPZSTEPSOLVER_H
#include "pzsolve.h"

#include "pzstream.h"

#include <list>

template<class TVar>
class TPZFMatrix;


/**
 * @brief Defines step solvers class. \ref solver "Solver"
 * @ingroup solver
 */
template<class TVar>
class TPZStepSolver: public TPZMatrixSolver<TVar>
{
public:
	TPZStepSolver(TPZAutoPointer<TPZMatrix<TVar> > refmat = 0);
	
	TPZStepSolver(const TPZStepSolver<TVar> & copy);
	
	virtual ~TPZStepSolver();
	
	void SetSOR(const int numiterations, const REAL overrelax, const REAL tol,
				const int FromCurrent);
	
	void SetSSOR(const int numiterations, const REAL overrelax, const REAL tol,
				 const int FromCurrent);
	
	void
	SetJacobi(const int numiterations, const REAL tol, const int FromCurrent);
	
	void SetCG(const int numiterations, const TPZMatrixSolver<TVar> &pre,
			   const REAL tol, const int FromCurrent);
	
	void SetGMRES(const int numiterations, const int numvectors,
				  const TPZMatrixSolver<TVar> &pre, const REAL tol, const int FromCurrent);
	
	void SetBiCGStab(const int numiterations, const TPZMatrixSolver<TVar> &pre,
					 const REAL tol, const int FromCurrent);
	
	void SetDirect(const DecomposeType decomp);
	
	void SetMultiply();
	
	virtual TPZSolver<TVar> *Clone() const
	{
		return new TPZStepSolver<TVar>(*this);
	}
	
	void SetTolerance(REAL tol)
	{
		fTol = tol;
	}
	
    /** @brief reset the data structure of the solver object */
	void ResetSolver();
    
    virtual typename TPZMatrixSolver<TVar>::MSolver Solver()
    {
        return fSolver;
    }
	
	/** @brief returns the equations for which the equations had zero pivot */
	std::list<int> &Singular()
	{
		return fSingular;
	}
	
	/** @brief This method will reset the matrix associated with the solver */
	/** This is useful when the matrix needs to be recomputed in a non linear problem */
	virtual void ResetMatrix();

	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
	{
		if (fPrecond)
			fPrecond->UpdateFrom(matrix);
		TPZMatrixSolver<TVar>::UpdateFrom(matrix);
	}
	
    
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual = 0);
    
    /** @brief Decompose the system of equations if a direct solver is used */
    virtual void Decompose();
    
    /** @brief Define the preconditioner as a solver object */
	void SetPreconditioner(TPZSolver<TVar> &solve);
    
    /** @brief access method to the preconditioner */
    TPZSolver<TVar> *PreConditioner()
    {
        return fPrecond;
    }
	
	/** @brief Serialization methods */
	virtual int ClassId() const;
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
	
	
private:
	typename TPZMatrixSolver<TVar>::MSolver fSolver;
	DecomposeType fDecompose;
	int fNumIterations;
	int fNumVectors;
	REAL fTol;
	REAL fOverRelax;
	
	/** @brief Solver using preconditioner matrix */
	TPZSolver<TVar> *fPrecond;
	int fFromCurrent;
	
	std::list<int> fSingular;
};

#endif
