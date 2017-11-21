/**
 * @file
 * @brief Contains TPZSolver class which defines a abstract class of solvers  which will be used by matrix classes.
 */

#ifndef TPZSOLVERH
#define TPZSOLVERH


#include "pzsave.h"
#include "pzfmatrix.h"


/**
 * @ingroup solver
 * @brief Defines a abstract class of solvers  which will be used by matrix classes. \ref solver "Solver"
 */
template<class TVar>
class TPZSolver: public TPZSaveable
{

public:
    
    /**
     * @brief Default constructor
     */
    TPZSolver() = default;
    
	/**
	 * @brief Solves the system of linear equations
	 * @param F contains Force vector
	 * @param result contains the solution
	 * @param residual contains the residual for that linear system
	 */
	virtual void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
					   TPZFMatrix<TVar>  *residual = 0) = 0;
    
    /** @brief Decompose the system of equations if a direct solver is used */
    virtual void Decompose()
    {
    }
	
	/** @brief Clones the current object returning a pointer of type TPZSolver */
	virtual TPZSolver *Clone() const = 0;
	/** @brief Destructor */
	virtual ~TPZSolver();
	
	/** @brief This method will reset the matrix associated with the solver */
	/** This is useful when the matrix needs to be recomputed in a non linear problem */
	virtual void ResetMatrix()
	{
	}
	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
	{
		std::cout << __PRETTY_FUNCTION__ << " called\n";
	}

};

#endif  // TPZSOLVERH
