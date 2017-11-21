

#ifndef TPZALGEBRAICSYSTEMSOLVERH
#define TPZALGEBRAICSYSTEMSOLVERH

#include "TPZSolver.h"
#include "pzfmatrix.h"

//@TODO: Document the class
/**
 * @ingroup solver
 * @brief  Defines a class of matrix solvers. \ref solver "Solver"
 */
template<class TVar>
class TPZAlgebraicSystemSolver : public TPZSolver<TVar>
{
    
public:
    /**
     * @enum MSolver
     * @brief Defines a series of solvers available in PZ
     * @param ENoSolver No solver selected
     * @param EJacobi Jacobi solver selected
     * @param ESOR Successive Over Relaxation solver selected
     * @param ESSOR Symmetric Successive Over Relaxation solver selected
     * @param ECG Conjugate Gradiente solver selected
     * @param EDirect LU, LDLt or Cholesky selected
     */
    enum MSolver
    {
        ENoSolver, EJacobi, ESOR, ESSOR, ECG, EGMRES, EBICGSTAB, EDirect, EMultiply
    };
    
    /**
     * @brief Constructor with initialization parameter
     * @param Refmat Sets reference matrix to 0
     */
    
    TPZAlgebraicSystemSolver(TPZAutoPointer<TPZMatrix<TVar> >  Refmat);
    
    TPZAlgebraicSystemSolver();
    
    /**
     * @brief Copy constructor
     * @param Source Model object to be copied from
     */
    TPZAlgebraicSystemSolver(const TPZAlgebraicSystemSolver<TVar> &Source);
    
    /** @brief Destructor */
    virtual ~TPZAlgebraicSystemSolver();
    
    /**
     * @brief Sets a matrix to the current object
     * @param Refmat Sets reference matrix to RefMat
     */
    virtual void SetMatrix(TPZAutoPointer<TPZMatrix<TVar> > Refmat)
    {
        fContainer = Refmat;
    }
    
    /** @brief Updates the values of the current matrix based on the values of the matrix */
    virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
    {
        if (fReferenceMatrix == matrix && matrix)
        {
            if(this->fContainer) this->fContainer->UpdateFrom(matrix);
        }
    }
    /** @brief Resets current object */
    void ResetMatrix();
    
    /** @brief This method gives a preconditioner to share a matrix with the referring solver object */
    virtual void SetReferenceMatrix(TPZAutoPointer<TPZMatrix<TVar> > matrix)
    {
        fReferenceMatrix = matrix;
    }
    
    /** @brief Returns a pointer to TPZMatrix<>*/
    TPZAutoPointer<TPZMatrix<TVar> > Matrix() const
    {
        return fContainer;
    }
    
    void ReallocMatrix() {
        fContainer.ReallocForNuma(0);
    }
    
    /**
     * @brief Shares the current matrix with another object of same type
     * @param other Object that will share current matrix
     */
    void ShareMatrix(TPZAlgebraicSystemSolver<TVar> & other);
    
    virtual MSolver Solver()
    {
        return ENoSolver;
    }
    
protected:
    
private:
    /** @brief Container classes */
    TPZAutoPointer<TPZMatrix<TVar> > fContainer;
protected:
    /** @brief Reference matrix used to update the current matrix */
    TPZAutoPointer<TPZMatrix<TVar> > fReferenceMatrix;
    
protected:
    /** @brief Manipulation matrix */
    TPZFMatrix<TVar>  fScratch;
public:
    /** @brief Saveable specific methods */
    virtual int ClassId() const
    {
        return 666; //@TODO: Reimplement using hash function
    }
    virtual void Write(TPZStream &buf, int withclassid);
    virtual void Read(TPZStream &buf, void *context);
};

#endif  // TPZALGEBRAICSYSTEMSOLVERH
