//
// Created by Francisco Teixeira Orlandini on 11/23/17.
//

#ifndef PZ_TPZEIGENSOLVER_H
#define PZ_TPZEIGENSOLVER_H

#include "TPZSolver.h"

/**
* @ingroup solver
* @brief  Defines a class of solvers for eigenvalue problems. \ref solver "Solver"
*/
template <typename TVar>
class TPZEigenSolver : public TPZSolver<TVar> {
public:
    /**
     * @brief Method for solving generalized(or not) eigenvalue problems.
     * @param B Rhs matrix. For generalized problems.
     * @param eigenVectors If calculated, eigenvectors are stored in a column-wise fashion
     * @param eigenValues eigenvalues are stored in the first (and only) column.
     */
//    void Solve(const TPZFMatrix<TVar> &B = 0, TPZFMatrix<TVar> &eigenVectors,
//                       TPZFMatrix<TVar>  *eigenValues) override;

};


#endif //PZ_TPZEIGENSOLVER_H
