//
// Created by Francisco Teixeira Orlandini on 11/27/17.
//

#ifndef PZ_TPZLAPACKWRAPPER_H
#define PZ_TPZLAPACKWRAPPER_H

#include <pzfmatrix.h>
#include <pzsbndmat.h>

template<class TVar>
class SPZAlwaysComplex;

template<class TVar>
class TPZLapackWrapper {
    friend class TPZFMatrix<TVar>;
    friend class TPZSBMatrix<TVar>;
    /*******************
    *    TPZFMATRIX    *
    *******************/
    static int SolveEigenProblem(TPZFMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors);

    /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors
     * @param w Stores the eigenvalues
     */
    static int SolveEigenProblem(TPZFMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w);

    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
     * @param w Stores the eigenvalues
     * @param Stores the correspondent eigenvectors
     */
    static int SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &A, TPZFMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors);
    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
     * @param w Stores the eigenvalues
     */
    static int SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &A, TPZFMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w);

    /*******************
    *    TPSBMATRIX    *
    *******************/
    static int SolveEigenProblem(TPZSBMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors);

    /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors
     * @param w Stores the eigenvalues
     */
    static int SolveEigenProblem(TPZSBMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w);

    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors
     * @param w Stores the eigenvalues
     * @param Stores the correspondent eigenvectors
     */
    static int SolveGeneralisedEigenProblem(TPZSBMatrix<TVar> &A, TPZSBMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors);
    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors
     * @param w Stores the eigenvalues
     */
    static int SolveGeneralisedEigenProblem(TPZSBMatrix<TVar> &A, TPZSBMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w);
};

#endif //PZ_TPZLAPACKWRAPPER_H
