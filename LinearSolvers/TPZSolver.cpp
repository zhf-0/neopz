/**
 * @file
 * @brief Contains the implementation of the TPZSolver methods.
 */

#include "TPZSolver.h"

/** Destructor */
template <class TVar>
TPZSolver<TVar>::~TPZSolver()
{
}

template class TPZSolver<float>;
template class TPZSolver<std::complex<float> >;

template class TPZSolver<double>;
template class TPZSolver<std::complex<double> >;

template class TPZSolver<long double>;
template class TPZSolver<std::complex<long double> >;
