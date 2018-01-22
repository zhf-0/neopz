
/**
 * @file
 * @brief Contains the implementation of the TPZAlgebraicSystemSolver methods.
 */

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzmatred"));
#endif

#include "TPZAlgebraicSystemSolver.h"

#include <stdlib.h>
using namespace std;

template <class TVar>
TPZAlgebraicSystemSolver<TVar>::TPZAlgebraicSystemSolver(TPZAutoPointer<TPZMatrix<TVar> > Refmat) :
fScratch()
{
    fContainer = Refmat;
}

template<class TVar>
TPZAlgebraicSystemSolver<TVar>::TPZAlgebraicSystemSolver() :
fScratch()
{
}

template <class TVar>
TPZAlgebraicSystemSolver<TVar>::TPZAlgebraicSystemSolver(const TPZAlgebraicSystemSolver<TVar> &Source) :
fScratch()
{
    fReferenceMatrix = Source.fReferenceMatrix;
    fContainer = Source.fContainer;
}

template <class TVar>
TPZAlgebraicSystemSolver<TVar>::~TPZAlgebraicSystemSolver()
{
}

template <class TVar>
void TPZAlgebraicSystemSolver<TVar>::ResetMatrix()
{
    TPZAutoPointer<TPZMatrix<TVar> > reset;
    fContainer = reset;
}

template <class TVar>
void TPZAlgebraicSystemSolver<TVar>::ShareMatrix(TPZAlgebraicSystemSolver<TVar> &other)
{
    if (this == &other)
        return;
    fContainer = other.fContainer;
}
template <class TVar>
void TPZAlgebraicSystemSolver<TVar>::Write(TPZStream &buf, int withclassid)
{
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Entering " << __PRETTY_FUNCTION__;
        LOGPZ_DEBUG(logger,sout.str());
    }
#endif
    TPZSolver<TVar>::Write(buf,withclassid);
    if(fContainer)
    {
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "fContainer AutoPointer valid on " << __PRETTY_FUNCTION__;
            LOGPZ_DEBUG(logger,sout.str());
        }
#endif
        
        fContainer->Write(buf, 1);
    }
    else
    {
        int flag = -1;
        buf.Write(&flag, 1);
    }
    if(fReferenceMatrix)
    {
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "fReferenceMatrix AutoPointer valid! It Shouldn't ! Expect Trouble " << __PRETTY_FUNCTION__;
            LOGPZ_WARN(logger,sout.str());
        }
#endif
        fReferenceMatrix->Write(buf, 1);
    }
    else
    {
        int flag = -1;
        buf.Write(&flag, 1);
    }
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Leaving" << __PRETTY_FUNCTION__;
        LOGPZ_DEBUG(logger,sout.str());
    }
#endif
}

template <class TVar>
void TPZAlgebraicSystemSolver<TVar>::Read(TPZStream &buf, void *context)
{
    TPZSolver<TVar>::Read(buf,context);
    fContainer = dynamic_cast<TPZMatrix<TVar> *>(TPZSaveable::Restore(buf, context));
    fReferenceMatrix = dynamic_cast<TPZMatrix<TVar> *>(TPZSaveable::Restore(buf, context));
}


template class TPZAlgebraicSystemSolver<float>;
template class TPZAlgebraicSystemSolver<std::complex<float> >;

template class TPZAlgebraicSystemSolver<double>;
template class TPZAlgebraicSystemSolver<std::complex<double> >;

template class TPZAlgebraicSystemSolver<long double>;
template class TPZAlgebraicSystemSolver<std::complex<long double> >;

