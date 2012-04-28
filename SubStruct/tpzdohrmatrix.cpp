/**
 * @file
 * @brief Contains the implementation of the TPZDohrMatrix methods. 
 */
/***************************************************************************
 *   Copyright (C) 2006 by Philippe Devloo   *
 *   phil@fec.unicamp.br   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "tpzdohrmatrix.h"
//#include "tpzdohrsubstruct.h"

#include "tpzdohrassembly.h"
#include "pzlog.h"

#include "TPZfTime.h"
#include "TPZTimeTemp.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("substruct.dohrsubstruct"));
#endif


template<class TVar, class TSubStruct>
TPZDohrMatrix<TVar,TSubStruct>::TPZDohrMatrix(TPZAutoPointer<TPZDohrAssembly<TVar> > assembly)
: TPZMatrix<TVar>(), fNumThreads(0), fAssembly(assembly)
{
}


template<class TVar, class TSubStruct>
TPZDohrMatrix<TVar,TSubStruct>::~TPZDohrMatrix()
{
}


template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
										const TVar alpha,const TVar beta,const int opt,const int stride) const
{
	TPZfTime mult;
	if ((!opt && this->Cols() != x.Rows()*stride) || this->Rows() != x.Rows()*stride)
		this->Error( "Operator* <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		this->Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}
	this->PrepareZ(y,z,beta,opt,stride);
	
	typename SubsList::const_iterator iter;
	int isub = 0;
	if (fNumThreads == 0) {
		for (iter=fGlobal.begin();iter!=fGlobal.end();iter++,isub++) {
			TPZFMatrix<TVar> xlocal,zlocal;
			fAssembly->Extract(isub,x,xlocal);
			zlocal.Redim(xlocal.Rows(),xlocal.Cols());
			(*iter)->ContributeKULocal(alpha,xlocal,zlocal);
			fAssembly->Assemble(isub,zlocal,z);
			//         z.Print("Resultado intermediario");
		}		
	}
	else {
		TPZAutoPointer<TPZDohrAssembleList<TVar> > assemblelist = new TPZDohrAssembleList<TVar>(fGlobal.size(),z,this->fAssembly);
		
		TPZDohrThreadMultList<TVar,TSubStruct> multwork(x,alpha,fAssembly,assemblelist);
		typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator iter;
		int isub=0;
		for (iter=fGlobal.begin(); iter!=fGlobal.end(); iter++,isub++) {
			TPZDohrThreadMultData<TSubStruct> data(isub,*iter);
			multwork.AddItem(data);
		}
		TPZVec<pthread_t> AllThreads(fNumThreads+1);
		int i;
		for (i=0; i<fNumThreads; i++) {
			pthread_create(&AllThreads[i+1], 0, TPZDohrThreadMultList<TVar,TSubStruct>::ThreadWork, &multwork);
		}
		pthread_create(&AllThreads[0], 0, TPZDohrAssembleList<TVar>::Assemble, assemblelist.operator->());
		
		for (i=0; i<fNumThreads+1; i++) {
			void *result;
			pthread_join(AllThreads[i], &result);
		}
	}
	tempo.fMultiply.Push(mult.ReturnTimeDouble());
}

template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::Initialize() 
{
	std::cout << "Number of substructures " << fGlobal.size() << std::endl;
	tempo.fNumSub = fGlobal.size();																// alimenta timeTemp com o numero de substruturas
	TPZFMatrix<TVar> diag(this->Rows(),1,0.);
	typename SubsList::iterator iter;
	int isub = 0;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++,isub++) {
        //Basic initialization for each substructure (compute the matrices)
        //(*iter)->Initialize();
		TPZFMatrix<TVar> diaglocal;
        (*iter)->ContributeDiagonalLocal(diaglocal);
		LOGPZ_DEBUG(logger,"Before assemble diagonal")
		this->fAssembly->Assemble(isub,diaglocal,diag);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Substructure " << isub << " ";
			diag.Print("Global Diagonal matrix",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
        std::cout << '*';
        std::cout.flush();
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		diag.Print("Global Diagonal matrix",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	std::cout << std::endl;
	for (iter=fGlobal.begin(),isub=0;iter!=fGlobal.end();iter++,isub++) {
        //Computes the Weights for each substructure
		TPZFMatrix<TVar> diaglocal;
		this->fAssembly->Extract(isub,diag,diaglocal);
        (*iter)->ComputeWeightsLocal(diaglocal);
		
	}
}

/**
 * Adjust the residual to zero the residual of the internal connects
 */
template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::AdjustResidual(TPZFMatrix<TVar> &res)
{
	typename SubsList::iterator iter;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++) {
		(*iter)->AdjustResidual(res);
	}
}

/**
 * Add the solution corresponding to the internal residual
 */
template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::AddInternalSolution(TPZFMatrix<TVar> &solution)
{
	typename SubsList::iterator iter;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++) {
		(*iter)->AddInternalSolution(solution);
	}
}

template<class TVar, class TSubStruct>
void *TPZDohrThreadMultList<TVar,TSubStruct>::ThreadWork(void *ptr)
{
	TPZDohrThreadMultList<TVar,TSubStruct> *myptr = (TPZDohrThreadMultList<TVar,TSubStruct> *) ptr;
	TPZDohrThreadMultData<TSubStruct> runner = myptr->PopItem();
	while (runner.IsValid()) {
		TPZFMatrix<TVar> xlocal;
		myptr->fAssembly->Extract(runner.fisub,*(myptr->fInput),xlocal);
		TPZAutoPointer<TPZDohrAssembleItem<TVar> > assembleItem = new TPZDohrAssembleItem<TVar>(runner.fisub,xlocal.Rows(),xlocal.Cols());
		runner.fSub->ContributeKULocal(myptr->fAlpha,xlocal,assembleItem->fAssembleData);
		myptr->fAssemblyStructure->AddItem(assembleItem);
		runner = myptr->PopItem();
	}
	return ptr;
}

template <>
int TPZDohrMatrix<double, TPZDohrSubstructCondense<double> >::ClassId() const {
    return TPZDOHRMATRIXSUBSTRUCTCONDENSE;
}

template <>
int TPZDohrMatrix<double, TPZDohrSubstruct<double> >::ClassId() const {
    return TPZDOHRMATRIXSUBSTRUCT;
}

/**
 * @brief Unpacks the object structure from a stream of bytes
 * @param buf The buffer containing the object in a packed form
 * @param context 
 */
template <>
void TPZDohrMatrix<double, TPZDohrSubstructCondense<double> >::Read(TPZStream &buf, void *context )
{
    TPZMatrix<double>::Read(buf, context);
    fAssembly = new TPZDohrAssembly<double>;
    fAssembly->Read(buf);
    buf.Read(&fNumCoarse);
    buf.Read(&fNumThreads);
    int sz;
    buf.Read(&sz);
    for (int i=0; i<sz; i++) {
        TPZAutoPointer<TPZDohrSubstructCondense<double> > sub = new TPZDohrSubstructCondense<double>;
        sub->Read(buf);
        fGlobal.push_back(sub);
    }
}
/**
 * @brief Packs the object structure in a stream of bytes
 * @param buf Buffer which will receive the bytes
 * @param withclassid
 */
template <>
void TPZDohrMatrix<double,TPZDohrSubstructCondense<double> >::Write( TPZStream &buf, int withclassid )
{
    TPZMatrix<double>::Write(buf, withclassid);
    fAssembly->Write(buf);
    buf.Write(&fNumCoarse);
    buf.Write(&fNumThreads);
    SubsList::iterator it;
    int size = fGlobal.size();
    buf.Write(&size);
    for (it=fGlobal.begin(); it != fGlobal.end(); it++) {
        (*it)->Write(buf);
    }
    size = 0;
}

template <>
void TPZDohrMatrix<double, TPZDohrSubstruct<double> >::Read(TPZStream &buf, void *context )
{
    DebugStop();
}

template <>
void TPZDohrMatrix<double, TPZDohrSubstruct<double> >::Write( TPZStream &buf, int withclassid )
{
    DebugStop();
}

template class TPZDohrMatrix<double, TPZDohrSubstruct<double> >;
template class TPZDohrMatrix<double, TPZDohrSubstructCondense<double> >;

template class TPZRestoreClass< TPZDohrMatrix<double, TPZDohrSubstructCondense<double> > , TPZDOHRMATRIXSUBSTRUCTCONDENSE>;
template class TPZRestoreClass< TPZDohrMatrix<double, TPZDohrSubstruct<double> > , TPZDOHRMATRIXSUBSTRUCT>;


