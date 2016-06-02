//
//  TPZAssembleCS.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#include "TPZAssembleCSTBB.h"

#include "pz_pthread.h"
#include "run_stats_table.h"
#include "pzlog.h"
#include "pzcompel.h"
#include "pzsubcmesh.h"
#include "pzmaterial.h"

#include "TPZTimer.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZAssembleCSTBB"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


void TPZAssembleCSTBB::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs)
{
#ifdef USING_TBB
    ThreadData threaddata(this,mat,rhs,fConfig->fMaterialIds);
    AssembleTask tasks(&threaddata);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, fConfig->fMesh->NElements()), tasks);
#else
    std::cout << __PRETTY_FUNCTION__ << " please link with TBB\n";
    TPZAssembleCS::MultiThread_Assemble(mat, rhs);

#endif
    
#ifdef LOG4CXX
    if(loggerCheck->isDebugEnabled())
    {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        mat.Print("Matriz de Rigidez: ",sout,EMathematicaInput);
        rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
#endif
}


void TPZAssembleCSTBB::MultiThread_Assemble(TPZFMatrix<STATE> & rhs)
{
    
#ifdef USING_TBB
    ThreadData threaddata(this,rhs,fConfig->fMaterialIds);
    AssembleTask tasks(&threaddata);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, fConfig->fMesh->NElements()), tasks);
#else
    std::cout << __PRETTY_FUNCTION__ << " please link with TBB\n";
    TPZAssembleCS::MultiThread_Assemble(rhs);

#endif
}




#ifdef USING_TBB
void TPZAssembleCSTBB::AssembleTask::operator()(const tbb::blocked_range<size_t>& range) const
{
    TPZCompMesh *cmesh = data->fAssemble->Config()->fMesh;
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    
    for(size_t iel=range.begin(); iel!=range.end(); ++iel )
    {
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "Computing element " << iel;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#ifdef LOG4CXX
        std::stringstream sout;
        sout << "Element " << iel << " elapsed time ";
        TPZTimer timeforel(sout.str());
        timeforel.start();
#endif
        
        TPZAutoPointer<TPZElementMatrix> ek;
        TPZAutoPointer<TPZElementMatrix> ef = new TPZElementMatrix(cmesh,TPZElementMatrix::EF);
        if (data->fGlobMatrix) {
            ek = new TPZElementMatrix(cmesh,TPZElementMatrix::EK);
        }
        else
        {
            ek = ef;
        }
        
        TPZCompEl *el = cmesh->ElementVec()[iel];
        
        if (el) {
            
            TPZElementMatrix *ekp = ek.operator->();
            TPZElementMatrix *efp = ef.operator->();
            TPZElementMatrix &ekr = *ekp;
            TPZElementMatrix &efr = *efp;
            
            if (data->fGlobMatrix) {
                el->CalcStiff(ekr,efr);
            }
            else
            {
                el->CalcResidual(efr);
            }
            
            if(guiInterface) if(guiInterface->AmIKilled()){
                break;
            }
            
            if(!el->HasDependency()) {
                ek->ComputeDestinationIndices();
                
                if(data->fAssemble->Config()->fEquationFilter.IsActive())
                {
                    data->fAssemble->Config()->fEquationFilter.Filter(ek->fSourceIndex,ek->fDestinationIndex);
                }
#ifdef LOG4CXX
                if(loggerel->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element index " << iel << std::endl;
                    ek->fMat.Print("Element stiffness matrix",sout);
                    ef->fMat.Print("Element right hand side", sout);
                    LOGPZ_DEBUG(loggerel,sout.str())
                }
#endif
            } else {
                // the element has dependent nodes
                if (data->fGlobMatrix) {
                    ek->ApplyConstraints();
                }
                ef->ApplyConstraints();
                ek->ComputeDestinationIndices();
                if(data->fAssemble->Config()->fEquationFilter.IsActive())
                {
                    data->fAssemble->Config()->fEquationFilter.Filter(ek->fSourceIndex,ek->fDestinationIndex);
                }
#ifdef LOG4CXX
                if(loggerel2->isDebugEnabled() && el->Reference() &&  el->Reference()->MaterialId() == 1 && el->IsInterface())
                {
                    std::stringstream sout;
                    el->Reference()->Print(sout);
                    el->Print(sout);
                    ek->Print(sout);
                    //			ef->Print(sout);
                    LOGPZ_DEBUG(loggerel2,sout.str())
                }
#endif
#ifdef LOG4CXX
                if(loggerel->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element index " << iel << std::endl;
                    ek->fConstrMat.Print("Element stiffness matrix",sout);
                    ef->fConstrMat.Print("Element right hand side", sout);
                    LOGPZ_DEBUG(loggerel,sout.str())
                }
#endif
            }
            
#ifdef LOG4CXX
            timeforel.stop();
            if (logger->isDebugEnabled())
            {
                std::stringstream sout;
                sout << timeforel.processName() <<  timeforel;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            // Assemble the matrix
            
            if(!ek->HasDependency())
            {
                if (data->fGlobMatrix) {
                    PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElementK,"");
                    data->fGlobMatrix->AddKel(ek->fMat,ek->fSourceIndex,ek->fDestinationIndex);
                    PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElementK,"");
                    
                }
                PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElementF,"");
                data->fGlobRhs->AddFel(ef->fMat,ek->fSourceIndex,ek->fDestinationIndex);
                PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElementF,"");
            }
            else
            {
                if (data->fGlobMatrix) {
                    PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElementK,"");
                    data->fGlobMatrix->AddKel(ek->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
                    PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElementK,"");
                }
                PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElementF,"");
                data->fGlobRhs->AddFel(ef->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
                PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElementF,"");
            }
        }
    }
}
#endif
