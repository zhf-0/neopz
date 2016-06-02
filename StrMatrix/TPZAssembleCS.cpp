//
//  TPZAssembleCS.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#include "TPZAssembleCS.h"

#include "pz_pthread.h"
#include "run_stats_table.h"
#include "pzlog.h"
#include "pzcompel.h"
#include "pzsubcmesh.h"
#include "pzmaterial.h"

#include "TPZTimer.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZAssembleCS"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


void TPZAssembleCS::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs)
{
    ThreadData threaddata(this,mat,rhs,fConfig->fMaterialIds);

    const int numthreads = this->fConfig->fNumThreads;
    TPZVec<pthread_t> allthreads(numthreads);
    int itr;
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_CREATE(&allthreads[itr], NULL,ThreadData::ThreadWork,
                          &threaddata, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
    }
    
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


void TPZAssembleCS::MultiThread_Assemble(TPZFMatrix<STATE> & rhs)
{
    ThreadData threaddata(this,rhs,fConfig->fMaterialIds);
    
    const int numthreads = this->fConfig->fNumThreads;
    TPZVec<pthread_t> allthreads(numthreads);
    int itr;
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_CREATE(&allthreads[itr], NULL,ThreadData::ThreadWork,
                          &threaddata, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
    }
}



TPZAssembleCS::ThreadData::ThreadData(TPZAssembleCS *assemble, TPZMatrix<STATE> &mat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds)
: fAssemble(assemble), fGlobMatrix(&mat), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZAssembleCS::ThreadData::ThreadData()");
    PZ_PTHREAD_MUTEX_INIT(&fAccessElementK,NULL,"");
    PZ_PTHREAD_MUTEX_INIT(&fAccessElementF,NULL,"");
    /*	sem_t *sem_open( ... );
     int sem_close(sem_t *sem);
     int sem_unlink(const char *name);
     */
    /*
     #ifdef MACOSX
     std::stringstream sout;
     static int counter = 0;
     sout << "AssemblySem" << counter++;
     fAssembly = sem_open(sout.str().c_str(), O_CREAT,777,1);
     if(fAssembly == SEM_FAILED)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     DebugStop();
     }
     #else
     int sem_result = sem_init(&fAssembly,0,0);
     if(sem_result != 0)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     }
     #endif
     */
}

TPZAssembleCS::ThreadData::ThreadData(TPZAssembleCS *assemble,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds)
: fAssemble(assemble), fGlobMatrix(0), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZAssembleCS::ThreadData::ThreadData()");
    PZ_PTHREAD_MUTEX_INIT(&fAccessElementK,NULL,"");
    PZ_PTHREAD_MUTEX_INIT(&fAccessElementF,NULL,"");
    /*	sem_t *sem_open( ... );
     int sem_close(sem_t *sem);
     int sem_unlink(const char *name);
     */
    /*
     #ifdef MACOSX
     std::stringstream sout;
     static int counter = 0;
     sout << "AssemblySem" << counter++;
     fAssembly = sem_open(sout.str().c_str(), O_CREAT,777,1);
     if(fAssembly == SEM_FAILED)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     DebugStop();
     }
     #else
     int sem_result = sem_init(&fAssembly,0,0);
     if(sem_result != 0)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     }
     #endif
     */
}

TPZAssembleCS::ThreadData::~ThreadData()
{
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"TPZAssembleCS::ThreadData::~ThreadData()");
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElementK,"");
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElementF,"");
    
    /*
     
     #ifdef MACOSX
     sem_close(fAssembly);
     #else
     sem_destroy(&fAssembly);
     #endif
     */
}

void *TPZAssembleCS::ThreadData::ThreadWork(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    // compute the next element (this method is threadsafe)
    long iel = data->NextElement();
    TPZCompMesh *cmesh = data->fAssemble->fConfig->fMesh;
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    long nel = cmesh->NElements();
    while(iel < nel)
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
            
            if(data->fAssemble->fConfig->fEquationFilter.IsActive())
            {
                data->fAssemble->fConfig->fEquationFilter.Filter(ek->fSourceIndex,ek->fDestinationIndex);
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
            if(data->fAssemble->fConfig->fEquationFilter.IsActive())
            {
                data->fAssemble->fConfig->fEquationFilter.Filter(ek->fSourceIndex,ek->fDestinationIndex);
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
#ifdef LOG4CXX
        timeforel.stop();
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << timeforel.processName() <<  timeforel;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        
        // compute the next element (this method is threadsafe)
        iel = data->NextElement();
    }
    
    return 0;
}

long TPZAssembleCS::ThreadData::NextElement()
{
    TPZCompMesh *cmesh = fAssemble->fConfig->fMesh;
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    
    long nel = elementvec.NElements();
    long my_el;
    
    while (1) {
        
        PZ_PTHREAD_MUTEX_LOCK(&fAccessElement,"TPZStructMatrix::ThreadData::NextElement()");
        my_el = fNextElement++;
        PZ_PTHREAD_MUTEX_UNLOCK(&fAccessElement,"TPZStructMatrix::ThreadData::NextElement()");
        
        if (my_el >= nel-1)
            break;
        
        TPZCompEl *el = elementvec[my_el];
        if(!el) continue;
        if(fAssemble->fConfig->fMaterialIds.size() == 0) break;
        TPZMaterial * mat = el->Material();
        TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
        if(!mat)
        {
            if(!submesh)
            {
                continue;
            }
            else if(submesh->NeedsComputing(fAssemble->fConfig->fMaterialIds) == false) continue;
        }
        else
        {
            int matid = mat->Id();
            if(this->ShouldCompute(matid) == false) continue;
        }
        break;
    }
    
    return my_el;
}

