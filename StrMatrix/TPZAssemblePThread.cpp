//
//  TPZAssemblePThread.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#include "TPZAssemblePThread.h"

#include "pzlog.h"
#include "pzsubcmesh.h"
#include "pzmaterial.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZAssemble"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
static LoggerPtr loggerGlobStiff(Logger::getLogger("pz.strmatrix.globalstiffness"));
#endif



void TPZAssemblePThread::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs)
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
    
    ThreadData::ThreadAssembly(&threaddata);
    
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


void TPZAssemblePThread::MultiThread_Assemble(TPZFMatrix<STATE> & rhs)
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
    
    ThreadData::ThreadAssembly(&threaddata);
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
    }
}



TPZAssemblePThread::ThreadData::ThreadData(TPZAssemblePThread *assemble, TPZMatrix<STATE> &mat,
                                        TPZFMatrix<STATE> &rhs,
                                        std::set<int> &MaterialIds)
: fAssemble(assemble), fGlobMatrix(&mat), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZStructMatrix::ThreadData::ThreadData()");
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

TPZAssemblePThread::ThreadData::ThreadData(TPZAssemblePThread *assemble,
                                        TPZFMatrix<STATE> &rhs,
                                        std::set<int> &MaterialIds)
: fAssemble(assemble), fGlobMatrix(0), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZStructMatrix::ThreadData::ThreadData()");
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

TPZAssemblePThread::ThreadData::~ThreadData()
{
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"TPZStructMatrix::ThreadData::~ThreadData()");
    /*
     #ifdef MACOSX
     sem_close(fAssembly);
     #else
     sem_destroy(&fAssembly);
     #endif
     */
}

//#define DRY_RUN

void *TPZAssemblePThread::ThreadData::ThreadWork(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    // compute the next element (this method is threadsafe)
    long iel = data->NextElement();
    TPZCompMesh *cmesh = data->fAssemble->fConfig->fMesh;
    long nel = cmesh->NElements();
    while(iel < nel)
    {
        
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
        
#ifndef DRY_RUN
        if (data->fGlobMatrix) {
            el->CalcStiff(ekr,efr);
        }
        else
        {
            el->CalcResidual(efr);
        }
#else
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(el);
            if (!intel) {
                DebugStop();
            }
            if (data->fGlobMatrix) {
                intel->InitializeElementMatrix(ekr, efr);
            }
            else
            {
                intel->InitializeElementMatrix(efr);
            }
        }
#endif
        
        if(!el->HasDependency()) {
            ek->ComputeDestinationIndices();
            
            if(data->fAssemble->fConfig->fEquationFilter.IsActive())
            {
                data->fAssemble->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
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
                data->fAssemble->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
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
        
        
        // put the elementmatrices on the stack to be assembled (threadsafe)
        data->ComputedElementMatrix(iel,ek,ef);
        // compute the next element (this method is threadsafe)
        iel = data->NextElement();
    }
    PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElement,"TPZStructMatrix::ThreadData::ThreadWork");
    data->fAssembly.Post();
    /*
     #ifdef MACOSX
     sem_post(data->fAssembly);
     #else
     sem_post(&data->fAssembly);
     #endif
     */
    PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElement,"TPZStructMatrix::ThreadData::ThreadWork");
    
    return 0;
}

// The function which will compute the assembly
void *TPZAssemblePThread::ThreadData::ThreadAssembly(void *threaddata)
{
    ThreadData *data = (ThreadData *) threaddata;
    TPZCompMesh *cmesh = data->fAssemble->fConfig->fMesh;
    long nel = cmesh->NElements();
    PZ_PTHREAD_MUTEX_LOCK(&(data->fAccessElement),"TPZStructMatrix::ThreadData::ThreadAssembly");
    long nextel = data->fNextElement;
    int numprocessed = data->fProcessed.size();
    while(nextel < nel || numprocessed)
    {
        std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > >::iterator itavail;
        std::set<int>::iterator itprocess;
        bool keeplooking = false;
        if(data->fSubmitted.size() && data->fProcessed.size())
        {
            itavail = data->fSubmitted.begin();
            itprocess = data->fProcessed.begin();
            if(itavail->first == *itprocess)
            {
                // make sure we come back to look for one more element
                keeplooking = true;
                // Get a hold of the data
#ifdef LOG4CXX
                int iel = *itprocess;
#endif
                data->fProcessed.erase(itprocess);
                TPZAutoPointer<TPZElementMatrix> ek = itavail->second.first;
                TPZAutoPointer<TPZElementMatrix> ef = itavail->second.second;
                data->fSubmitted.erase(itavail);
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Assembling element " << iel;
                    LOGPZ_DEBUG(logger,sout.str())
                }
#endif
#ifdef CHECKCONSISTENCY
                //static TPZCheckConsistency stiffconsist("ElementStiff");
                stiffconsist.SetOverWrite(true);
                bool result;
                result = stiffconsist.CheckObject(ek->fMat);
                if(!result)
                {
                    globalresult = false;
                    std::stringstream sout;
                    sout << "element " << iel << " computed differently";
                    LOGPZ_ERROR(loggerCheck,sout.str())
                }
#endif
                
                // Release the mutex
                PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElement,"TPZStructMatrix::ThreadData::ThreadAssembly");
                
#ifndef DRY_RUN
                // Assemble the matrix
                if(!ek->HasDependency())
                {
                    if (data->fGlobMatrix) {
                        data->fGlobMatrix->AddKel(ek->fMat,ek->fSourceIndex,ek->fDestinationIndex);
                    }
                    data->fGlobRhs->AddFel(ef->fMat,ek->fSourceIndex,ek->fDestinationIndex);
                }
                else
                {
                    if (data->fGlobMatrix) {
                        data->fGlobMatrix->AddKel(ek->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
                    }
                    data->fGlobRhs->AddFel(ef->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
                }
#endif
                // acquire the mutex
                PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElement,"TPZStructMatrix::ThreadData::ThreadAssembly");
            }
        }
        if(!keeplooking)
        {
            PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElement,"TPZStructMatrix::ThreadData::ThreadAssembly");
#ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                LOGPZ_DEBUG(logger,"Going to sleep within assembly")
            }
#endif
            // wait for a signal
            data->fAssembly.Wait();
            /*
             #ifdef MACOSX
             sem_wait(data->fAssembly);
             #else
             sem_wait(&data->fAssembly);
             #endif
             */
#ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                LOGPZ_DEBUG(logger,"Waking up for assembly")
            }
#endif
            PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElement,"TPZStructMatrix::ThreadData::ThreadAssembly");
        }
        nextel = data->fNextElement;
        numprocessed = data->fProcessed.size();
        
    }
    //	std::cout << std::endl;
#ifdef LOG4CXX
    if(loggerCheck->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "nextel = " << nextel << " numprocessed = " << numprocessed << " submitted " << data->fSubmitted.size() << std::endl;
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
#endif
    PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElement,"TPZStructMatrix::ThreadData::ThreadAssembly");
    return 0;
}

long TPZAssemblePThread::ThreadData::NextElement()
{
    PZ_PTHREAD_MUTEX_LOCK(&fAccessElement,"TPZStructMatrix::ThreadData::NextElement()");
    long iel;
    long nextel = fNextElement;
    TPZCompMesh *cmesh = fAssemble->fConfig->fMesh;
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    long nel = elementvec.NElements();
    for(iel=fNextElement; iel < nel; iel++)
    {
        TPZCompEl *el = elementvec[iel];
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
    fNextElement = iel+1;
    nextel = iel;
    if(iel<nel) fProcessed.insert(iel); //AQUIBORIN pelo que percebi, aqui tem que acontecer antes do unlock no caso sem Critical Section
    PZ_PTHREAD_MUTEX_UNLOCK(&fAccessElement,"TPZStructMatrix::ThreadData::NextElement()");
#ifdef LOG4CXX
    {
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " returning " << nextel << " fNextElement " << fNextElement;
            LOGPZ_DEBUG(logger,sout.str())
        }
    }
#endif
    return nextel;
}


// put the computed element matrices in the map
void TPZAssemblePThread::ThreadData::ComputedElementMatrix(long iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef)
{
    PZ_PTHREAD_MUTEX_LOCK(&fAccessElement,"TPZStructMatrix::ThreadData::ComputedElementMatrix()");
    std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > el(ek,ef);
    fSubmitted[iel] = el;
    fAssembly.Post();
    /*
     #ifdef MACOSX
     sem_post(fAssembly);
     #else
     sem_post(&fAssembly);
     #endif
     */
    PZ_PTHREAD_MUTEX_UNLOCK(&fAccessElement,"TPZStructMatrix::ThreadData::ComputedElementMatrix()");
}
