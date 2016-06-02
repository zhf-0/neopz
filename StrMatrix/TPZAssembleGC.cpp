//
//  TPZAssemble.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#include "TPZAssembleGC.h"

#include "run_stats_table.h"
#include "TPZTimer.h"

#include "pzlog.h"
#include "pzsubcmesh.h"
#include "pzmaterial.h"
#include "TPZThreadTools.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZAsemble"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
static LoggerPtr loggerGlobStiff(Logger::getLogger("pz.strmatrix.globalstiffness"));
#endif


RunStatsTable stat_ass_graph("-ass_graph", "Run statistics table for the graph creation and coloring TPZAssembleGlobalColor.");



void TPZAssembleGlobalColor::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs)
{
    
    stat_ass_graph.start();
    if (felSequenceColor.NElements() != fConfig->fMesh->NElements())
    {
        TPZManVector<long> ElementOrder;
        OrderElement(fConfig->fMesh, ElementOrder);
        ElementColoring(fConfig->fMesh, ElementOrder, felSequenceColor, fnextBlocked);
        stat_ass_graph.stop();
    }
    ThreadData threaddata(this,mat,rhs);
    
    threaddata.fnextBlocked=&fnextBlocked;
    threaddata.felSequenceColor=&felSequenceColor;
    
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


void TPZAssembleGlobalColor::MultiThread_Assemble(TPZFMatrix<STATE> & rhs)
{
    if (felSequenceColor.NElements() != fConfig->fMesh->NElements())
    {
        TPZManVector<long> ElementOrder;
        OrderElement(fConfig->fMesh, ElementOrder);
        ElementColoring(fConfig->fMesh, ElementOrder, felSequenceColor, fnextBlocked);
        stat_ass_graph.stop();
    }

    ThreadData threaddata(this,rhs);
    
    threaddata.fnextBlocked=&fnextBlocked;
    threaddata.felSequenceColor=&felSequenceColor;
    
    const int numthreads = this->fConfig->fNumThreads;
    TPZVec<pthread_t> allthreads(numthreads);
    int itr;
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_CREATE(&allthreads[itr], NULL,ThreadData::ThreadWorkResidual,
                          &threaddata, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
    }
}



TPZAssembleGlobalColor::ThreadData::ThreadData(TPZAssembleGlobalColor *strmat, TPZMatrix<STATE> &mat,
                                          TPZFMatrix<STATE> &rhs)
: fStruct(strmat), fGlobMatrix(&mat), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZAssembleGlobalColor::ThreadData::ThreadData()");
    pthread_cond_init(&fCondition, NULL);
}

TPZAssembleGlobalColor::ThreadData::ThreadData(TPZAssembleGlobalColor *strmat,
                                          TPZFMatrix<STATE> &rhs)
: fStruct(strmat), fGlobMatrix(0), fGlobRhs(&rhs), fNextElement(0)
{
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZAssembleGlobalColor::ThreadData::ThreadData()");
    pthread_cond_init(&fCondition, NULL);
}

TPZAssembleGlobalColor::ThreadData::~ThreadData()
{
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"TPZAssembleGlobalColor::ThreadData::~ThreadData()");
    pthread_cond_destroy(&fCondition);
}

void *TPZAssembleGlobalColor::ThreadData::ThreadWork(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    TPZCompMesh *cmesh = data->fStruct->fConfig->fMesh;
    int nel = cmesh->NElements();
    bool hasWork = true;
    int iel = 0;
    TPZElementMatrix ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    while(hasWork)
    {
        tht::EnterCriticalSection(data->fAccessElement);
        // nextelement is a protected value
        long localiel = data->fNextElement;
        long blockedel = -1;
        if(data->felBlocked.size()) blockedel = data->felBlocked.begin()->first;
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Checking out " << localiel << " next blocked element " << blockedel;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        if(localiel < nel)
        {
            // The elblocked data structure indicates the elements blocked by the elements being processed (why is this needed?)
            if (!data->felBlocked.size() || data->felBlocked.begin()->first > localiel)
            {
                
                if(localiel==-1) DebugStop();
                // this is the next element that will be computed
                iel = (*data->felSequenceColor)[localiel];
                
                // identify the element which will be blocked by iel
                int elBl = (*data->fnextBlocked)[localiel];
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element can be computed iel = " << localiel << " element blocked " << elBl;
                    LOGPZ_DEBUG(logger, sout.str())
                }
                
#endif
                // update the datastructure with the highest element that can be processed while iel hasn't been computed yet
                if (elBl >= 0){
                    data->felBlocked[elBl]++;
                    hasWork = true;
                }
                // the next thread will monitor a higher numbered element
                data->fNextElement++;
            }
            else if (data->felBlocked.size() || data->felBlocked.begin()->first <= localiel){
#ifdef LOG4CXX
                if (logger->isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "Going to sleep cannot do " << localiel << " because of " << data->felBlocked.begin()->first << " has to be computed first";
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                // localiel cannot be processed yet
                data->fSleeping=true;
                while(data->fSleeping){
                    pthread_cond_wait(&data->fCondition, &data->fAccessElement);
                }
                iel = -1;
                hasWork = true;
            }
            else{
                DebugStop();
            }
        }
        else{
            hasWork = false;
            iel = -1;
        }
        
        tht::LeaveCriticalSection(data->fAccessElement);
        
#ifdef LOG4CXX
        std::stringstream sout;
        sout << "Element " << localiel << " elapsed time ";
        TPZTimer timeforel(sout.str());
        timeforel.start();
#endif
        
        if (iel >= 0){
            TPZCompEl *el = cmesh->ElementVec()[iel];
            el->CalcStiff(ek,ef);
            if(!el->HasDependency()) {
                
                ek.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                
#ifdef LOG4CXX
                if(loggerel->isDebugEnabled())
                {
                    std::stringstream sout;
                    ek.fMat.Print("Element stiffness matrix",sout);
                    ef.fMat.Print("Element right hand side", sout);
                    LOGPZ_DEBUG(loggerel,sout.str())
                }
#endif
            } else {
                // the element has dependent nodes
                ek.ApplyConstraints();
                ef.ApplyConstraints();
                ek.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                
#ifdef LOG4CXX
                if(loggerel2->isDebugEnabled() && el->Reference() &&  el->Reference()->MaterialId() == 1 && el->IsInterface())
                {
                    std::stringstream sout;
                    el->Reference()->Print(sout);
                    el->Print(sout);
                    ek.Print(sout);
                    ef.Print(sout);
                    LOGPZ_DEBUG(loggerel2,sout.str())
                }
#endif
            }
            
            
            if(data->fGlobMatrix){
                // Assemble the matrix
                if(!ek.HasDependency())
                {
                    data->fGlobMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                    data->fGlobRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
                else
                {
                    data->fGlobMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                    data->fGlobRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
                
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
            
            tht::EnterCriticalSection( data->fAccessElement );
            
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Computed and Assembled " << localiel;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            // clean up the data structure
            int elBl = (*data->fnextBlocked)[localiel];
            if (elBl >= 0 && data->felBlocked.find(elBl) != data->felBlocked.end())
            {
                data->felBlocked[elBl]--;
                int dataelbl = data->felBlocked[elBl];
                if (dataelbl == 0){
#ifdef LOG4CXX
                    if(logger->isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Element " << elBl << " released";
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    data->felBlocked.erase(elBl);
                    if(data->fSleeping) {
                        data->fSleeping=false;
                        
#ifdef LOG4CXX
                        if(logger->isDebugEnabled()){
                            LOGPZ_DEBUG(logger, "Waking up everybody")
                        }
#endif
                        // wake up everybody
                        pthread_cond_broadcast(&data->fCondition);
                    }
                }
                else if (dataelbl < 0)
                {
                    DebugStop();
                }
                else
                {
#ifdef LOG4CXX
                    if (logger->isDebugEnabled()) {
                        std::stringstream sout;
                        sout << "Not freeing the blocked element " << elBl;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
            }
            else if(elBl >= 0){
                DebugStop();
            }
            
            tht::LeaveCriticalSection( data->fAccessElement );
        }
    }
    return 0;
}

void *TPZAssembleGlobalColor::ThreadData::ThreadWorkResidual(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    TPZCompMesh *cmesh = data->fStruct->fConfig->fMesh;
    int nel = cmesh->NElements();
    bool hasWork = true;
    int iel = 0;
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    while(hasWork)
    {
        tht::EnterCriticalSection(data->fAccessElement);
        int localiel = data->fNextElement;
        if(data->fNextElement < nel){
            if (!data->felBlocked.size() || data->felBlocked.begin()->first > localiel){
                iel = (*data->felSequenceColor)[localiel];
                
                if(localiel==-1) DebugStop();
                
                int elBl = (*data->fnextBlocked)[localiel];
                if (elBl >= 0){
                    data->felBlocked[elBl]++;
                    hasWork = true;
                }
                
                data->fNextElement++;
            }
            else if (data->felBlocked.size() || data->felBlocked.begin()->first <= localiel){
                data->fSleeping=true;
                while(data->fSleeping){
                    pthread_cond_wait(&data->fCondition, &data->fAccessElement);
                }
                iel = -1;
                hasWork = true;
            }
            else{
                DebugStop();
            }
        }
        else{
            hasWork = false;
            iel = -1;
        }
        
        tht::LeaveCriticalSection(data->fAccessElement);
        if (iel >= 0){
            TPZCompEl *el = cmesh->ElementVec()[iel];
            
            el->CalcResidual(ef);
            
            if(!el->HasDependency()) {
                ef.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                
#ifdef LOG4CXX
                if(loggerel->isDebugEnabled())
                {
                    std::stringstream sout;
                    ef.fMat.Print("Element right hand side", sout);
                    LOGPZ_DEBUG(loggerel,sout.str())
                }
#endif
            } else {
                // the element has dependent nodes
                ef.ApplyConstraints();
                ef.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
#ifdef LOG4CXX
                if(loggerel2->isDebugEnabled() && el->Reference() &&  el->Reference()->MaterialId() == 1 && el->IsInterface())
                {
                    std::stringstream sout;
                    el->Reference()->Print(sout);
                    el->Print(sout);
                    LOGPZ_DEBUG(loggerel2,sout.str())
                }
#endif
            }
            
            if(data->fGlobRhs){
                // Assemble the matrix
                if(!ef.HasDependency())
                {
                    data->fGlobRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
                }
                else
                {
                    data->fGlobRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
                }
            }
            tht::EnterCriticalSection( data->fAccessElement );
            int elBl = (*data->fnextBlocked)[localiel];
            if (elBl >= 0 && data->felBlocked.find(elBl) != data->felBlocked.end()){
                data->felBlocked[elBl]--;
                if (data->felBlocked[elBl] == 0){
                    data->felBlocked.erase(elBl);
                    if(data->fSleeping) {
                        data->fSleeping=false;
                    }
                    pthread_cond_broadcast(&data->fCondition);
                }
            }
            else if(elBl >= 0){
                DebugStop();
            }
            tht::LeaveCriticalSection(data->fAccessElement);
        }
    }
    
    return 0;
}

