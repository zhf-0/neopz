//
//  TPZAssembleOptColor.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#include "TPZAssembleOptColor.h"

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

RunStatsTable stat_ass_graph_ot("-ass_graph_ot", "Run statistics table for the graph creation and coloring TPZStructMatrixOT.");

void TPZAssembleOptimizedColor::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs)
{
    if(fElSequenceColor.NElements() != fConfig->fMesh->NElements())
    {
        stat_ass_graph_ot.start();
        TPZManVector<long> ElementOrder;
        TPZStructMatrixOT::OrderElement(this->fConfig->fMesh, ElementOrder);
        TPZVec<long> elcolors;
        TPZStructMatrixOT::ElementColoring(this->fConfig->fMesh, ElementOrder, fElSequenceColor, fElBlocked, elcolors);
        stat_ass_graph_ot.stop();
    }
    
    const int numthreads = this->fConfig->fNumThreads;
    std::cout << "Assemble numthreads = " << numthreads << std::endl;
    TPZVec<pthread_t> allthreads(numthreads);
    int itr;
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZAssembleOptimizedColor::ThreadData::ThreadData()");
    pthread_cond_init(&fCondition, NULL);
    
#ifdef USING_BOOST
    this->fCurrentIndex = 0;
#endif
    
    fElementCompleted = -1;
    fElementsComputed.Resize(fConfig->fMesh->NElements());
    fElementsComputed.Fill(0);
    fSomeoneIsSleeping = 0;
    TPZManVector<ThreadData*> allthreaddata(numthreads);
#ifdef PZDEBUG
    {
        for (long i=1; i<fElBlocked.size(); i++) {
            if (fElBlocked[i] < fElBlocked[i-1]) {
                std::cout << "i = " << i << " fElBlocked[i-1] " << fElBlocked[i-1] << " fElBlocked[i] " << fElBlocked[i] << std::endl;
            }
        }
    }
#endif
    
    for(itr=0; itr<numthreads; itr++)
    {
        allthreaddata[itr] = new ThreadData(this, itr, mat, rhs);
        ThreadData &threaddata = *allthreaddata[itr];
        
        threaddata.fElBlocked=&fElBlocked;
        threaddata.fElSequenceColor=&fElSequenceColor;
        threaddata.fElementCompleted = &fElementCompleted;
        threaddata.fComputedElements = &fElementsComputed;
        threaddata.fSomeoneIsSleeping = &fSomeoneIsSleeping;
        threaddata.fCondition = &fCondition;
        threaddata.fAccessElement = &fAccessElement;
        PZ_PTHREAD_CREATE(&allthreads[itr], NULL,ThreadData::ThreadWork,
                          &threaddata, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        delete allthreaddata[itr];
    }
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"TPZAssembleOptimizedColor::ThreadData::~ThreadData()");
    pthread_cond_destroy(&fCondition);
    
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


void TPZAssembleOptimizedColor::MultiThread_Assemble(TPZFMatrix<STATE> & rhs)
{
    if(fElSequenceColor.NElements() != fConfig->fMesh->NElements())
    {
        stat_ass_graph_ot.start();
        TPZManVector<long> ElementOrder;
        TPZStructMatrixOT::OrderElement(this->fConfig->fMesh, ElementOrder);
        TPZVec<long> elcolors;
        TPZStructMatrixOT::ElementColoring(this->fConfig->fMesh, ElementOrder, fElSequenceColor, fElBlocked, elcolors);
        stat_ass_graph_ot.stop();
    }
    

    const int numthreads = this->fConfig->fNumThreads;
    TPZVec<pthread_t> allthreads(numthreads);
    int itr;
    
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZAssembleOptimizedColor::ThreadData::ThreadData()");
    pthread_cond_init(&fCondition, NULL);
    
#ifdef USING_BOOST
    this->fCurrentIndex = 0;
#endif
    fElementCompleted = -1;
    fElementsComputed.Resize(fConfig->fMesh->NElements());
    fElementsComputed.Fill(0);
    fSomeoneIsSleeping = 0;
    TPZManVector<ThreadData*> allthreaddata(numthreads);
    
    for(itr=0; itr<numthreads; itr++)
    {
        allthreaddata[itr] = new ThreadData(this, itr, rhs);
        ThreadData &threaddata = *allthreaddata[itr];
        threaddata.fElBlocked=&fElBlocked;
        threaddata.fElSequenceColor=&fElSequenceColor;
        threaddata.fElementCompleted = &fElementCompleted;
        threaddata.fComputedElements = &fElementsComputed;
        threaddata.fSomeoneIsSleeping = &fSomeoneIsSleeping;
        threaddata.fCondition = &fCondition;
        threaddata.fAccessElement = &fAccessElement;
        
        PZ_PTHREAD_CREATE(&allthreads[itr], NULL,ThreadData::ThreadWorkResidual,
                          &threaddata, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        delete allthreaddata[itr];
    }
    
    PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"TPZAssembleOptimizedColor::ThreadData::~ThreadData()");
    pthread_cond_destroy(&fCondition);
    
#ifdef LOG4CXX
    if(loggerCheck->isDebugEnabled())
    {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
#endif
    
}



TPZAssembleOptimizedColor::ThreadData::ThreadData(TPZAssembleOptimizedColor *strmat, int seqnum, TPZMatrix<STATE> &mat,
                                          TPZFMatrix<STATE> &rhs)
: fStruct(strmat), fGlobMatrix(&mat), fGlobRhs(&rhs), fThreadSeqNum(seqnum)
{
#ifdef USING_BOOST
    fCurrentIndex = &strmat->fCurrentIndex;
#endif
    
}

TPZAssembleOptimizedColor::ThreadData::ThreadData(TPZAssembleOptimizedColor *strmat, int seqnum,
                                          TPZFMatrix<STATE> &rhs)
: fStruct(strmat), fGlobMatrix(0), fGlobRhs(&rhs), fThreadSeqNum(seqnum)
{
#ifdef USING_BOOST
    this->fCurrentIndex = &strmat->fCurrentIndex;
#endif
}

TPZAssembleOptimizedColor::ThreadData::~ThreadData()
{
}

//#define DRY_RUN

void *TPZAssembleOptimizedColor::ThreadData::ThreadWork(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    //    TPZAssembleOptimizedColor *strmat = data->fStruct;
    TPZVec<long> &ComputedElements = *(data->fComputedElements);
    TPZVec<long> &ElBlocked = *(data->fElBlocked);
    
    int &SomeoneIsSleeping = *(data->fSomeoneIsSleeping);
    
    TPZCompMesh *cmesh = data->fStruct->fConfig->fMesh;
    TPZElementMatrix ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    long numelements = data->fElSequenceColor->size();
#ifdef USING_BOOST
    long index = data->fCurrentIndex->fetch_add(1);
#else
    long index = data->fThreadSeqNum;
#endif
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
        sout << "index = " << index << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef HUGEDEBUG
    tht::EnterCriticalSection(*data->fAccessElement);
    std::cout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
    std::cout << "index = " << index << std::endl;
    std::cout.flush();
    tht::LeaveCriticalSection(*data->fAccessElement);
#endif
    
#ifndef USING_BOOST
    int nthreads = data->fStruct->GetNumThreads();
    for (index = data->fThreadSeqNum; index < numelements; index += nthreads)
#else
        while (index < numelements)
#endif
        {
            
            long iel = data->fElSequenceColor->operator[](index);
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Computing element " << index;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
#ifdef LOG4CXX
            std::stringstream sout;
            sout << "Element " << index << " elapsed time ";
            TPZTimer timeforel(sout.str());
            timeforel.start();
#endif
            
            if (iel >= 0){
                TPZCompEl *el = cmesh->ElementVec()[iel];
#ifndef DRY_RUN
                el->CalcStiff(ek,ef);
#else
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(el);
                intel->InitializeElementMatrix(ek,ef);
#endif
                
                
#ifndef DRY_RUN
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
                }
#endif
                
#ifdef LOG4CXX
                timeforel.stop();
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << timeforel.processName() <<  timeforel;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
                long localcompleted = *(data->fElementCompleted);
                bool localupdated = false;
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                long needscomputed = ElBlocked[index];
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    if (localupdated) {
                        sout << "Localcompleted updated without thread lock\n";
                    }
                    
                    sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
#ifdef HUGEDEBUG
                tht::EnterCriticalSection(*data->fAccessElement);
                std::cout << "threadEK " << data->fThreadSeqNum << " index " << index << " localcompleted " << localcompleted << " needscomputed " << needscomputed << std::endl;
                tht::LeaveCriticalSection( *data->fAccessElement );
#endif
                
                bool hadtowait = false;
                while (needscomputed > localcompleted) {
                    // block the thread till the element needed has been assembled
                    tht::EnterCriticalSection(*data->fAccessElement);
                    SomeoneIsSleeping = 1;
                    hadtowait = true;
#ifdef HUGEDEBUG
                    std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " going to sleep waiting for " << needscomputed << std::endl;
                    std::cout.flush();
#endif
#ifdef LOG4CXX
                    if (logger->isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Element " << index << " cannot be assembled - going to sleep";
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    pthread_cond_wait(data->fCondition, data->fAccessElement);
                    tht::LeaveCriticalSection( *data->fAccessElement );
                    
                    localcompleted = *data->fElementCompleted;
                    localupdated = false;
                    while (ComputedElements[localcompleted+1] == 1) {
                        localcompleted++;
                        localupdated = true;
                    }
#ifdef LOG4CXX
                    if (logger->isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "thread wakeup for element index " << index << std::endl;
                        if (localupdated) {
                            sout << "Localcompleted updated without thread lock\n";
                        }
                        
                        sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
                
#ifdef HUGEDEBUG
                if (hadtowait) {
                    tht::EnterCriticalSection(*data->fAccessElement);
                    std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " continuing\n";
                    tht::LeaveCriticalSection( *data->fAccessElement );
                }
#endif
                
#ifndef DRY_RUN
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
#endif
                
                localupdated = false;
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                if (localcompleted == index-1) {
                    localcompleted++;
                }
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                bool elementcompletedupdate = false;
                if (*data->fElementCompleted < localcompleted) {
                    //                std::cout << "Updating element completed " << localcompleted << std::endl;
                    *data->fElementCompleted = localcompleted;
                    elementcompletedupdate = true;
                }
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element " << index << " has been assembled ";
                    if (localupdated) {
                        sout << "\nLocalcompleted updated without thread lock\n";
                    }
                    if (elementcompletedupdate) {
                        sout << "\nfElementCompleted updated to localcompleted\n";
                    }
                    sout << "local completed " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                ComputedElements[index] = 1;
                if (SomeoneIsSleeping) {
                    tht::EnterCriticalSection( *data->fAccessElement );
#ifdef HUGEDEBUG
                    std::cout << "threadEK " <<data->fThreadSeqNum <<  " Computed index " << index << " Waking up ElementsCompleted " << *data->fElementCompleted << std::endl;
                    std::cout.flush();
#endif
                    
#ifdef LOG4CXX
                    if (logger->isDebugEnabled())
                    {
                        LOGPZ_DEBUG(logger, "condition broadcast")
                    }
#endif
                    SomeoneIsSleeping = 0;
                    pthread_cond_broadcast(data->fCondition);
                    tht::LeaveCriticalSection( *data->fAccessElement );
                }
                
            }
            else
            {
                std::cout << "the element in ElColorSequence is negative???\n";
                DebugStop();
            }
#ifdef USING_BOOST
            index = data->fCurrentIndex->fetch_add(1);
#endif
            
        }
    // just make sure threads that were accidentally blocked get woken up
    tht::EnterCriticalSection( *data->fAccessElement );
    bool localupdated = false;
    long localcompleted = *(data->fElementCompleted);
    while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
        localcompleted++;
        localupdated = true;
    }
    if (localcompleted == index-1) {
        localcompleted++;
    }
    bool elementcompletedupdate = false;
    if (*data->fElementCompleted < localcompleted) {
        //                std::cout << "Updating element completed " << localcompleted << std::endl;
        *data->fElementCompleted = localcompleted;
        elementcompletedupdate = true;
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        LOGPZ_DEBUG(logger, "Finishing up")
        if (localupdated) {
            LOGPZ_DEBUG(logger, "updated localcompleted")
        }
        if (elementcompletedupdate) {
            LOGPZ_DEBUG(logger, "updated fElementCompleted")
        }
        if (localupdated || elementcompletedupdate) {
            std::stringstream sout;
            sout << "localcompleted " << localcompleted;
            LOGPZ_DEBUG(logger, sout.str())
        }
        LOGPZ_DEBUG(logger, "finishing and condition broadcast")
    }
#endif
    pthread_cond_broadcast(data->fCondition);
    SomeoneIsSleeping = 0;
    tht::LeaveCriticalSection( *data->fAccessElement );
    return 0;
}

void *TPZAssembleOptimizedColor::ThreadData::ThreadWorkResidual(void *datavoid)
{
#ifdef LOG4CXX
    logger->setLevel(log4cxx::Level::getInfo());
#endif
    ThreadData *data = (ThreadData *) datavoid;
    //    TPZAssembleOptimizedColor *strmat = data->fStruct;
    TPZVec<long> &ComputedElements = *(data->fComputedElements);
    TPZVec<long> &ElBlocked = *(data->fElBlocked);
    
    int &SomeoneIsSleeping = *(data->fSomeoneIsSleeping);
    
    TPZCompMesh *cmesh = data->fStruct->fConfig->fMesh;
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    long numelements = data->fElSequenceColor->size();
#ifdef USING_BOOST
    long index = data->fCurrentIndex->fetch_add(1);
#else
    long index = data->fThreadSeqNum;
#endif
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
        sout << "index = " << index << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef HUGEDEBUG
    tht::EnterCriticalSection(*data->fAccessElement);
    std::cout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
    std::cout << "index = " << index << std::endl;
    std::cout.flush();
    tht::LeaveCriticalSection(*data->fAccessElement);
#endif
    
#ifndef USING_BOOST
    int nthreads = data->fStruct->GetNumThreads();
    for (index = data->fThreadSeqNum; index < numelements; index += nthreads)
#else
        while (index < numelements)
#endif
        {
            
            long iel = data->fElSequenceColor->operator[](index);
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Computing element " << index;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
#ifdef LOG4CXX
            std::stringstream sout;
            sout << "Element " << index << " elapsed time ";
            TPZTimer timeforel(sout.str());
            timeforel.start();
#endif
            
            if (iel >= 0){
                TPZCompEl *el = cmesh->ElementVec()[iel];
#ifndef DRY_RUN
                el->CalcResidual(ef);
#else
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(el);
                intel->InitializeElementMatrix(ef);
#endif
                
#ifndef DRY_RUN
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
                }
#endif
                
#ifdef LOG4CXX
                timeforel.stop();
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << timeforel.processName() <<  timeforel;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
                long localcompleted = *(data->fElementCompleted);
                bool localupdated = false;
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                long needscomputed = ElBlocked[index];
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    if (localupdated) {
                        sout << "Localcompleted updated without thread lock\n";
                    }
                    
                    sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
#ifdef HUGEDEBUG
                tht::EnterCriticalSection(*data->fAccessElement);
                std::cout << "threadEK " << data->fThreadSeqNum << " index " << index << " localcompleted " << localcompleted << " needscomputed " << needscomputed << std::endl;
                tht::LeaveCriticalSection( *data->fAccessElement );
#endif
                
                bool hadtowait = false;
#ifdef LOG4CXX
                if (logger->isInfoEnabled() && needscomputed > localcompleted)
                {
                    std::stringstream sout;
                    sout << "Element " << index << " cannot be assembled - going to sleep";
                    LOGPZ_INFO(logger, sout.str())
                }
#endif
                while (needscomputed > localcompleted) {
                    // block the thread till the element needed has been assembled
                    tht::EnterCriticalSection(*data->fAccessElement);
                    SomeoneIsSleeping = 1;
                    hadtowait = true;
#ifdef HUGEDEBUG
                    std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " going to sleep waiting for " << needscomputed << std::endl;
                    std::cout.flush();
#endif
#ifdef LOG4CXX
                    if (logger->isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Element " << index << " cannot be assembled - going to sleep";
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    pthread_cond_wait(data->fCondition, data->fAccessElement);
                    tht::LeaveCriticalSection( *data->fAccessElement );
                    
                    localcompleted = *data->fElementCompleted;
                    localupdated = false;
                    while (ComputedElements[localcompleted+1] == 1) {
                        localcompleted++;
                        localupdated = true;
                    }
#ifdef LOG4CXX
                    if (logger->isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "thread wakeup for element index " << index << std::endl;
                        if (localupdated) {
                            sout << "Localcompleted updated without thread lock\n";
                        }
                        
                        sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
                
                
#ifdef LOG4CXX
                if (logger->isInfoEnabled() && hadtowait)
                {
                    std::stringstream sout;
                    sout << "thread wakeup for element index " << index << std::endl;
                    LOGPZ_INFO(logger, sout.str())
                }
#endif
                
                
#ifdef HUGEDEBUG
                if (hadtowait) {
                    tht::EnterCriticalSection(*data->fAccessElement);
                    std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " continuing\n";
                    tht::LeaveCriticalSection( *data->fAccessElement );
                }
#endif
                
#ifndef DRY_RUN
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
#endif
                
                localupdated = false;
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                if (localcompleted == index-1) {
                    localcompleted++;
                }
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                bool elementcompletedupdate = false;
                if (*data->fElementCompleted < localcompleted) {
                    //                std::cout << "Updating element completed " << localcompleted << std::endl;
                    *data->fElementCompleted = localcompleted;
                    elementcompletedupdate = true;
                }
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element " << index << " has been assembled ";
                    if (localupdated) {
                        sout << "\nLocalcompleted updated without thread lock\n";
                    }
                    if (elementcompletedupdate) {
                        sout << "\nfElementCompleted updated to localcompleted\n";
                    }
                    sout << "local completed " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                ComputedElements[index] = 1;
                if (SomeoneIsSleeping) {
                    tht::EnterCriticalSection( *data->fAccessElement );
#ifdef HUGEDEBUG
                    std::cout << "threadEK " <<data->fThreadSeqNum <<  " Computed index " << index << " Waking up ElementsCompleted " << *data->fElementCompleted << std::endl;
                    std::cout.flush();
#endif
                    
#ifdef LOG4CXX
                    if (logger->isDebugEnabled())
                    {
                        LOGPZ_DEBUG(logger, "condition broadcast")
                    }
#endif
                    SomeoneIsSleeping = 0;
                    pthread_cond_broadcast(data->fCondition);
                    tht::LeaveCriticalSection( *data->fAccessElement );
                }
                
            }
            else
            {
                std::cout << "the element in ElColorSequence is negative???\n";
                DebugStop();
            }
#ifdef USING_BOOST
            index = data->fCurrentIndex->fetch_add(1);
#endif
            
        }
    // just make sure threads that were accidentally blocked get woken up
    tht::EnterCriticalSection( *data->fAccessElement );
    bool localupdated = false;
    long localcompleted = *(data->fElementCompleted);
    while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
        localcompleted++;
        localupdated = true;
    }
    if (localcompleted == index-1) {
        localcompleted++;
    }
    bool elementcompletedupdate = false;
    if (*data->fElementCompleted < localcompleted) {
        //                std::cout << "Updating element completed " << localcompleted << std::endl;
        *data->fElementCompleted = localcompleted;
        elementcompletedupdate = true;
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        LOGPZ_DEBUG(logger, "Finishing up")
        if (localupdated) {
            LOGPZ_DEBUG(logger, "updated localcompleted")
        }
        if (elementcompletedupdate) {
            LOGPZ_DEBUG(logger, "updated fElementCompleted")
        }
        if (localupdated || elementcompletedupdate) {
            std::stringstream sout;
            sout << "localcompleted " << localcompleted;
            LOGPZ_DEBUG(logger, sout.str())
        }
        LOGPZ_DEBUG(logger, "finishing and condition broadcast")
    }
#endif
    pthread_cond_broadcast(data->fCondition);
    SomeoneIsSleeping = 0;
    tht::LeaveCriticalSection( *data->fAccessElement );
    return 0;
}

//static bool CanAssemble(TPZStack<long> &connectlist, TPZVec<int> &elContribute)
//{
//    for (int i = 0 ; i < connectlist.NElements() ; i++)
//    {
//        if (elContribute[connectlist[i]] >= 0){
//            return false;
//        }
//    }
//    return true;
//}

static void AssembleColor(int el,TPZStack<long> &connectlist, TPZVec<int> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        elContribute[connectlist[i]] = el;
    }
}

static int WhoBlockedMe(TPZStack<long> &connectlist, TPZVec<int> &elContribute, TPZVec<int> &elSeqinv)
{
    int el = -1;
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        int elBlocked = elContribute[connectlist[i]];
        if (elBlocked == -1) continue;
        int elBlockedIndex = elSeqinv[elBlocked];
        if (el == -1) el = elBlockedIndex;
        if (elBlockedIndex > el) el = elBlockedIndex;
    }
    return el;
}

//static void RemoveEl(int el,TPZCompMesh *cmesh,TPZVec<int> &elContribute,int elSequence)
//{
//    TPZCompEl *cel = cmesh->ElementVec()[el];
//    if(!cel) DebugStop();
//    TPZStack<long> connectlist;
//    cel->BuildConnectList(connectlist);
//    for (int i = 0 ; i < connectlist.NElements() ; i++)
//    {
//        int conindex = connectlist[i];
//        if (elContribute[conindex] != elSequence){
//            DebugStop();
//        }
//        elContribute[conindex] = -1;
//    }
//}

static int MinPassIndex(TPZStack<long> &connectlist,TPZVec<int> &elContribute, TPZVec<int> &passIndex)
{
    int minPassIndex = -1;
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        int elcont = elContribute[connectlist[i]];
        int passindex = -1;
        if (elcont != -1){
            passindex = passIndex[elcont];
            if (minPassIndex == -1) minPassIndex = passindex;
        }
        if (minPassIndex < passindex) minPassIndex = passindex;
    }
    return minPassIndex;
}

// elSequence (input) element sequence acording to the connect sequence numbers
// elSequenceColor (output) the coloured element sequence
// elBlocked the element index which needs to have been computed before assembling the element
// elColors (output) number of elements in each color
void TPZAssembleOptimizedColor::ElementColoring(TPZCompMesh *cmesh, TPZVec<long> &elSequence, TPZVec<long> &elSequenceColor,
                                        TPZVec<long> &elBlocked, TPZVec<long> &NumelColors)
{
    
    const int nnodes = cmesh->NConnects();
    const int nel = cmesh->ElementVec().NElements();
    
    if (!nel) return;
    
    NumelColors.Resize(nel, -1);
    
    // elContribute contains the element index which last contributed to the node
    // passIndex essentially contains the color of the element (?)
    TPZManVector<int> elContribute(nnodes,-1), passIndex(nel,-1), elSequenceColorInv(nel,-1);
    long elsequencesize = elSequence.size();
    elSequenceColor.Resize(elsequencesize);
    elSequenceColor.Fill(-1);
    elBlocked.Resize(elsequencesize);
    elBlocked.Fill(-1);
    int nelProcessed = 0;
    int currentEl = 0;
    int currentPassIndex = 0;
    while (nelProcessed < elSequence.NElements()){
        
        int elindex = elSequence[currentEl];
        
        // if this element hasn t been computed in a previous pass
        if(elSequenceColorInv[elindex] == -1)
        {
            TPZCompEl *cel = cmesh->ElementVec()[elindex];
            
            
            if(!cel) continue;
            TPZStack<long> connectlist;
            cel->BuildConnectList(connectlist);
            //      std::cout << "elcontribute " << elContribute << std::endl;
            //      std::cout << "connectlist " << connectlist << std::endl;
            
            // compute the lowest color (pass index) of the elements that have contributed to this set of nodes
            int minPass = MinPassIndex(connectlist,elContribute,passIndex);
            // no element has ever seen any of these nodes
            if (minPass == -1){
                passIndex[elindex] = currentPassIndex;
                // the element index is put into elContribute (from there it is possible to know the colour as well)
                AssembleColor(elindex,connectlist,elContribute);
                // initialize the data structures
                elSequenceColor[nelProcessed] = elindex;
                elSequenceColorInv[elindex] = nelProcessed;
                nelProcessed++;
            }
            // this element cannot be computed as it is connected to another element of the current colour
            else if (minPass == currentPassIndex){
            }
            // the elements connected to this node are from a previous colour
            else if (minPass < currentPassIndex){
                // the element with largest index which contributes to the set of nodes
                // el is given in the new sequence order
                const int el = WhoBlockedMe(connectlist,elContribute, elSequenceColorInv);
                // elblocked means the future element the element el will block
                if (elBlocked[nelProcessed] != -1) DebugStop();
                elBlocked[nelProcessed] = el;
                passIndex[elindex] = currentPassIndex;
                AssembleColor(elindex,connectlist,elContribute);
                elSequenceColor[nelProcessed] = elindex;
                elSequenceColorInv[elindex] = nelProcessed;
                nelProcessed++;
            }
            else{
                DebugStop();
            }
        }
        currentEl++;
        if (currentEl == elSequence.NElements()){
            currentEl = 0;
            NumelColors[currentPassIndex]=nelProcessed;
            currentPassIndex++;
        }
    }
    
    NumelColors[currentPassIndex]=NumelColors[currentPassIndex-1]+1;
    
    NumelColors.Resize(currentPassIndex+1);
    
    
#ifdef PZDEBUG
    std::ofstream toto("../ColorMeshDebug.txt");
    toto << "elSequence\n" << elSequence << std::endl;
    toto << "elSequenceColor\n" << elSequenceColor << std::endl;
    toto << "elSequenceColorInv\n" << elSequenceColorInv << std::endl;
    toto << "elBlocked\n" << elBlocked << std::endl;
    toto << "elContribute\n" << elContribute << std::endl;
    toto << "passIndex\n" << passIndex << std::endl;
    toto << "NumelColors\n" << NumelColors << std::endl;
    toto.close();
#endif
}

