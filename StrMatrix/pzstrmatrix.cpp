/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrix methods.
 */

#include "pzstrmatrix.h"

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"
#include "pzsfulmat.h"

#include "pzgnode.h"
#include "TPZTimer.h"
#include "TPZThreadTools.h"

#include "TPZAssemble.h"
#include "TPZAssemblePThread.h"

#include "pzcheckconsistency.h"
#include "pzmaterial.h"
#include "run_stats_table.h"

using namespace std;

#include "pzlog.h"

#include "pz_pthread.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZStructMatrix"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
static LoggerPtr loggerGlobStiff(Logger::getLogger("pz.strmatrix.globalstiffness"));
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


TPZStructMatrix::TPZStructMatrix(TPZCompMesh *mesh) : fAssembleConfig()
{
    fAssembleConfig.fMesh = mesh;
    long neq = mesh->NEquations();
    fAssembleConfig.fEquationFilter.SetNumEq(neq);
    this->SetNumThreads(0);
    fStiffAssemble = new TPZAssemblePThread(fAssembleConfig);
    fRhsAssemble = new TPZAssemblePThread(fAssembleConfig);
}

TPZStructMatrix::TPZStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh) : fAssembleConfig()
{
    fAssembleConfig.fMesh = cmesh.operator->();
    fAssembleConfig.fCompMesh = cmesh;
    fAssembleConfig.fEquationFilter.SetNumEq(cmesh->NEquations());
    this->SetNumThreads(0);
}

TPZStructMatrix::TPZStructMatrix(const TPZStructMatrix &copy) : fAssembleConfig(copy.fAssembleConfig)
{
    if (!fStiffAssemble || ! fRhsAssemble) {
        DebugStop();
    }
    fStiffAssemble = fStiffAssemble->Clone(fAssembleConfig);
    fRhsAssemble = fRhsAssemble->Clone(fAssembleConfig);
}

TPZStructMatrix::~TPZStructMatrix(){
    delete fStiffAssemble;
    delete fRhsAssemble;
}



TPZMatrix<STATE> *TPZStructMatrix::Create() {
    cout << "TPZStructMatrix::Create should never be called\n";
    return 0;
}

TPZStructMatrix *TPZStructMatrix::Clone() {
    cout << "TPZStructMatrix::Clone should never be called\n";
    DebugStop();
    return 0;
}


void TPZStructMatrix::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs)
{
    fStiffAssemble->Assemble(stiffness,rhs);
}

void TPZStructMatrix::Assemble(TPZFMatrix<STATE> & rhs){
    fRhsAssemble->Assemble(rhs);
}


TPZMatrix<STATE> * TPZStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs)
{
    TPZMatrix<STATE> *stiff = Create();
    
    long cols = MAX(1, rhs.Cols());
    rhs.Redim(fAssembleConfig.fEquationFilter.NEqExpand(),cols);
    Assemble(*stiff,rhs);
    
#ifdef LOG4CXX2
    if(loggerel->isDebugEnabled())
    {
        std::stringstream sout;
        stiff->Print("Stiffness matrix",sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel,sout.str())
    }
#endif
    return stiff;
    
}

/// Set the set of material ids which will be considered when assembling the system
void TPZStructMatrix::SetMaterialIds(const std::set<int> &materialids)
{
    fAssembleConfig.fMaterialIds = materialids;
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::set<int>::const_iterator it;
        std::stringstream sout;
        sout << "setting input material ids ";
        for(it=materialids.begin(); it!= materialids.end(); it++)
        {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    if(!fAssembleConfig.fMesh)
    {
        LOGPZ_WARN(logger,"SetMaterialIds called without mesh")
        return;
    }
    long iel;
    TPZAdmChunkVector<TPZCompEl*> &elvec = fAssembleConfig.fMesh->ElementVec();
    long nel = elvec.NElements();
    for(iel=0; iel<nel; iel++)
    {
        TPZCompEl *cel = elvec[iel];
        if(!cel) continue;
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (cel);
        if(!subcmesh) continue;
        TPZAutoPointer<TPZAnalysis> anal = subcmesh->Analysis();
        if(!anal)
        {
            LOGPZ_ERROR(logger,"SetMaterialIds called for substructure without analysis object")
            DebugStop();
        }
        TPZAutoPointer<TPZStructMatrix> str = anal->StructMatrix();
        if(!str)
        {
            LOGPZ_WARN(logger,"SetMaterialIds called for substructure without structural matrix")
            continue;
        }
        str->SetMaterialIds(materialids);
    }
}


