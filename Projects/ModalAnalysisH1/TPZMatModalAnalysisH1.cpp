#include "TPZMatModalAnalysisH1.h"
#include "TPZMatHCurlProjection.h"
#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif




TPZMatModalAnalysisH1::TPZMatModalAnalysisH1(int id, REAL freq, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &) ) :
TPZMaterial(id), fUr(ur), fEr(er)
{
    assembling = NDefined;
    fW = 2.*M_PI*freq;
    fSolvingForTE = true;
}

TPZMatModalAnalysisH1::TPZMatModalAnalysisH1(int id) : TPZMaterial(id), fUr(urDefault),
fEr(erDefault)
{
    assembling = NDefined;
    fW = 2.*M_PI*1e+9;
    fSolvingForTE = true;
}

/** @brief Default constructor */
TPZMatModalAnalysisH1::TPZMatModalAnalysisH1() : TPZMaterial(), fUr(urDefault),
fEr(erDefault)
{
    assembling = NDefined;
    fW=2.*M_PI*1e+9;
    fSolvingForTE = true;
}


TPZMatModalAnalysisH1::TPZMatModalAnalysisH1(const TPZMatModalAnalysisH1 &mat) : TPZMaterial(mat), fUr(mat.fUr),
fEr(mat.fEr)
{
    assembling = NDefined;
    fW = mat.fW;
    fSolvingForTE = mat.fSolvingForTE;
}

TPZMatModalAnalysisH1::~TPZMatModalAnalysisH1()
{
    
}


void TPZMatModalAnalysisH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatModalAnalysisH1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

void TPZMatModalAnalysisH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatModalAnalysisH1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

int TPZMatModalAnalysisH1::IntegrationRuleOrder(int elPMaxOrder) const
{
    return 2+elPMaxOrder*2;
}


int TPZMatModalAnalysisH1::VariableIndex(const std::string &name)
{
    if( strcmp(name.c_str(), "Ez") == 0) return 0;
    DebugStop();
    return 1;
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatModalAnalysisH1::NSolutionVariables(int var)
{
    switch (var) {
        case 0: //Ez
            return 1;
            break;
        default:
            DebugStop();
            break;
    }
    return 1;
}





