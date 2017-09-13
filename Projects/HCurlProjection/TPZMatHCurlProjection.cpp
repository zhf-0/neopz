#include "TPZMatHCurlProjection.h"

#include "pzbndcond.h"
#include "pzlog.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif



TPZMatHCurlProjection::TPZMatHCurlProjection(int id) : TPZVecL2(id)
{

}

/** @brief Default constructor */
TPZMatHCurlProjection::TPZMatHCurlProjection() : TPZVecL2()
{

}


TPZMatHCurlProjection::TPZMatHCurlProjection(const TPZMatHCurlProjection &mat) : TPZVecL2(mat)
{
    
}

TPZMatHCurlProjection::~TPZMatHCurlProjection()
{
    
}

void TPZMatHCurlProjection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    /*********************CREATE HCURL FUNCTIONS*****************************/
    TPZFNMatrix<36, REAL> phiHCurlAxes = data.phi;
    
    TPZFNMatrix<40, REAL> phiHCurl;
    
    TPZAxesTools<REAL>::Axes2XYZ(phiHCurlAxes, phiHCurl, data.axes, false);
    
    TPZManVector<REAL, 3> ax1(3), ax2(3), elNormal(3);
    for (int i = 0; i < 3; i++) {
        ax1[i] = data.axes(0, i);
        ax2[i] = data.axes(1, i);
    }
    
    //*****************GET FORCING FUNCTION****************//
    TPZManVector<STATE, 3> force(3);
    force.Fill(0.);
    if (fForcingFunction) {
        fForcingFunction->Execute(data.x, force);
    } else {
        DebugStop(); // RHS not set!
    }
    
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHCurlFunctions = phiHCurl.Rows();
    
    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        STATE load = 0.;
        load += phiHCurl(iVec, 0) * force[0];
        load += phiHCurl(iVec, 1) * force[1];
        load += phiHCurl(iVec, 2) * force[2];
        ef(iVec) += load * weight;
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE stiff = 0.;
            
            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += phiHCurl(iVec, 0) * phiHCurl(jVec, 0);
            phiIdotPhiJ += phiHCurl(iVec, 1) * phiHCurl(jVec, 1);
            phiIdotPhiJ += phiHCurl(iVec, 2) * phiHCurl(jVec, 2);
            
            stiff = phiIdotPhiJ;
            ek(iVec, jVec) += stiff * weight;
        }
    }
}

void TPZMatHCurlProjection::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatHCurlProjection::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatHCurlProjection::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    TPZFMatrix<REAL> &phiHCurl = data.phi;
    
    int nHCurlFunctions = phiHCurl.Rows();
    REAL BIG = TPZMaterial::gBigNumber;
    
    // const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de
    // condicao mista
    const STATE v2 = bc.Val2()(0, 0); // sera posto no vetor F
    
    switch (bc.Type()) {
        case 0:
            for (int i = 0; i < nHCurlFunctions; i++) {
                const STATE rhs = phiHCurl(i, 0) * BIG * v2;
                ef(i, 0) += rhs * weight;
                for (int j = 0; j < nHCurlFunctions; j++) {
                    const STATE stiff = phiHCurl(i, 0) * phiHCurl(j, 0) * BIG;
                    ek(i, j) += stiff * weight;
                }
            }
            break;
        case 1:
            DebugStop();
            break;
        case 2:
            DebugStop();
            break;
    }
}

void TPZMatHCurlProjection::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();

}

void TPZMatHCurlProjection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatHCurlProjection::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

int TPZMatHCurlProjection::IntegrationRuleOrder(int elPMaxOrder) const
{
    return 2+elPMaxOrder*2;
}


int TPZMatHCurlProjection::VariableIndex(const std::string &name)
{
    if (strcmp(name.c_str(), "E") == 0)
        return 0;
    if (strcmp(name.c_str(), "curlE") == 0)
        return 1;
    DebugStop();
    return 1;
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatHCurlProjection::NSolutionVariables(int var)
{
    switch (var) {
    case 0: // E
        return 2;
        break;
    case 1: // curlE
        return 2;
        break;
    default:
        DebugStop();
        break;
    }
    return 1;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatHCurlProjection::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    
    switch (var) {
        case 0: // E
        {
            Solout = data.sol[0];
        } break;
        case 1: // curlE
        {
            Solout[0] = data.dsol[0](2,0);
        } break;
        default:
            DebugStop();
            break;
    }
}


/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatHCurlProjection::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    
    DebugStop();
}






