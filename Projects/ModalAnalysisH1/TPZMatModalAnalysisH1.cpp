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
    whichMode = modeType::NDefined;
    fKtSquared = -999;
}

TPZMatModalAnalysisH1::TPZMatModalAnalysisH1(int id) : TPZMaterial(id), fUr(urDefault),
fEr(erDefault)
{
    assembling = NDefined;
    fW = 2.*M_PI*1e+9;
    whichMode = modeType::NDefined;
    fKtSquared = -999;
}

/** @brief Default constructor */
TPZMatModalAnalysisH1::TPZMatModalAnalysisH1() : TPZMaterial(), fUr(urDefault),
fEr(erDefault)
{
    assembling = NDefined;
    fW=2.*M_PI*1e+9;
    whichMode = modeType::NDefined;
    fKtSquared = -999;
}


TPZMatModalAnalysisH1::TPZMatModalAnalysisH1(const TPZMatModalAnalysisH1 &mat) : TPZMaterial(mat), fUr(mat.fUr),
fEr(mat.fEr)
{
    assembling = NDefined;
    fW = mat.fW;
    fKtSquared = -999;
}

TPZMatModalAnalysisH1::~TPZMatModalAnalysisH1()
{
    
}


void TPZMatModalAnalysisH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix< REAL > &phi = data.phi;
    
    TPZFNMatrix<36,REAL> dphidaxes = data.dphix;
    TPZFNMatrix<3,REAL> dphi;
    TPZAxesTools<REAL>::Axes2XYZ(dphidaxes, dphi, data.axes);
    
    int nshape=phi.Rows();
    
    for(int i = 0 ; i<nshape ; i++)
    {
        for(int j=0;j<nshape;j++)
        {
            STATE stiffA = 0.;
            stiffA += dphi(0,i)*dphi(0,j);
            //stiffA -= k0Squared * phi(i,0) * phi(j,0);
            
            STATE stiffB = phi(i,0) * phi(j,0);
            
            switch (assembling) {
                case A:
                    ek(i,j) += stiffA*weight;
                    break;
                case B:
                    ek(i,j) += stiffB*weight;
                    break;
                default:
                    DebugStop();
                    break;
            }
        }
    }
}

void TPZMatModalAnalysisH1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    TPZFMatrix<REAL> &phi = data.phi;
    int nshape=phi.Rows();
    
    REAL BIG = TPZMaterial::gBigNumber;//sera posto na matriz K
    STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
    
    if (bc.Type() == 0) {//TM modes have dirichlet boundary conditions
        whichMode = modesTM;
    }else{
        whichMode = modesTE;
    }
    
    if (assembling == whichMatrix::B) {
        return;
    }
    switch ( bc.Type() ) {
        case 0: // Dirichlet
            
            for(int i = 0 ; i<nshape ; i++)
            {
                const STATE rhs = phi(i,0)*BIG*v2;
                ef(i,0) += rhs*weight;
                for(int j=0;j<nshape;j++)
                {
                    const STATE stiff = phi(i,0)*phi(j,0);
                    ek(i,j) += stiff*weight*BIG;
                }
            }
            break;
            
        case 1: // Neumann
            
            for(int i = 0 ; i<nshape ; i++)
            {
                const STATE rhs = phi(i,0)*v2;
                ef(i,0) += rhs*weight;
            }
            break;
            
        case 2: // Mista
            DebugStop();
        default:
            DebugStop();
            break;
    }
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
    if( strcmp(name.c_str(), "Hz") == 0) return 0;
    if( strcmp(name.c_str(), "Et") == 0) return 1;
    if( strcmp(name.c_str(), "Ht") == 0) return 2;
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
        case 0: //Ez ou Hz
            return 1;
            break;
        case 1: //Et
            return 2;
            break;
        case 2: //Ht
            return 2;
            break;
        default:
            DebugStop();
            break;
    }
    return 1;
}

void TPZMatModalAnalysisH1::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    TPZVec<STATE> axialField(1,0.) , gradAxialField(2,0.) , curlAxialField(2,0.);
    axialField = data.sol[0];
    TPZGradSolVec datasol = data.dsol;
    
    const STATE muR = fUr(data.x);
    const STATE epsilonR = fEr(data.x);
    const STATE k0Squared = muR * epsilonR * fW * fW / ( M_C * M_C );
    const STATE betaZ = sqrt(k0Squared - fKtSquared);
    
    gradAxialField[0] = -betaZ/(fKtSquared) * data.dsol[0](0,0);
    gradAxialField[1] = -betaZ/(fKtSquared) * data.dsol[0](1,0);
    
    
    switch (whichMode) {
        case modesTM:
            curlAxialField[0] = fW * muR /fKtSquared * ( 1.) * data.dsol[0](1,0);
            curlAxialField[1] = fW * muR /fKtSquared * (-1.) * data.dsol[0](0,0);
            break;
        case modesTE:
            
            curlAxialField[0] = fW * epsilonR /fKtSquared * ( 1.) * data.dsol[0](1,0);
            curlAxialField[1] = fW * epsilonR /fKtSquared * (-1.) * data.dsol[0](0,0);
            break;
        default:
            DebugStop();
    }
    
    switch (var) {
        case 0: //Ez ou Hz
        {
            Solout = axialField;
        }
            break;
        case 1: //Et
        {
            if(whichMode == modeType::modesTE) {
                Solout = curlAxialField;
            }
            else{
                Solout = gradAxialField;
            }
        }
            break;
        case 2: //Ht
        {
            if(whichMode == modeType::modesTE) {
                Solout = gradAxialField;
            }
            else{
                Solout = curlAxialField;
            }
        }
            break;
        default:
            DebugStop();
            break;
    }
}




