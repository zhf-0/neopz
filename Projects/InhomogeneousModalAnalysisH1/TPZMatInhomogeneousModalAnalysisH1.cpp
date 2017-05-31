#include "TPZMatInhomogeneousModalAnalysisH1.h"
#include "TPZMatHCurlProjection.h" //c0, u0, e0, definidos lá
#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif




TPZMatInhomogeneousModalAnalysisH1::TPZMatInhomogeneousModalAnalysisH1(int id, REAL kzOverk02, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &) ) :
TPZMaterial(id), fUr(ur), fEr(er)
{
    fIsTesting = false;
    fAssembling = NDefined;
    fKzOverK02 = kzOverk02;
}


TPZMatInhomogeneousModalAnalysisH1::~TPZMatInhomogeneousModalAnalysisH1()
{
    
}
void TPZMatInhomogeneousModalAnalysisH1::ContributeValidateFunctions(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if( fIsTesting == false ){
        DebugStop();
    }
    enum whichTest {dotSca, gradSca};
    whichTest test = dotSca;
    /*********************CREATE HZ FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiHz = datavec[ fHzIndex ].phi;
    TPZFNMatrix<36,REAL> dphiHzdaxes = datavec[ fHzIndex ].dphix;
    TPZFNMatrix<3,REAL> dphiHz;
    TPZAxesTools<REAL>::Axes2XYZ(dphiHzdaxes, dphiHz, datavec[ fHzIndex ].axes);
    TPZFNMatrix<3,REAL> gradPhiHz(phiHz.Rows() , 3 , 0.);
    
    for ( int iFunc = 0 ; iFunc < phiHz.Rows(); iFunc++ ) {
        
        gradPhiHz ( iFunc , 0 ) = dphiHz ( 0 , iFunc );
        gradPhiHz ( iFunc , 1 ) = dphiHz ( 1 , iFunc );
        
    }
    /*********************CREATE EZ FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiEz = datavec[ fEzIndex ].phi;
    TPZFNMatrix<36,REAL> dphiEzdaxes = datavec[ fEzIndex ].dphix;
    TPZFNMatrix<3,REAL> dphiEz;
    TPZAxesTools<REAL>::Axes2XYZ(dphiEzdaxes, dphiEz, datavec[ fEzIndex ].axes);
    TPZFNMatrix<3,REAL> gradPhiEz(phiEz.Rows() , 3 , 0.);
    
    for ( int iFunc = 0 ; iFunc < phiEz.Rows(); iFunc++ ) {
        
        gradPhiEz ( iFunc , 0 ) = dphiEz ( 0 , iFunc );
        gradPhiEz ( iFunc , 1 ) = dphiEz ( 1 , iFunc );
        
    }
    
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHzFunctions  = phiHz.Rows();
    const int nEzFunctions  = phiEz.Rows();
    const int firstHz = fHzIndex * nEzFunctions;
    const int firstEz = fEzIndex * nHzFunctions;
    
    for (int iHz = 0; iHz < nHzFunctions; iHz++) {
        for (int jHz = 0; jHz < nHzFunctions; jHz++) {
            STATE gradPhiScaDotGradPhiSca = 0.;
            gradPhiScaDotGradPhiSca += gradPhiHz(iHz , 0) * gradPhiHz(jHz , 0);
            gradPhiScaDotGradPhiSca += gradPhiHz(iHz , 1) * gradPhiHz(jHz , 1);
            gradPhiScaDotGradPhiSca += gradPhiHz(iHz , 2) * gradPhiHz(jHz , 2);
            
            if (test == dotSca) {
                ek( firstHz + iHz , firstHz + jHz) += phiHz( iHz , 0 ) * phiHz( jHz , 0 ) * weight;
            }
            else if (test == gradSca){
                ek( firstHz + iHz , firstHz + jHz) += gradPhiScaDotGradPhiSca * weight ;
            }
        }
        for (int jEz = 0; jEz < nHzFunctions; jEz++) {
            STATE gradPhiScaDotGradPhiSca = 0.;
            gradPhiScaDotGradPhiSca += gradPhiHz( iHz , 0) * gradPhiEz( jEz , 0);
            gradPhiScaDotGradPhiSca += gradPhiHz( iHz , 1) * gradPhiEz( jEz , 1);
            gradPhiScaDotGradPhiSca += gradPhiHz( iHz , 2) * gradPhiEz( jEz , 2);
            
            if (test == dotSca) {
                ek( firstHz + iHz , firstEz + jEz) += phiHz( iHz , 0 ) * phiEz( jEz , 0 ) * weight;
            }
            else if (test == gradSca){
                ek( firstHz + iHz , firstEz + jEz) += gradPhiScaDotGradPhiSca * weight ;
            }
        }
    }
    for (int iEz = 0; iEz < nEzFunctions; iEz++) {
        for (int jHz = 0; jHz < nHzFunctions; jHz++) {
            STATE gradPhiScaDotGradPhiSca = 0.;
            gradPhiScaDotGradPhiSca += gradPhiEz( iEz , 0) * gradPhiHz( jHz , 0);
            gradPhiScaDotGradPhiSca += gradPhiEz( iEz , 1) * gradPhiHz( jHz , 1);
            gradPhiScaDotGradPhiSca += gradPhiEz( iEz , 2) * gradPhiHz( jHz , 2);
            
            if (test == dotSca) {
                ek( firstEz + iEz , firstHz + jHz) += phiEz( iEz , 0 ) * phiHz( jHz , 0 ) * weight;
            }
            else if (test == gradSca){
                ek( firstEz + iEz , firstHz + jHz) += gradPhiScaDotGradPhiSca * weight ;
            }
        }
        for (int jEz = 0; jEz < nEzFunctions; jEz++) {
            STATE gradPhiScaDotGradPhiSca = 0.;
            gradPhiScaDotGradPhiSca += gradPhiEz( iEz , 0) * gradPhiEz( jEz , 0);
            gradPhiScaDotGradPhiSca += gradPhiEz( iEz , 1) * gradPhiEz( jEz , 1);
            gradPhiScaDotGradPhiSca += gradPhiEz( iEz , 2) * gradPhiEz( jEz , 2);
            
            if (test == dotSca) {
                ek( firstEz + iEz , firstEz + jEz) += phiEz( iEz , 0 ) * phiEz( jEz , 0 ) * weight;
            }
            else if (test == gradSca){
                ek( firstEz + iEz , firstEz + jEz) += gradPhiScaDotGradPhiSca * weight ;
            }
        }
    }
}


void TPZMatInhomogeneousModalAnalysisH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatInhomogeneousModalAnalysisH1::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    fIsTesting = false;
    if( fIsTesting == true ){
        ContributeValidateFunctions(datavec, weight, ek, ef);
        return;
    }
    
    /*********************CREATE HZ FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiHz = datavec[ fHzIndex ].phi;
    TPZFNMatrix<36,REAL> dphiHzdaxes = datavec[ fHzIndex ].dphix;
    TPZFNMatrix<3,REAL> dphiHz;
    TPZAxesTools<REAL>::Axes2XYZ(dphiHzdaxes, dphiHz, datavec[ fHzIndex ].axes);
    TPZFNMatrix<2,REAL> gradPhiHz(phiHz.Rows() , 2 , 0.);
    
    for ( int iFunc = 0 ; iFunc < phiHz.Rows(); iFunc++ ) {
        
        gradPhiHz ( iFunc , 0 ) = dphiHz ( 0 , iFunc );
        gradPhiHz ( iFunc , 1 ) = dphiHz ( 1 , iFunc );
        
    }
    /*********************CREATE EZ FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiEz = datavec[ fEzIndex ].phi;
    TPZFNMatrix<36,REAL> dphiEzdaxes = datavec[ fEzIndex ].dphix;
    TPZFNMatrix<3,REAL> dphiEz;
    TPZAxesTools<REAL>::Axes2XYZ(dphiEzdaxes, dphiEz, datavec[ fEzIndex ].axes);
    TPZFNMatrix<2,REAL> gradPhiEz(phiEz.Rows() , 2 , 0.);
    
    for ( int iFunc = 0 ; iFunc < phiEz.Rows(); iFunc++ ) {
        
        gradPhiEz ( iFunc , 0 ) = dphiEz ( 0 , iFunc );
        gradPhiEz ( iFunc , 1 ) = dphiEz ( 1 , iFunc );
        
    }
    
    //*****************MEDIA INFORMATION****************//
    const STATE muR =  fUr(datavec[0].x);
    const STATE epsilonR = fEr(datavec[0].x);
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHzFunctions  = phiHz.Rows();
    const int nEzFunctions  = phiEz.Rows();
    const int firstHz = fHzIndex * nEzFunctions;
    const int firstEz = fEzIndex * nHzFunctions;
    
    for (int iEz = 0; iEz < nEzFunctions; iEz++) {
        for (int jEz = 0; jEz < nEzFunctions; jEz++) {
            if (this->fAssembling == whichMatrix::A) {
                STATE gradPhiScaDotGradPhiSca = 0.;
                gradPhiScaDotGradPhiSca += gradPhiEz( iEz , 0) * gradPhiEz( jEz , 0);
                gradPhiScaDotGradPhiSca += gradPhiEz( iEz , 1) * gradPhiEz( jEz , 1);
                const STATE cte = (epsilonR*M_EZERO)/(muR*epsilonR-fKzOverK02)/sqrt(M_UZERO*M_EZERO);
                ek( firstEz + iEz , firstEz + jEz) += cte*gradPhiScaDotGradPhiSca * weight ;
            }
            else {
                const STATE cte = epsilonR*M_EZERO/sqrt(M_UZERO*M_EZERO);
                ek( firstEz + iEz , firstEz + jEz) += cte * phiEz( iEz , 0 ) * phiEz( jEz , 0 ) * weight ;
            }
        }
        for (int jHz = 0; jHz < nHzFunctions; jHz++) {
            if (this->fAssembling == whichMatrix::A) {
                STATE gradPhiScaVecGradPhiSca = 0.;
                gradPhiScaVecGradPhiSca += gradPhiEz( iEz , 0) * gradPhiHz( jHz , 1);
                gradPhiScaVecGradPhiSca -= gradPhiEz( iEz , 1) * gradPhiHz( jHz , 0);
                
                const STATE cte = (sqrt(fKzOverK02) / M_C)/(muR*epsilonR-fKzOverK02)/sqrt(M_UZERO*M_EZERO);
                ek( firstEz + iEz , firstHz + jHz) += cte * gradPhiScaVecGradPhiSca * weight ;
            }
        }
    }
    
    for (int iHz = 0; iHz < nHzFunctions; iHz++) {
        for (int jHz = 0; jHz < nHzFunctions; jHz++) {
            
            if (this->fAssembling == whichMatrix::A) {
                STATE gradPhiScaDotGradPhiSca = 0.;
                gradPhiScaDotGradPhiSca += gradPhiHz(iHz , 0) * gradPhiHz(jHz , 0);
                gradPhiScaDotGradPhiSca += gradPhiHz(iHz , 1) * gradPhiHz(jHz , 1);
                
                const STATE cte = (muR*M_UZERO)/(muR*epsilonR-fKzOverK02) / (sqrt(M_UZERO*M_EZERO));
                ek( firstHz + iHz , firstHz + jHz) += cte*gradPhiScaDotGradPhiSca * weight ;
            }
            else {
                const STATE cte = muR*M_UZERO / (sqrt(M_UZERO*M_EZERO));
                ek( firstHz + iHz , firstHz + jHz) += cte * phiHz( iHz , 0 ) * phiHz( jHz , 0 ) * weight ;
            }
        }
        for (int jEz = 0; jEz < nHzFunctions; jEz++) {
            if (this->fAssembling == whichMatrix::A) {
                STATE gradPhiScaVecGradPhiSca = 0.;
                gradPhiScaVecGradPhiSca += gradPhiHz( iHz , 0) * gradPhiEz( jEz , 1);
                gradPhiScaVecGradPhiSca -= gradPhiHz( iHz , 1) * gradPhiEz( jEz , 0);
                
                const STATE cte = (sqrt(fKzOverK02) / M_C)/(muR*epsilonR-fKzOverK02)/(sqrt(M_UZERO*M_EZERO));
                ek( firstHz + iHz , firstEz + jEz) += cte * gradPhiScaVecGradPhiSca * weight ;
            }
        }
    }
}

void TPZMatInhomogeneousModalAnalysisH1::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatInhomogeneousModalAnalysisH1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}


void TPZMatInhomogeneousModalAnalysisH1::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    return;
}

void TPZMatInhomogeneousModalAnalysisH1::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatInhomogeneousModalAnalysisH1::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

int TPZMatInhomogeneousModalAnalysisH1::IntegrationRuleOrder(int elPMaxOrder) const
{
    return 2+elPMaxOrder*2;
}


int TPZMatInhomogeneousModalAnalysisH1::VariableIndex(const std::string &name)
{
    if( strcmp(name.c_str(), "Ez") == 0) return 0;
    if( strcmp(name.c_str(), "Hz") == 0) return 1;
    if( strcmp(name.c_str(), "Et") == 0) return 2;
    if( strcmp(name.c_str(), "Ht") == 0) return 3;
    DebugStop();
    return 1;
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatInhomogeneousModalAnalysisH1::NSolutionVariables(int var)
{
    
    switch (var) {
        case 0: //Ez
            return 1;
            break;
        case 1: //Hz
            return 1;
            break;
        case 2: //Et
            return 2;
            break;
        case 3: //Ht
            return 2;
            break;
        default:
            DebugStop();
            break;
    }
    return 1;
    
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatInhomogeneousModalAnalysisH1::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    
    DebugStop();
}


/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatInhomogeneousModalAnalysisH1::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    
    TPZVec<STATE> hz(1,0.);
    TPZVec<STATE> ez(1,0.);
    
    ez = datavec[ fEzIndex ].sol[0];
    hz = datavec[ fHzIndex ].sol[0];
    
    const STATE muR = fUr(datavec[0].x);
    const STATE epsilonR = fEr(datavec[0].x);
//    const STATE k0Squared = muR * epsilonR * fW * fW / ( M_C * M_C );
//    const STATE betaZ = sqrt(k0Squared - fKtSquared);
    
    switch (var) {
        case 0: //Ez
        {
            Solout = ez;
        }
            break;
        case 1://Hz
        {
            Solout = hz;
        }
            break;
        default:
            DebugStop();
            break;
    }
}
