#include "TPZMatModalAnalysisHDiv.h"
#include "TPZMatHCurlProjection.h"
#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif




TPZMatModalAnalysisHDiv::TPZMatModalAnalysisHDiv(int id, REAL freq, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &) ) :
TPZVecL2(id), fUr(ur), fEr(er)
{
    assembling = NDefined;
    fW = 2.*M_PI*freq;
    whichMode = modeType::NDefined;
    fKtSquared = -999;
}

TPZMatModalAnalysisHDiv::TPZMatModalAnalysisHDiv(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault)
{
    assembling = NDefined;
    fW = 2.*M_PI*1e+9;
    whichMode = modeType::NDefined;
    fKtSquared = -999;
}

/** @brief Default constructor */
TPZMatModalAnalysisHDiv::TPZMatModalAnalysisHDiv() : TPZVecL2(), fUr(urDefault),
fEr(erDefault)
{
    assembling = NDefined;
    fW=2.*M_PI*1e+9;
    whichMode = modeType::NDefined;
    fKtSquared = -999;
}


TPZMatModalAnalysisHDiv::TPZMatModalAnalysisHDiv(const TPZMatModalAnalysisHDiv &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr)
{
    assembling = NDefined;
    fW = mat.fW;
    fKtSquared = -999;
}

TPZMatModalAnalysisHDiv::~TPZMatModalAnalysisHDiv()
{
    
}


void TPZMatModalAnalysisHDiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    /*********************CREATE L2 FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiL2 = datavec[ l2index ].phi;
    TPZFNMatrix<36,REAL> dphiL2daxes = datavec[ l2index ].dphix;
    TPZFNMatrix<3,REAL> dphiL2;
    TPZAxesTools<REAL>::Axes2XYZ(dphiL2daxes, dphiL2, datavec[ l2index ].axes);
    TPZFNMatrix<3,REAL> gradPhiL2(phiL2.Rows() , 3 , 0.);
    for ( int iFunc = 0 ; iFunc < phiL2.Rows(); iFunc++ ) {
        
        gradPhiL2 ( iFunc , 0 ) = dphiL2 ( 0 , iFunc );
        gradPhiL2 ( iFunc , 1 ) = dphiL2 ( 1 , iFunc );
    }
    
    /*********************CREATE HDIV FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiScaHDiv = datavec[ hdivindex ].phi;
    
    int phrq = datavec[ hdivindex ].fVecShapeIndex.NElements();
    //  std::cout<<"x"<<std::endl<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
    
    TPZFNMatrix< 36 , REAL > phiVecHDiv(phrq , 3 , 0.) , divPhiVecHDiv(phrq, 1, 0.);
    TPZFMatrix<REAL> &dphiQdaxes = datavec[ hdivindex ].dphix;
    TPZFNMatrix<3,REAL> dphiQ;
    TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, datavec[ hdivindex ].axes);
    for (int iq = 0 ; iq < phrq ; iq++) {
        int ivecind = datavec[ hdivindex ].fVecShapeIndex[iq].first;
        int ishapeind = datavec[ hdivindex ].fVecShapeIndex[iq].second;
        
        
        phiVecHDiv(iq , 0) = phiScaHDiv(ishapeind , 0) * datavec[ hdivindex ].fNormalVec(0,ivecind);
        phiVecHDiv(iq , 1) = phiScaHDiv(ishapeind , 0) * datavec[ hdivindex ].fNormalVec(1,ivecind);
        phiVecHDiv(iq , 2) = phiScaHDiv(ishapeind , 0) * datavec[ hdivindex ].fNormalVec(2,ivecind);
        
        divPhiVecHDiv(iq , 0) += dphiQ(0,ishapeind) * datavec[ hdivindex ].fNormalVec(0,ivecind);
        divPhiVecHDiv(iq , 0) += dphiQ(1,ishapeind) * datavec[ hdivindex ].fNormalVec(1,ivecind);
        divPhiVecHDiv(iq , 0) += dphiQ(2,ishapeind) * datavec[ hdivindex ].fNormalVec(2,ivecind);
    }

    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHDivFunctions  = phrq;
    const int nL2Functions  = phiL2.Rows();
    const int firstL2 = l2index * nHDivFunctions;
    const int firstHDiv = hdivindex * nL2Functions;
    
    for (int iVec = 0; iVec < nHDivFunctions; iVec++) {
        for (int jVec = 0; jVec < nHDivFunctions; jVec++) {
            STATE stiffAtt = 0.;
            STATE stiffBtt = 0.;
            
            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += std::conj( phiVecHDiv(iVec , 0) ) * phiVecHDiv(jVec , 0);
            phiIdotPhiJ += std::conj( phiVecHDiv(iVec , 1) ) * phiVecHDiv(jVec , 1);
            phiIdotPhiJ += std::conj( phiVecHDiv(iVec , 2) ) * phiVecHDiv(jVec , 2);
            
            stiffAtt = phiIdotPhiJ;
            stiffBtt = 0;
            
            if (this->assembling == A) {
                ek( firstHDiv + iVec , firstHDiv + jVec ) += stiffAtt * weight ;
            }
            else if (this->assembling == B){
                ek( firstHDiv + iVec , firstHDiv + jVec ) += stiffBtt * weight ;
            }
            else{
                DebugStop();
            }
            
        }
        for (int jSca = 0; jSca < nL2Functions; jSca++) {
            STATE stiffAtz = 0.;
            STATE phiDivV = 0.;
            
            phiDivV += std::conj( divPhiVecHDiv(iVec , 0) ) * phiL2(jSca , 0);
            
            stiffAtz = phiDivV;
            if (this->assembling == A) {
                ek( firstHDiv + iVec , firstL2 + jSca ) += stiffAtz * weight ;
            }
            else if (this->assembling == B){
                ek( firstHDiv + iVec , firstL2 + jSca ) += 0.;
            }
            else{
                DebugStop();
            }
            //ek( firstHCurl + iVec , firstH1 + jSca ) += stiff * weight ;
        }
    }
    for (int iSca = 0; iSca < nL2Functions; iSca++) {
        for (int jVec = 0; jVec < nHDivFunctions; jVec++) {
            
            STATE phiDivV = 0.;
            
            phiDivV += std::conj( phiL2(iSca , 0) ) * divPhiVecHDiv(jVec , 0);
            
            STATE stiffAzt = 0.;
            
            stiffAzt = phiDivV;
            if (this->assembling == A) {
                ek( firstL2 + iSca , firstHDiv +  jVec) += stiffAzt * weight ;
            }
            else if (this->assembling == B){
                ek( firstL2 + iSca , firstHDiv +  jVec ) += 0. ;
            }
            else{
                DebugStop();
            }
            
            //ek( firstH1 + iSca , firstHCurl +  jVec ) += stiff * weight ;
        }
        for (int jSca = 0; jSca < nL2Functions; jSca++) {
            STATE stiffBzz = 0.;
            
            stiffBzz =  -1. * std::conj( phiL2( iSca , 0 ) ) * phiL2( jSca , 0 );
            //ek( firstH1 + iSca , firstH1 + jSca ) += stiff * weight ;
            //            if( iSca == jSca){
            //                std::cout<<"stiffBzz "<<iSca<<" "<<jSca<<":"<<stiffBzz<<std::endl;
            //            }
            if (this->assembling == A) {
                ek( firstL2 + iSca , firstL2 + jSca) += 0. ;
            }
            else if (this->assembling == B){
                ek( firstL2 + iSca , firstL2 + jSca) += stiffBzz * weight ;
            }
            else{
                DebugStop();
            }
        }
    }
}

void TPZMatModalAnalysisHDiv::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    if (bc.Type() == 0) {//TM modes have dirichlet boundary conditions
        whichMode = modesTM;
        return;
    }else{
        whichMode = modesTE;
    }
    return;    
}

int TPZMatModalAnalysisHDiv::IntegrationRuleOrder(int elPMaxOrder) const
{
    return 2+elPMaxOrder*2;
}


int TPZMatModalAnalysisHDiv::VariableIndex(const std::string &name)
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
int TPZMatModalAnalysisHDiv::NSolutionVariables(int var)
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



void TPZMatModalAnalysisHDiv::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    
    TPZVec<STATE> axialField(1,0.);
    TPZVec<STATE> transversalField(2,0.);
    TPZVec<STATE> curlAxialField(2,0.);
    
    axialField = datavec[ l2index ].sol[0];
    transversalField = datavec[ hdivindex ].sol[0];
    
    const STATE muR = fUr(datavec[l2index].x);
    const STATE epsilonR = fEr(datavec[l2index].x);
    const STATE k0Squared = muR * epsilonR * fW * fW / ( M_C * M_C );
    const STATE betaZ = sqrt(k0Squared - fKtSquared);
    
    
    transversalField[0] = -1. * imaginary * betaZ / fKtSquared * transversalField[0];
    transversalField[1] = -1. * imaginary * betaZ / fKtSquared * transversalField[1];
    
    if (whichMode == modeType::modesTE) {
        curlAxialField[0] =  1. * transversalField[1] * fW * muR / betaZ;
        curlAxialField[1] = -1. * transversalField[0] * fW * muR / betaZ;
    }
    else{
        curlAxialField[0] =  1. * transversalField[1] * fW * epsilonR / betaZ;
        curlAxialField[1] = -1. * transversalField[0] * fW * epsilonR / betaZ;
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
                Solout = transversalField;
            }
        }
            break;
        case 2: //Ht
        {
            if(whichMode == modeType::modesTE) {
                Solout = transversalField;
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




