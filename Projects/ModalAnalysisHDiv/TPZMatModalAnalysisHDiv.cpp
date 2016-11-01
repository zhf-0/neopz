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
    TPZFNMatrix<12,REAL> phiScaHCurl = datavec[ hdivindex ].phi;
    TPZManVector<REAL,3> x = datavec[ l2index ].x;
    TPZManVector<REAL,3> xParametric = datavec[ l2index ].xParametric;
    
    int phrq = datavec[ hdivindex ].fVecShapeIndex.NElements();
    //  std::cout<<"x"<<std::endl<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
    
    TPZFNMatrix< 36 , REAL > phiVecHDiv(phrq , 3 , 0.) , divPhiVecHDiv(phrq, 1, 0.);
    TPZFMatrix<REAL> &dphiQdaxes = datavec[ hdivindex ].dphix;
    TPZFNMatrix<3,REAL> dphiQ;
    TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, datavec[ hdivindex ].axes);
    for (int iq = 0 ; iq < phrq ; iq++) {
        int ivecind = datavec[ hdivindex ].fVecShapeIndex[iq].first;
        int ishapeind = datavec[ hdivindex ].fVecShapeIndex[iq].second;
        
        
        phiVecHDiv(iq , 0) = phiScaHCurl(ishapeind , 0) * datavec[ hdivindex ].fNormalVec(0,ivecind);
        phiVecHDiv(iq , 1) = phiScaHCurl(ishapeind , 0) * datavec[ hdivindex ].fNormalVec(1,ivecind);
        phiVecHDiv(iq , 2) = phiScaHCurl(ishapeind , 0) * datavec[ hdivindex ].fNormalVec(2,ivecind);
        
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
    TPZFMatrix<REAL> &phiHDiv = datavec[hdivindex].phi;
    TPZFMatrix<REAL> &phiL2 = datavec[l2index].phi;
    const int nHDivFunctions  = phiHDiv.Rows();
    const int nL2Functions  = phiL2.Rows();
    const int firstL2 = l2index * nHDivFunctions;
    const int firstHDiv = hdivindex * nL2Functions;
    
    if (bc.Type() == 0) {//TM modes have dirichlet boundary conditions
        whichMode = modesTM;
        return;
    }else{
        whichMode = modesTE;
    }
    return;
    
    
    if (assembling == whichMatrix::B) {
        return;
    }
    
    {
        
        int nshape=phiL2.Rows();
        REAL BIG = TPZMaterial::gBigNumber;
        
        //const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
        const STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
        
        switch ( bc.Type() )
        {
            case 0:
                for(int i = 0 ; i<nshape ; i++)
                {
                    const STATE rhs = phiL2(i,0) * BIG  * v2;
                    ef(firstL2+i,0) += rhs*weight;
                    for(int j=0;j<nshape;j++)
                    {
                        const STATE stiff = phiL2(i,0) * phiL2(j,0) * BIG ;
                        ek(firstL2+i,firstL2+j) += stiff*weight;
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
    {
        
        int nshape=phiHDiv.Rows();
        REAL BIG = TPZMaterial::gBigNumber;
        
        //const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
        const STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
        
        switch ( bc.Type() )
        {
            case 0:
                DebugStop();
                break;
            case 1:
                for(int i = 0 ; i<nshape ; i++)
                {
                    const STATE rhs = phiHDiv(i,0) * BIG  * v2;
                    ef(firstHDiv+i,0) += rhs*weight;
                    for(int j=0;j<nshape;j++)
                    {
                        const STATE stiff = phiHDiv(i,0) * phiHDiv(j,0) * BIG ;
                        ek(firstHDiv+i,firstHDiv+j) += stiff*weight;
                    }
                }
                break;
            case 2:
                DebugStop();
                break;
        }
    }
    
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

    TPZVec<STATE> axialField(1,0.) , gradAxialField(2,0.) , curlAxialField(2,0.);
    axialField = datavec[l2index].sol[0];
    TPZManVector<STATE> datasoldiv = datavec[hdivindex].sol[0];
    const STATE muR =  fUr(datavec[l2index].x);
    const STATE epsilonR = fEr(datavec[l2index].x);
    const STATE k0Squared = muR * epsilonR * fW * fW / ( M_C * M_C );
    const STATE betaZ = sqrt(k0Squared - fKtSquared);
    
    gradAxialField[0] = -betaZ/(fKtSquared) * datasoldiv[0];
    gradAxialField[1] = -betaZ/(fKtSquared) * datasoldiv[1];
    
    
    switch (whichMode) {
        case modesTM:
            curlAxialField[0] = fW * muR /fKtSquared * ( 1.) * datasoldiv[1];
            curlAxialField[1] = fW * muR /fKtSquared * (-1.) * datasoldiv[0];
            break;
        case modesTE:
            curlAxialField[0] = fW * epsilonR /fKtSquared * ( 1.) * datasoldiv[1];
            curlAxialField[1] = fW * epsilonR /fKtSquared * (-1.) * datasoldiv[0];
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




