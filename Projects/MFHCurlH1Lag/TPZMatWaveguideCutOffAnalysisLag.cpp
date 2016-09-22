#include "TPZMatWaveguideCutOffAnalysisLag.h"

#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran2"));
#endif


TPZMatWaveguideCutOffAnalysisLag::TPZMatWaveguideCutOffAnalysisLag(int id, REAL f0, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &))
: TPZMatMFHCurlH1Lag::TPZMatMFHCurlH1Lag(id,f0,ur,er)
{
    
}

TPZMatWaveguideCutOffAnalysisLag::~TPZMatWaveguideCutOffAnalysisLag()
{
    
}

void TPZMatWaveguideCutOffAnalysisLag::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    isTesting = false;
    if( isTesting == true ){
        ContributeValidateFunctions(datavec, weight, ek, ef);
        return;
    }
    /*********************CREATE H1 FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiH1 = datavec[ h1meshindex ].phi;
    TPZFNMatrix<36,REAL> dphiH1daxes = datavec[ h1meshindex ].dphix;
    TPZFNMatrix<3,REAL> dphiH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiH1daxes, dphiH1, datavec[ hcurlmeshindex ].axes);
    TPZFNMatrix<3,REAL> gradPhiH1(phiH1.Rows() , 3 , 0.);
    for ( int iFunc = 0 ; iFunc < phiH1.Rows(); iFunc++ ) {
        
        gradPhiH1 ( iFunc , 0 ) = dphiH1 ( 0 , iFunc );
        gradPhiH1 ( iFunc , 1 ) = dphiH1 ( 1 , iFunc );
    }
    
    /*********************CREATE HDIV FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiScaHCurl = datavec[ hcurlmeshindex ].phi;
    TPZManVector<REAL,3> x = datavec[ h1meshindex ].x;
    TPZManVector<REAL,3> xParametric = datavec[ h1meshindex ].xParametric;
    
    int phrq = datavec[ hcurlmeshindex ].fVecShapeIndex.NElements();
    //  std::cout<<"x"<<std::endl<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
    
    TPZFNMatrix< 36 , REAL > phiVecHDiv(phrq , 3 , 0.);
    for (int iq = 0 ; iq < phrq ; iq++) {
        int ivecind = datavec[ hcurlmeshindex ].fVecShapeIndex[iq].first;
        int ishapeind = datavec[ hcurlmeshindex ].fVecShapeIndex[iq].second;
        
        phiVecHDiv(iq , 0) = phiScaHCurl(ishapeind , 0) * datavec[ hcurlmeshindex ].fNormalVec(0,ivecind);
        phiVecHDiv(iq , 1) = phiScaHCurl(ishapeind , 0) * datavec[ hcurlmeshindex ].fNormalVec(1,ivecind);
        phiVecHDiv(iq , 2) = phiScaHCurl(ishapeind , 0) * datavec[ hcurlmeshindex ].fNormalVec(2,ivecind);
    }
    
    /*********************CALCULATE NORMAL VECTOR****************************/
    TPZManVector<REAL,3> ax1(3),ax2(3), elNormal(3);
    for (int i=0; i<3; i++) {
        ax1[i] = datavec[ hcurlmeshindex ].axes(0,i);
        ax2[i] = datavec[ hcurlmeshindex ].axes(1,i);
    }
    Cross(ax1, ax2, elNormal);
    
    /*********************CREATE HCURL FUNCTIONS*****************************/
    TPZFNMatrix< 36 , REAL > phiVecHCurl(phrq , 3 , 0.);
    RotateForHCurl(elNormal , phiVecHDiv , phiVecHCurl);
    /*********************COMPUTE CURL****************************/
    TPZFMatrix<REAL> &dphiQdaxes = datavec[ hcurlmeshindex ].dphix;
    TPZFNMatrix<3,REAL> dphiQ;
    TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, datavec[ hcurlmeshindex ].axes);
    TPZFNMatrix<3,REAL> gradPhiForHCurl(phrq , 3 , 0.);
    TPZFNMatrix<3,REAL> ivecHCurl(phrq , 3 , 0.);
    TPZManVector<REAL,3> iVecHDiv(3,0.), ivecForCurl(3,0.);
    for (int iPhi = 0; iPhi < phrq; iPhi++) {
        int ivecind = datavec[ hcurlmeshindex ].fVecShapeIndex[iPhi].first;
        int ishapeind = datavec[ hcurlmeshindex ].fVecShapeIndex[iPhi].second;
        iVecHDiv[0] = datavec[ hcurlmeshindex ].fNormalVec(0,ivecind);
        iVecHDiv[1] = datavec[ hcurlmeshindex ].fNormalVec(1,ivecind);
        iVecHDiv[2] = datavec[ hcurlmeshindex ].fNormalVec(2,ivecind);
        Cross(elNormal, iVecHDiv, ivecForCurl);
        for (int i = 0; i<dphiQ.Rows(); i++) {
            gradPhiForHCurl(iPhi,i) = dphiQ(i,ishapeind);
            ivecHCurl(iPhi,i) = ivecForCurl[i];
        }
    }
    TPZFNMatrix<40,REAL> curlPhi;
    ComputeCurl(gradPhiForHCurl, ivecHCurl, curlPhi);
    
    const STATE muR =  fUr(x);
    const STATE epsilonR = fEr(x);
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHCurlFunctions  = phrq;
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = h1meshindex * nHCurlFunctions;
    const int firstHCurl = hcurlmeshindex * nH1Functions;

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE stiffAtt = 0.;
            STATE stiffBtt = 0.;
            STATE curlIdotCurlJ = 0.;
            curlIdotCurlJ += std::conj( curlPhi(iVec , 0) ) * curlPhi(jVec , 0);
            curlIdotCurlJ += std::conj( curlPhi(iVec , 1) ) * curlPhi(jVec , 1);
            curlIdotCurlJ += std::conj( curlPhi(iVec , 2) ) * curlPhi(jVec , 2);
            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 0) ) * phiVecHCurl(jVec , 0);
            phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 1) ) * phiVecHCurl(jVec , 1);
            phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 2) ) * phiVecHCurl(jVec , 2);
            
            stiffAtt = 1./muR * curlIdotCurlJ;
            stiffBtt = epsilonR * phiIdotPhiJ;;
//            if( iVec == jVec ){
//                std::cout<<"stiffBtt "<<iVec<<" "<<jVec<<":"<<stiffBtt<<std::endl;
//            }
            //ek( firstHCurl + iVec , firstHCurl + jVec ) += curlIdotCurlJ * weight ;
            if (this->assembling == A) {
              ek( firstHCurl + iVec , firstHCurl + jVec ) += stiffAtt * weight ;
            }
            else if (this->assembling == B){
              ek( firstHCurl + iVec , firstHCurl + jVec ) += stiffBtt * weight ;
            }
            else{
                DebugStop();
            }
            
        }
    }
    for (int iSca = 0; iSca < nH1Functions; iSca++) {
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE gradPhiScaDotGradPhiSca = 0.;
            STATE stiffAzz = 0.;
            STATE stiffBzz = 0.;
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 0) ) * gradPhiH1(jSca , 0);
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 1) ) * gradPhiH1(jSca , 1);
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 2) ) * gradPhiH1(jSca , 2);
            
            stiffAzz =  1./muR * gradPhiScaDotGradPhiSca;
            stiffBzz = epsilonR * std::conj( phiH1( iSca , 0 ) ) * phiH1( jSca , 0 );
            //ek( firstH1 + iSca , firstH1 + jSca ) += stiff * weight ;
//            if( iSca == jSca){
//                std::cout<<"stiffBzz "<<iSca<<" "<<jSca<<":"<<stiffBzz<<std::endl;
//            }
            if (this->assembling == A) {
                ek( firstH1 + iSca , firstH1 + jSca) += stiffAzz * weight ;
            }
            else if (this->assembling == B){
                ek( firstH1 + iSca , firstH1 + jSca) += stiffBzz * weight ;
            }
            else{
                DebugStop();
            }
        }
    }
}

void TPZMatWaveguideCutOffAnalysisLag::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}











