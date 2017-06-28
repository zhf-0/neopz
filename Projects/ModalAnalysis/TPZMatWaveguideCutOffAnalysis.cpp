#include "TPZMatWaveguideCutOffAnalysis.h"

#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran2"));
#endif


TPZMatWaveguideCutOffAnalysis::TPZMatWaveguideCutOffAnalysis(int id, REAL f0, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &))
: TPZMatModalAnalysis::TPZMatModalAnalysis(id,f0,ur,er)
{
    
}

TPZMatWaveguideCutOffAnalysis::~TPZMatWaveguideCutOffAnalysis()
{
    
}

void TPZMatWaveguideCutOffAnalysis::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	
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
    
    /*********************CREATE HCURL FUNCTIONS*****************************/
    TPZFNMatrix< 36 , REAL > phiHCurlAxes = datavec[ hcurlmeshindex ].phi;
    TPZFNMatrix<40,REAL> curlPhiDAxes = datavec[ hcurlmeshindex ].dphix;
    
    TPZFNMatrix<40,REAL> curlPhi, phiHCurl;
    
    TPZAxesTools<REAL>::Axes2XYZ(phiHCurlAxes , phiHCurl , datavec[hcurlmeshindex].axes , false);
    
    TPZManVector<REAL,3> ax1(3),ax2(3), elNormal(3);
    for (int i=0; i<3; i++) {
        ax1[i] = datavec[ hcurlmeshindex ].axes(0,i);
        ax2[i] = datavec[ hcurlmeshindex ].axes(1,i);
    }
    Cross(ax1, ax2, elNormal);
    TPZFNMatrix<3,REAL> normalVec(1,3);
    normalVec(0,0) = elNormal[0];
    normalVec(0,1) = elNormal[1];
    normalVec(0,2) = elNormal[2];
    TPZAxesTools<REAL>::Axes2XYZ(curlPhiDAxes, curlPhi, normalVec);
    
    TPZManVector<REAL,3> x = datavec[ h1meshindex ].x;
    const STATE muR =  fUr(x);
    const STATE epsilonR = fEr(x);
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHCurlFunctions  = phiHCurl.Rows();
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = h1meshindex * nHCurlFunctions;
    const int firstHCurl = hcurlmeshindex * nH1Functions;

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE stiffAtt = 0.;
            STATE stiffBtt = 0.;
            STATE curlIdotCurlJ = 0.;
            curlIdotCurlJ += curlPhi(0 , iVec) * curlPhi(0 , jVec);
            curlIdotCurlJ += curlPhi(1 , iVec) * curlPhi(1 , jVec);
            curlIdotCurlJ += curlPhi(2 , iVec) * curlPhi(2 , jVec);
            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += phiHCurl(iVec , 0) * phiHCurl(jVec , 0);
            phiIdotPhiJ += phiHCurl(iVec , 1) * phiHCurl(jVec , 1);
            phiIdotPhiJ += phiHCurl(iVec , 2) * phiHCurl(jVec , 2);
            
            stiffAtt = 1./muR * curlIdotCurlJ;
            stiffBtt = epsilonR * phiIdotPhiJ;
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
            gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 0) * gradPhiH1(jSca , 0);
            gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 1) * gradPhiH1(jSca , 1);
            gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 2) * gradPhiH1(jSca , 2);
            
            stiffAzz =  1./muR * gradPhiScaDotGradPhiSca;
            stiffBzz = epsilonR * phiH1( iSca , 0 ) * phiH1( jSca , 0 );
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











