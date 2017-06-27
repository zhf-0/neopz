#include "TPZMatMFHDivRotH1.h"

#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif




TPZMatMFHDivRotH1::TPZMatMFHDivRotH1(int id, REAL freq, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &) ) :
TPZMatModalAnalysis(id,freq,ur,er)
{
}

TPZMatMFHDivRotH1::TPZMatMFHDivRotH1(int id) : TPZMatModalAnalysis(id)
{
}

/** @brief Default constructor */
TPZMatMFHDivRotH1::TPZMatMFHDivRotH1() : TPZMatModalAnalysis()
{
}


TPZMatMFHDivRotH1::TPZMatMFHDivRotH1(const TPZMatMFHDivRotH1 &mat) : TPZMatModalAnalysis(mat){
}

TPZMatMFHDivRotH1::~TPZMatMFHDivRotH1()
{
    
}
void TPZMatMFHDivRotH1::RotateForHCurl(TPZVec<REAL> normal , TPZFMatrix<REAL> vHdiv , TPZFMatrix<REAL> &vHcurl ){
    int nFunctions = vHdiv.Rows();
    vHcurl.Resize( vHdiv.Rows(), vHdiv.Cols());
    vHcurl.Zero();
    
    for (int i = 0 ; i < nFunctions; i++) {
        
        vHcurl(i,0) = normal[1]*vHdiv(i,2) - vHdiv(i,1)*normal[2];
        vHcurl(i,1) = normal[2]*vHdiv(i,0) - vHdiv(i,2)*normal[0];
        vHcurl(i,2) = normal[0]*vHdiv(i,1) - vHdiv(i,0)*normal[1];
    }
}
void TPZMatMFHDivRotH1::ComputeCurl(TPZFMatrix<REAL> gradScalarPhi , TPZFMatrix<REAL> ivecHCurl , TPZFMatrix<REAL> &curlPhi ){
    int nFunctions = gradScalarPhi.Rows();
    curlPhi.Resize( nFunctions , 3);
    curlPhi.Zero();
    TPZManVector<REAL,3> result(3,0.);
    
    for (int i = 0 ; i < nFunctions; i++) {
        curlPhi(i,0) = gradScalarPhi(i,1)*ivecHCurl(i,2) - ivecHCurl(i,1)*gradScalarPhi(i,2);
        curlPhi(i,1) = gradScalarPhi(i,2)*ivecHCurl(i,0) - ivecHCurl(i,2)*gradScalarPhi(i,0);
        curlPhi(i,2) = gradScalarPhi(i,0)*ivecHCurl(i,1) - ivecHCurl(i,0)*gradScalarPhi(i,1);
    }
}

void TPZMatMFHDivRotH1::ContributeValidateFunctions(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if( isTesting == false ){
        DebugStop();
    }
    enum whichTest {curl = 0 , dotVec, dotSca, mixed, gradSca};
    whichTest test = mixed;
    /*********************CREATE H1 FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiH1 = datavec[ h1meshindex ].phi;
    TPZFNMatrix<36,REAL> dphiH1daxes = datavec[ h1meshindex ].dphix;
    TPZFNMatrix<3,REAL> dphiH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiH1daxes, dphiH1, datavec[ h1meshindex ].axes);
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
    
    
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHCurlFunctions  = phrq;
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = h1meshindex * nHCurlFunctions;
    const int firstHCurl = hcurlmeshindex * nH1Functions;
    
    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE curlIdotCurlJ = 0.;
            curlIdotCurlJ += std::conj( curlPhi(iVec , 0) ) * curlPhi(jVec , 0);
            curlIdotCurlJ += std::conj( curlPhi(iVec , 1) ) * curlPhi(jVec , 1);
            curlIdotCurlJ += std::conj( curlPhi(iVec , 2) ) * curlPhi(jVec , 2);
            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 0) ) * phiVecHCurl(jVec , 0);
            phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 1) ) * phiVecHCurl(jVec , 1);
            phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 2) ) * phiVecHCurl(jVec , 2);
            if (test == curl) {
                ek( firstHCurl + iVec , firstHCurl + jVec ) += curlIdotCurlJ * weight ;
            }
            else if( test == dotVec){
                ek( firstHCurl + iVec , firstHCurl + jVec ) += phiIdotPhiJ * weight ;
            }
        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE phiVecDotGradPhiSca = 0.;
            
            phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 0) ) * gradPhiH1(jSca , 0);
            phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 1) ) * gradPhiH1(jSca , 1);
            phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 2) ) * gradPhiH1(jSca , 2);
            if ( test == mixed) {
                ek( firstHCurl + iVec , firstH1 + jSca ) += phiVecDotGradPhiSca * weight ;
            }
        }
    }
    for (int iSca = 0; iSca < nH1Functions; iSca++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE phiVecDotGradPhiSca = 0.;
            phiVecDotGradPhiSca += phiVecHCurl(jVec , 0) * std::conj( gradPhiH1(iSca , 0) );
            phiVecDotGradPhiSca += phiVecHCurl(jVec , 1) * std::conj( gradPhiH1(iSca , 1) );
            phiVecDotGradPhiSca += phiVecHCurl(jVec , 2) * std::conj( gradPhiH1(iSca , 2) );
            
            if ( test == mixed) {
                ek( firstH1 + iSca , firstHCurl +  jVec ) += phiVecDotGradPhiSca * weight ;
            }
        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE gradPhiScaDotGradPhiSca = 0.;
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 0) ) * gradPhiH1(jSca , 0);
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 1) ) * gradPhiH1(jSca , 1);
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 2) ) * gradPhiH1(jSca , 2);
            
            if (test == dotSca) {
                ek( firstH1 + iSca , firstH1 + jSca) += std::conj( phiH1( iSca , 0 ) ) * phiH1( jSca , 0 ) * weight;
            }
            else if (test == gradSca){
                ek( firstH1 + iSca , firstH1 + jSca) += gradPhiScaDotGradPhiSca * weight ;
            }
        }
    }
}

void TPZMatMFHDivRotH1::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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
    TPZAxesTools<REAL>::Axes2XYZ(dphiH1daxes, dphiH1, datavec[ h1meshindex ].axes);
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
    REAL k0 = fW*sqrt(M_EZERO*M_UZERO);
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
            stiffAtt -= k0 * k0 * epsilonR * phiIdotPhiJ;
            stiffBtt = 1./muR * phiIdotPhiJ;
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
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            //      std::cout<<"function: "<<jSca<<std::endl
            //      <<std::setw(10)<<phiSca( jSca , 0 )<<std::endl;
            //      std::cout<<"grad: "<<std::endl
            //      <<std::setw(10)<<gradPhiSca( jSca , 0 )<<" "
            //      <<std::setw(10)<<gradPhiSca( jSca , 1 )<<" "
            //      <<std::setw(10)<<gradPhiSca( jSca , 2 )<<std::endl;
            STATE stiffBzt = 0.;
            STATE phiVecDotGradPhiSca = 0.;
            
            phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 0) ) * gradPhiH1(jSca , 0);
            phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 1) ) * gradPhiH1(jSca , 1);
            phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 2) ) * gradPhiH1(jSca , 2);
            
            stiffBzt = 1./muR * phiVecDotGradPhiSca;
            if (this->assembling == A) {
                ek( firstHCurl + iVec , firstH1 + jSca ) += 0. ;
            }
            else if (this->assembling == B){
                ek( firstHCurl + iVec , firstH1 + jSca ) += stiffBzt * weight ;
            }
            else{
                DebugStop();
            }
            //ek( firstHCurl + iVec , firstH1 + jSca ) += stiff * weight ;
        }
    }
    for (int iSca = 0; iSca < nH1Functions; iSca++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE phiVecDotGradPhiSca = 0.;
            STATE stiffBtz = 0.;
            phiVecDotGradPhiSca += phiVecHCurl(jVec , 0) * std::conj( gradPhiH1(iSca , 0) );
            phiVecDotGradPhiSca += phiVecHCurl(jVec , 1) * std::conj( gradPhiH1(iSca , 1) );
            phiVecDotGradPhiSca += phiVecHCurl(jVec , 2) * std::conj( gradPhiH1(iSca , 2) );
            stiffBtz = 1./muR * phiVecDotGradPhiSca;
            if (this->assembling == A) {
                ek( firstH1 + iSca , firstHCurl +  jVec) += 0. ;
            }
            else if (this->assembling == B){
                ek( firstH1 + iSca , firstHCurl +  jVec ) += stiffBtz * weight ;
            }
            else{
                DebugStop();
            }
            
            //ek( firstH1 + iSca , firstHCurl +  jVec ) += stiff * weight ;
        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE gradPhiScaDotGradPhiSca = 0.;
            STATE stiffBzz = 0.;
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 0) ) * gradPhiH1(jSca , 0);
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 1) ) * gradPhiH1(jSca , 1);
            gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 2) ) * gradPhiH1(jSca , 2);
            
            stiffBzz =  1./muR * gradPhiScaDotGradPhiSca;
            stiffBzz -=  k0 * k0 * epsilonR * std::conj( phiH1( iSca , 0 ) ) * phiH1( jSca , 0 );
            //ek( firstH1 + iSca , firstH1 + jSca ) += stiff * weight ;
//            if( iSca == jSca){
//                std::cout<<"stiffBzz "<<iSca<<" "<<jSca<<":"<<stiffBzz<<std::endl;
//            }
            if (this->assembling == A) {
                ek( firstH1 + iSca , firstH1 + jSca) += 0. ;
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

void TPZMatMFHDivRotH1::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    return;
    if( isTesting ) return;
    //if( this->assembling == B ) return;
    
    
    TPZFMatrix<REAL> &phiHCurl = datavec[hcurlmeshindex].phi;
    TPZFMatrix<REAL> &phiH1 = datavec[h1meshindex].phi;
    const int nHCurlFunctions  = phiHCurl.Rows();
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = h1meshindex * nHCurlFunctions;
    const int firstHCurl = hcurlmeshindex * nH1Functions;
    {
        
        int nshape=phiH1.Rows();
        REAL BIG = TPZMaterial::gBigNumber;
        
        //const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
        const STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
        
        switch ( bc.Type() )
        {
            case 0:
                if(this->assembling == B || this->assembling == A){
                    for(int i = 0 ; i<nshape ; i++)
                    {
                        const STATE rhs = phiH1(i,0) * BIG  * v2;
                        ef(firstH1+i,0) += rhs*weight;
                        for(int j=0;j<nshape;j++)
                        {
                            const STATE stiff = phiH1(i,0) * phiH1(j,0) * BIG ;
                            ek(firstH1+i,firstH1+j) += stiff*weight;
                        }
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
        
        int nshape=phiHCurl.Rows();
        REAL BIG = TPZMaterial::gBigNumber;
        
        //const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
        const STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
        
        switch ( bc.Type() )
        {
            case 0:
                for(int i = 0 ; i<nshape ; i++)
                {
                    const STATE rhs = phiHCurl(i,0) * BIG  * v2;
                    ef(firstHCurl+i,0) += rhs*weight;
                    for(int j=0;j<nshape;j++)
                    {
                        const STATE stiff = phiHCurl(i,0) * phiHCurl(j,0) * BIG ;
                        ek(firstHCurl+i,firstHCurl+j) += stiff*weight;
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

}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatMFHDivRotH1::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    
    TPZVec<STATE> et(3,0.);
    TPZVec<STATE> ez(1,0.);
    TPZManVector<STATE,3> ax1(3),ax2(3), normal(3);
    for (int i=0; i<3; i++) {
        ax1[i] = datavec[hcurlmeshindex].axes(0,i);
        ax2[i] = datavec[hcurlmeshindex].axes(1,i);
    }
    //ROTATE FOR HCURL
    Cross(ax1, ax2, normal);
    
    Cross(normal,datavec[ hcurlmeshindex ].sol[0], et);
    ez = datavec[ h1meshindex ].sol[0];
    switch (var) {
        case 0: //Et
        {
            Solout = et;
        }
            break;
        case 1://Ez
        {
            Solout = ez;
        }
            break;
        default:
            DebugStop();
            break;
    }
}






