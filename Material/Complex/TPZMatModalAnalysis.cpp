#include "TPZMatModalAnalysis.h"
#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "pzbndcond.h"

#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif




TPZMatModalAnalysis::TPZMatModalAnalysis(int id, REAL lambda, const STATE &ur, const STATE &er ,
                                        const REAL &scale) :
TPZVecL2(id), fUr(ur), fEr(er), fScaleFactor(scale)
{
    isTesting = false;
    fAssembling = NDefined;
    fLambda = lambda;
}

TPZMatModalAnalysis::TPZMatModalAnalysis(int id) : TPZVecL2(id), fUr(1.0),
fEr(1.0) , fScaleFactor(1.)
{
    isTesting = false;
    fAssembling = NDefined;
    fLambda = 1.55e-9;
}

/** @brief Default constructor */
TPZMatModalAnalysis::TPZMatModalAnalysis() : TPZVecL2(), fUr(1.0),
fEr(1.0) , fScaleFactor(1.)
{
    isTesting = false;
    fAssembling = NDefined;
    fLambda=1.55e-9;
}


TPZMatModalAnalysis::TPZMatModalAnalysis(const TPZMatModalAnalysis &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr), fScaleFactor(mat.fScaleFactor)
{
    isTesting = false;
    fAssembling = NDefined;
    fLambda = mat.fLambda;
}

TPZMatModalAnalysis::~TPZMatModalAnalysis()
{
    
}
#ifdef PZDEBUG
void TPZMatModalAnalysis::ContributeValidateFunctions(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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
    
    
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHCurlFunctions  = phiHCurl.Rows();
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = h1meshindex * nHCurlFunctions;
    const int firstHCurl = hcurlmeshindex * nH1Functions;
    
    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE curlIdotCurlJ = 0.;
            curlIdotCurlJ += std::conj( curlPhi(0 , iVec) ) * curlPhi(0 , jVec);
            curlIdotCurlJ += std::conj( curlPhi(1 , iVec) ) * curlPhi(1 , jVec);
            curlIdotCurlJ += std::conj( curlPhi(2 , iVec) ) * curlPhi(2 , jVec);
            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += std::conj( phiHCurl(iVec , 0) ) * phiHCurl(jVec , 0);
            phiIdotPhiJ += std::conj( phiHCurl(iVec , 1) ) * phiHCurl(jVec , 1);
            phiIdotPhiJ += std::conj( phiHCurl(iVec , 2) ) * phiHCurl(jVec , 2);
            if (test == curl) {
                ek( firstHCurl + iVec , firstHCurl + jVec ) += curlIdotCurlJ * weight ;
            }
            else if( test == dotVec){
                ek( firstHCurl + iVec , firstHCurl + jVec ) += phiIdotPhiJ * weight ;
            }
        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE phiVecDotGradPhiSca = 0.;
            
            phiVecDotGradPhiSca += std::conj( phiHCurl(iVec , 0) ) * gradPhiH1(jSca , 0);
            phiVecDotGradPhiSca += std::conj( phiHCurl(iVec , 1) ) * gradPhiH1(jSca , 1);
            phiVecDotGradPhiSca += std::conj( phiHCurl(iVec , 2) ) * gradPhiH1(jSca , 2);
            if ( test == mixed) {
                ek( firstHCurl + iVec , firstH1 + jSca ) += phiVecDotGradPhiSca * weight ;
            }
        }
    }
    for (int iSca = 0; iSca < nH1Functions; iSca++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE phiVecDotGradPhiSca = 0.;
            phiVecDotGradPhiSca += phiHCurl(jVec , 0) * std::conj( gradPhiH1(iSca , 0) );
            phiVecDotGradPhiSca += phiHCurl(jVec , 1) * std::conj( gradPhiH1(iSca , 1) );
            phiVecDotGradPhiSca += phiHCurl(jVec , 2) * std::conj( gradPhiH1(iSca , 2) );
            
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
#endif

void TPZMatModalAnalysis::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatModalAnalysis::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
#ifdef PZDEBUG
    isTesting = false;
    if( isTesting == true ){
        ContributeValidateFunctions(datavec, weight, ek, ef);
        return;
    }
#endif
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
    const REAL k0 = fScaleFactor * 2*M_PI/fLambda;
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    
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
            
            stiffAtt = 1./fUr * curlIdotCurlJ;
            stiffAtt -= k0 * k0 * fEr * phiIdotPhiJ;
            stiffBtt = 1./fUr * phiIdotPhiJ;
            if (this->fAssembling == A) {
              ek( firstHCurl + iVec , firstHCurl + jVec ) += stiffAtt * weight ;
            }
            else if (this->fAssembling == B){
              ek( firstHCurl + iVec , firstHCurl + jVec ) += stiffBtt * weight ;
            }
            else{
                DebugStop();
            }
            
        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE stiffBzt = 0.;
            STATE phiVecDotGradPhiSca = 0.;
            
            phiVecDotGradPhiSca += phiHCurl(iVec , 0) * gradPhiH1(jSca , 0);
            phiVecDotGradPhiSca += phiHCurl(iVec , 1) * gradPhiH1(jSca , 1);
            phiVecDotGradPhiSca += phiHCurl(iVec , 2) * gradPhiH1(jSca , 2);
            
            stiffBzt = 1./fUr * phiVecDotGradPhiSca;
            if (this->fAssembling == A) {
                ek( firstHCurl + iVec , firstH1 + jSca ) += 0. ;
            }
            else if (this->fAssembling == B){
                ek( firstHCurl + iVec , firstH1 + jSca ) += stiffBzt * weight ;
            }
            else{
                DebugStop();
            }
        }
    }
    for (int iSca = 0; iSca < nH1Functions; iSca++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE phiVecDotGradPhiSca = 0.;
            STATE stiffBtz = 0.;
            phiVecDotGradPhiSca += phiHCurl(jVec , 0) * gradPhiH1(iSca , 0);
            phiVecDotGradPhiSca += phiHCurl(jVec , 1) * gradPhiH1(iSca , 1);
            phiVecDotGradPhiSca += phiHCurl(jVec , 2) * gradPhiH1(iSca , 2);
            stiffBtz = 1./fUr * phiVecDotGradPhiSca;
            if (this->fAssembling == A) {
                ek( firstH1 + iSca , firstHCurl +  jVec) += 0. ;
            }
            else if (this->fAssembling == B){
                ek( firstH1 + iSca , firstHCurl +  jVec ) += stiffBtz * weight ;
            }
            else{
                DebugStop();
            }
        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE gradPhiScaDotGradPhiSca = 0.;
            STATE stiffBzz = 0.;
            gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 0) * gradPhiH1(jSca , 0);
            gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 1) * gradPhiH1(jSca , 1);
            gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 2) * gradPhiH1(jSca , 2);
            
            stiffBzz =  1./fUr * gradPhiScaDotGradPhiSca;
            stiffBzz -=  k0 * k0 * fEr * phiH1( iSca , 0 ) * phiH1( jSca , 0 );
			
            if (this->fAssembling == A) {
                ek( firstH1 + iSca , firstH1 + jSca) += 0. ;
            }
            else if (this->fAssembling == B){
                ek( firstH1 + iSca , firstH1 + jSca) += stiffBzz * weight ;
            }
            else{
                DebugStop();
            }
        }
    }
}

void TPZMatModalAnalysis::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatModalAnalysis::ContributeForcingRTBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    DebugStop();
}

void TPZMatModalAnalysis::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

void TPZMatModalAnalysis::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
#ifdef PZDEBUG
    if( isTesting ) return;
#endif
	
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
                if(this->fAssembling == B || this->fAssembling == A){
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

void TPZMatModalAnalysis::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatModalAnalysis::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

int TPZMatModalAnalysis::IntegrationRuleOrder(int elPMaxOrder) const
{
    return elPMaxOrder*2;
}


int TPZMatModalAnalysis::VariableIndex(const std::string &name)
{
    if( strcmp(name.c_str(), "Et") == 0) return 0;
    if( strcmp(name.c_str(), "Ez") == 0) return 1;
    if( strcmp(name.c_str(), "Material") == 0) return 2;
    DebugStop();
    return 1;
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatModalAnalysis::NSolutionVariables(int var)
{
    switch (var) {
        case 0: //Et
            return 2;
            break;
        case 1://Ez
            return 1;
        case 2://material
            return 2;
        default:
            DebugStop();
            break;
    }
    return 1;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatModalAnalysis::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    
    DebugStop();
}


/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatModalAnalysis::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    
    TPZVec<STATE> et(3,0.);
    TPZVec<STATE> ez(1,0.);
	
	et = datavec[ hcurlmeshindex ].sol[0];
    ez = datavec[ h1meshindex ].sol[0];
    switch (var) {
        case 0:{//et
            Solout = et;
            break;
        }
        case 1:{//ez
            Solout = ez;
            break;
        }

        case 2:{//material
            Solout.Resize(2);
            Solout[0] = fEr;
            Solout[1] = fEr;
            break;
        }
        default:
            DebugStop();
            break;
    }
}






