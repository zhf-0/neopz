//
// Created by Francisco Teixeira Orlandini on 3/6/18.
//

#include <pzaxestools.h>
#include <pzvec_extras.h>
#include "TPZMatWaveguidePml.h"

TPZMatWaveguidePml::TPZMatWaveguidePml(const int id,const TPZMatModalAnalysis &mat,
                                       const bool &att_x, REAL &pmlBeginX,
                                       const bool &att_y, REAL &pmlBeginY,
                                       const REAL &alphaMax, const REAL &d) :
        TPZMatModalAnalysis(mat), fAttX (att_x), fAttY (att_y),
        fPmlBeginX (pmlBeginX), fPmlBeginY (pmlBeginY),
        fAlphaMax (alphaMax), fD (d)
{
    this->SetId(id);
    if(fAlphaMax < 0) DebugStop(); //for the attenuation to happen
                                   // this value must be positive
    if(fD < 0) DebugStop(); // pml width must be positive
    if(!fAttX && !fAttY) DebugStop();//a pml that doesnt attenuate
                                     // in any direction?
}

TPZMatWaveguidePml::~TPZMatWaveguidePml(){
}

void TPZMatWaveguidePml::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /*****************CALCULATE S PML PARAMETERS*************************
     * In the current application, the waveguide's cross section is always
     * in the xy-plane. Therefore, sz will always be unity, and omitted for
     * the folllowing calculations. The same principle applies, for instance,
     * for the z-component of the hcurl functions, the x and y components of
     * their curl and so on.
     */
    TPZManVector<REAL,3> x = datavec[ h1meshindex ].XCenter;
    STATE sx = 1, sy = 1;
    if(fAttX){
        sx = 1. - imaginary * fAlphaMax * ((x[0]-fPmlBeginX) / fD ) * ((x[0]-fPmlBeginX) / fD );
    }
    if(fAttY){
        sy = 1. - imaginary * fAlphaMax * ((x[1]-fPmlBeginY) / fD ) * ((x[1]-fPmlBeginY) / fD );
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

    const REAL k0 = fScaleFactor * 2*M_PI/fLambda;
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/

    const int nHCurlFunctions  = phiHCurl.Rows();
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = h1meshindex * nHCurlFunctions;
    const int firstHCurl = hcurlmeshindex * nH1Functions;

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE curlIzdotCurlJz = 0.;
            curlIzdotCurlJz += curlPhi(2 , iVec) * curlPhi(2 , jVec);
            STATE phiIdotPhiJx = phiHCurl(iVec , 0) * phiHCurl(jVec , 0);
            STATE phiIdotPhiJy = phiHCurl(iVec , 1) * phiHCurl(jVec , 1);

            STATE stiffAtt = 0.;
            stiffAtt = 1./(sx * sy * fUr) * curlIzdotCurlJz;
            stiffAtt -= k0 * k0 * (fEr * sy/sx) * phiIdotPhiJx;
            stiffAtt -= k0 * k0 * (fEr * sx/sy) * phiIdotPhiJy;
            STATE stiffBtt = 0.;
            stiffBtt += (sy/(sx * fUr)) * phiIdotPhiJx;
            stiffBtt += (sx/(sy * fUr)) * phiIdotPhiJy;
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
            STATE phiVecDotGradPhiScax = phiHCurl(iVec , 0) * gradPhiH1(jSca , 0);
            STATE phiVecDotGradPhiScay = phiHCurl(iVec , 1) * gradPhiH1(jSca , 1);

            STATE stiffBzt = 0.;
            stiffBzt += (sy/(sx * fUr)) * phiVecDotGradPhiScax;
            stiffBzt += (sx/(sy * fUr)) * phiVecDotGradPhiScay;
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
            STATE phiVecDotGradPhiScax = phiHCurl(jVec , 0) * gradPhiH1(iSca , 0);
            STATE phiVecDotGradPhiScay = phiHCurl(jVec , 1) * gradPhiH1(iSca , 1);

            STATE stiffBtz = 0.;
            stiffBtz += (sy/(sx * fUr)) * phiVecDotGradPhiScax;
            stiffBtz += (sx/(sy * fUr)) * phiVecDotGradPhiScay;
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
            STATE gradPhiScaDotGradPhiScax = gradPhiH1(iSca , 0) * gradPhiH1(jSca , 0);
            STATE gradPhiScaDotGradPhiScay = gradPhiH1(iSca , 1) * gradPhiH1(jSca , 1);

            STATE stiffBzz = 0.;
            stiffBzz +=  (sy/(sx * fUr)) * gradPhiScaDotGradPhiScax;
            stiffBzz +=  (sx/(sy * fUr)) * gradPhiScaDotGradPhiScay;
            stiffBzz -=  k0 * k0 * (fEr * sx * sy) * phiH1( iSca , 0 ) * phiH1( jSca , 0 );

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



void TPZMatWaveguidePml::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
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
            /*****************CALCULATE S PML PARAMETERS*************************/
            TPZManVector<REAL,3> x = datavec[ h1meshindex ].x;
            REAL sx = 1, sy = 1;
            Solout.Resize(2);
            if(fAttX){
                sx = 1 + fAlphaMax * ((x[0]-fPmlBeginX) / fD ) * ((x[0]-fPmlBeginX) / fD );
            }
            if(fAttY){
                sy = 1 + fAlphaMax * ((x[1]-fPmlBeginY) / fD ) * ((x[1]-fPmlBeginY) / fD );
            }
            Solout[0] = fEr * sx;
            Solout[1] = fEr * sy;
            break;
        }
        default:
            DebugStop();
            break;
    }
}

int TPZMatWaveguidePml::IntegrationRuleOrder(int elPMaxOrder) const
{
    return 4+elPMaxOrder*2;
}
