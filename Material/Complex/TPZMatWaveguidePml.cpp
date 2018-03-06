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
    /*****************CALCULATE S PML PARAMETERS*************************/
    TPZManVector<REAL,3> x = datavec[ h1meshindex ].x;
    STATE sx = 1, sy = 1 , sz = 1;
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
    TPZFNMatrix<3,STATE> gradPhiH1(phiH1.Rows() , 3 , 0.);
    for ( int iFunc = 0 ; iFunc < phiH1.Rows(); iFunc++ ) {

        gradPhiH1 ( iFunc , 0 ) = dphiH1 ( 0 , iFunc );
        gradPhiH1 ( iFunc , 1 ) = dphiH1 ( 1 , iFunc );
    }

    /*********************CREATE HCURL FUNCTIONS*****************************/
    TPZFNMatrix< 36 , REAL > phiHCurlAxes = datavec[ hcurlmeshindex ].phi;
    TPZFNMatrix<40,REAL> curlPhiDAxes = datavec[ hcurlmeshindex ].dphix;

    TPZFNMatrix<20,REAL> phiHCurl,curlPhi;

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
            STATE stiffAtt = 0.;
            STATE stiffBtt = 0.;
            STATE curlIdotCurlJ = 0.;
            //curlIdotCurlJ += curlPhi(0 , iVec) * curlPhi(0 , jVec) * sx / ( sy * sz );//ok
            //curlIdotCurlJ += curlPhi(1 , iVec) * curlPhi(1 , jVec) * sy / ( sz * sx );//ok
            curlIdotCurlJ += curlPhi(2 , iVec) * curlPhi(2 , jVec) * sz / ( sx * sy );//ok
            STATE phiIdotPhiJx = 0.;
            phiIdotPhiJx += phiHCurl(iVec , 0) * phiHCurl(jVec , 0);
            STATE phiIdotPhiJy = 0.;
            phiIdotPhiJy += phiHCurl(iVec , 1) * phiHCurl(jVec , 1);
            STATE phiIdotPhiJz = 0.;
            phiIdotPhiJz += phiHCurl(iVec , 2) * phiHCurl(jVec , 2);

            stiffAtt = 1./fUr * curlIdotCurlJ;
            stiffAtt -= k0 * k0 * fEr * (phiIdotPhiJx * ( sy * sz ) / sx);
            stiffAtt -= k0 * k0 * fEr * (phiIdotPhiJy * ( sz * sx ) / sy);
            //stiffAtt -= k0 * k0 * fEr * (phiIdotPhiJz * ( sx * sy ) / sz);
            stiffBtt = 1./fUr * phiIdotPhiJx* sx / ( sy * sz );
            stiffBtt = 1./fUr * phiIdotPhiJy* sy / ( sz * sx );
            //stiffBtt = 1./fUr * phiIdotPhiJz* sz / ( sx * sy );
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

            phiVecDotGradPhiSca += phiHCurl(iVec , 0) * gradPhiH1(jSca , 0) * sx / ( sy * sz );
            phiVecDotGradPhiSca += phiHCurl(iVec , 1) * gradPhiH1(jSca , 1) * sy / ( sz * sx );
            //phiVecDotGradPhiSca += phiHCurl(iVec , 2) * gradPhiH1(jSca , 2) * sz / ( sx * sy );

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
            phiVecDotGradPhiSca += phiHCurl(jVec , 0) * gradPhiH1(iSca , 0) * sx / ( sy * sz );
            phiVecDotGradPhiSca += phiHCurl(jVec , 1) * gradPhiH1(iSca , 1) * sy / ( sz * sx );
            //phiVecDotGradPhiSca += phiHCurl(jVec , 2) * gradPhiH1(iSca , 2) * sz / ( sx * sy );
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
            gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 0) * gradPhiH1(jSca , 0) * sx / ( sy * sz );
            gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 1) * gradPhiH1(jSca , 1) * sy / ( sz * sx );
            //gradPhiScaDotGradPhiSca += gradPhiH1(iSca , 2) * gradPhiH1(jSca , 2) * sz / ( sx * sy );

            stiffBzz =  1./fUr * gradPhiScaDotGradPhiSca;
            stiffBzz -=  k0 * k0 * fEr * phiH1( iSca , 0 ) * phiH1( jSca , 0 ) * (sx * sy) / sz;

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
            REAL sx = 1, sy = 1 , sz = 1;
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
