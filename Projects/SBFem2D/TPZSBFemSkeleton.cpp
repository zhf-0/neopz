//
//  TPZSBFemSkeleton.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/27/16.
//
//

#include "TPZSBFemSkeleton.hpp"
#include "pzmaterial.h"



void TPZSBFemSkeleton::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                         TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                         REAL &detjac, TPZFMatrix<REAL> &jacinv,
                                         TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidaxes){
    DebugStop();
}

void TPZSBFemSkeleton::ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data){
    
    TPZGeoEl * ref = this->Reference();
    if (!ref){
        PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
        return;
    }//if
    
    ref->X(intpoint,data.x);
    long nshape = data.phi.Rows()/2;
    TPZFNMatrix<30,REAL> phi(nshape,1),dphi(1,nshape);
    this->Shape(intpoint,phi,dphi);
    for (int ish=0; ish<nshape; ish++) {
        data.phi(ish,0) = phi(ish,0);
        data.phi(ish+nshape,0) = phi(ish,0);
        data.dphi(0,ish) = dphi(0,ish);
        data.dphi(1,ish) = 0.;
        data.dphi(0,ish+nshape) = 0.;
        data.dphi(1,ish+nshape) = data.phi(ish,0);
    }
    this->Convert2Axes(data.dphi, data.jacinv, data.dphix);
    
    
}
