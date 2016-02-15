/**
 * @file
 * @brief Afirst attempt at a multi physics hcurl/h1 formulation
 * @author Francisco Orlandini
 * @since 2015
 */


#include <iostream>
#include <fstream>
#include <string>
#include "pzquad.h"
#include "pzmaterialdata.h"
#include "TPZMatMFHCurlFran.h"

void FillPhi( const int &nScaFunc , const int &nVecFunc , const REAL detjac , const TPZVec<REAL> &pos , TPZFMatrix<REAL> &phiSca , TPZFMatrix<REAL> &dphiSca , TPZFMatrix<REAL> &phiVec);

int main(int argc, char *argv[])
{
  //integration rule settings
  const int order = 10;
  TPZIntQuad intRule;
  intRule.SetType(order, order); //gauss-legendre
  TPZVec<REAL> pos(3 , 0.);
  REAL w = 0.;
  
  
  //material settings
  TPZMatMFHCurlFran *mat = new TPZMatMFHCurlFran(1);
  TPZVec<TPZMaterialData> datavec(2);
  datavec[0].fShapeType = TPZMaterialData::EVecShape;
  datavec[0].detjac = 1;
  const int nVecFunc = 4;
  datavec[1].fShapeType = TPZMaterialData::EScalarShape;
  datavec[1].detjac = 1;
  const int nScaFunc = 4;
  
  //create matrices
  TPZFMatrix<STATE> ek , ef;
  ek.Resize(nVecFunc + nScaFunc , nVecFunc + nScaFunc);
  ef.Resize(nVecFunc + nScaFunc , 1);
  
  for (int i = 0; i < intRule.NPoints(); i++) {
    intRule.Point(i, pos, w);
    FillPhi(nScaFunc, nVecFunc, datavec[0].detjac , pos, datavec[0].phi, datavec[0].dphix, datavec[1].phi);
    mat->Contribute(datavec, w, ek, ef);
  }

  return 0;
}

void FillPhi( const int &nScaFunc , const int &nVecFunc , const REAL detjac , const TPZVec<REAL> &pos , TPZFMatrix<REAL> &phiSca , TPZFMatrix<REAL> &dphiSca , TPZFMatrix<REAL> &phiVec)
{
  phiSca.Resize(nScaFunc , 1);
  phiSca.Zero();
  dphiSca.Resize(3 , nScaFunc);
  dphiSca.Zero();
  phiVec.Resize(nVecFunc , 3);
  phiVec.Zero();
  ///////////
  TPZVec<REAL> phiXsi(2 , 0.) , phiEta ( 2 , 0. ) , dphiXsi(2 , 0.) , dphiEta ( 2 , 0. );
  phiXsi[0] = ( 1 - pos[0] ) / 2.;
  dphiXsi[0] = -0.5;
  phiXsi[1] = ( 1 + pos[0] ) / 2.;
  dphiXsi[1] = 0.5;
  phiEta[0] = ( 1 - pos[1] ) / 2.;
  dphiEta[0] = -0.5;
  phiEta[1] = ( 1 + pos[1] ) / 2.;
  dphiEta[1] = 0.5;
  
  for (int i = 0 ; i < nScaFunc ; i++) {
    phiSca = 1.;
    phiSca *= ( i % 3 ) ? phiXsi[1] : phiXsi[0];
    phiSca *= ( i / 2 ) ? phiEta[1] : phiEta[0];
  }
  
  
  return;
}
