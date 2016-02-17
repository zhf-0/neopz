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

void FillMatDataSca( const int &nScaFunc ,  TPZVec<REAL> &pos ,  TPZMaterialData &matData);

void FillMatDataVec( const int &nScaFunc , const int &nVecFunc ,  TPZVec<REAL> &pos ,  TPZMaterialData &matData);

STATE ur( const TPZVec<REAL> & x){
  //return ( 2.-imaginary*0.1 );
  return 1.;
}
STATE er( const TPZVec<REAL> & x){
  return 4.26;
}
int main(int argc, char *argv[])
{
  //integration rule settings
  const int order = 10;
  TPZIntQuad intRule;
  intRule.SetType(order, order); //gauss-legendre
  TPZVec<REAL> pos(3 , 0.);
  REAL w = 0.;
  
  
  //material settings
  REAL lambda = 1.;//1550 * 1e-6;
  REAL k0 = 2 * M_PI * 1. / lambda * M_C * sqrt(M_UZERO * M_EZERO);
  REAL kz = k0 / 4 ;
  REAL e0 = 1;
  REAL t = 0;
  REAL scale = 1;
  
  TPZMatMFHCurlFran *mat = new TPZMatMFHCurlFran(1, lambda, kz, ur,er, e0, t, scale);
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
    FillMatDataSca(nScaFunc, pos , datavec[0]);
    FillMatDataVec(nScaFunc , nVecFunc ,  pos ,  datavec[1]);
    mat->Contribute(datavec, w, ek, ef);
  }
  ek.Print("ek = ", std::cout, EMathematicaInput);
  return 0;
}

void FillMatDataSca( const int &nScaFunc ,  TPZVec<REAL> &pos , TPZMaterialData &matData)
{
  
  matData.x = pos;
  TPZFMatrix<REAL> &phiSca = matData.phi , &dphiSca = matData.dphix;
  phiSca.Resize(nScaFunc , 1);
  phiSca.Zero();
  dphiSca.Resize(3 , nScaFunc);
  dphiSca.Zero();
  ///////////
  TPZVec<REAL> phiXsi(2 , 0.) , phiEta ( 2 , 0. ) , dphiXsi(2 , 0.) , dphiEta ( 2 , 0. );
  phiXsi[0] = ( 1 - matData.x[0] ) / 2.;
  dphiXsi[0] = -0.5;
  phiXsi[1] = ( 1 + matData.x[0] ) / 2.;
  dphiXsi[1] = 0.5;
  phiEta[0] = ( 1 - matData.x[1] ) / 2.;
  dphiEta[0] = -0.5;
  phiEta[1] = ( 1 + matData.x[1] ) / 2.;
  dphiEta[1] = 0.5;
  
  phiSca(0 , 0) = phiEta[0] * phiXsi[0];
  dphiSca(0 , 0) = phiEta[0] * dphiXsi[0];
  dphiSca(1 , 0) = dphiEta[0] * phiXsi[0];
  /************************/
  phiSca(1 , 0) = phiEta[0] * phiXsi[1];
  dphiSca(0 , 1) = phiEta[0] * dphiXsi[1];
  dphiSca(1 , 1) = dphiEta[0] * phiXsi[1];
  /************************/
  phiSca(2 , 0) = phiEta[1] * phiXsi[1];
  dphiSca(0 , 2) = phiEta[1] * dphiXsi[1];
  dphiSca(1 , 2) = dphiEta[1] * phiXsi[1];
  /************************/
  phiSca(3 , 0) = phiEta[1] * phiXsi[0];
  dphiSca(0 , 3) = phiEta[1] * dphiXsi[0];
  dphiSca(1 , 3) = dphiEta[1] * phiXsi[0];
  
  return;
}

void FillMatDataVec( const int &nScaFunc , const int &nVecFunc ,  TPZVec<REAL> &pos , TPZMaterialData &matData)
{
  
  matData.x = pos;
  TPZFMatrix<REAL> &phiSca = matData.phi , &dphiSca = matData.dphix;
  phiSca.Resize(nScaFunc , 1);
  phiSca.Zero();
  dphiSca.Resize(3 , nScaFunc);
  dphiSca.Zero();
  ///////////
  TPZVec<REAL> phiXsi(2 , 0.) , phiEta ( 2 , 0. ) , dphiXsi(2 , 0.) , dphiEta ( 2 , 0. );
  phiXsi[0] = ( 1 - matData.x[0] ) / 2.;
  dphiXsi[0] = -0.5;
  phiXsi[1] = ( 1 + matData.x[0] ) / 2.;
  dphiXsi[1] = 0.5;
  phiEta[0] = ( 1 - matData.x[1] ) / 2.;
  dphiEta[0] = -0.5;
  phiEta[1] = ( 1 + matData.x[1] ) / 2.;
  dphiEta[1] = 0.5;
  
  
  phiSca(0 , 0) = phiEta[0];
  dphiSca(1 , 0) = dphiEta[0];
  /************************/
  phiSca(1 , 0) = phiXsi[1];
  dphiSca(0 , 1) = dphiXsi[1];
  /************************/
  phiSca(2 , 0) = phiEta[1];
  dphiSca(1 , 2) = dphiEta[1];
  /************************/
  phiSca(3 , 0) = phiXsi[0];
  dphiSca(0 , 3) = dphiXsi[0];
  
  matData.fNormalVec.Resize( 3, 4);
  matData.fNormalVec( 0 , 0) = 1;
  matData.fNormalVec( 1 , 0) = 0;
  matData.fNormalVec( 2 , 0) = 0;
  /*****************************/
  matData.fNormalVec( 0 , 1) = 0;
  matData.fNormalVec( 1 , 1) = 1;
  matData.fNormalVec( 2 , 1) = 0;
  /*****************************/
//  matData.fNormalVec( 0 , 2) =-1;
  matData.fNormalVec( 0 , 2) = 1;
  matData.fNormalVec( 1 , 2) = 0;
  matData.fNormalVec( 2 , 2) = 0;
  /*****************************/
  matData.fNormalVec( 0 , 3) = 0;
//  matData.fNormalVec( 1 , 3) =-1;
  matData.fNormalVec( 1 , 3) = 1;
  matData.fNormalVec( 2 , 3) = 0;
  std::pair<int , long> index;
  matData.fVecShapeIndex.Resize( nVecFunc , std::make_pair(-1,-1));
  int scaIndexes[4] = { 0 , 1 , 2 , 3};
  int vecIndexes[4] = { 0 , 1 , 2 , 3};
  
  for (int i = 0 ; i < 4; i ++) {
    matData.fVecShapeIndex[i] = std::make_pair( scaIndexes[i] , vecIndexes[i] );
  }
  
  //axes_{i,j} =  dx_i / dxchapeu_j
  matData.axes.Resize(3 , 2);
  matData.axes( 0 , 0 ) = 1;
  matData.axes( 0 , 1 ) = 0;
  
  matData.axes( 1 , 0 ) = 0;
  matData.axes( 1 , 1 ) = 1;

  matData.axes( 2 , 0 ) = 0;
  matData.axes( 2 , 1 ) = 0;
  return;
}
