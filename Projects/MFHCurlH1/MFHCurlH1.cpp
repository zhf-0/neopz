/**
 * @file
 * @brief Afirst attempt at a multi physics hcurl/h1 formulation
 * @details Adequadte for problems with longitudinal axis symmetry
 * such as some setions of waveguides (closed waveguides).
 * it uses hcurl for tranverse componentes and h1 for longitudinal components
 *
 * @author Francisco Orlandini
 * @since 2015
 */


#include <iostream>
#include <fstream>
#include <string>
#include "TPZMatMFHCurlFran.h"
#include "TPZMatHCurl2D.h"
#include "TPZMatComplexH12D.h"
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzlog.h"
#include "TPZTimer.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"
#include "TPZSBMatrixLapack.h"
#include "pzbuildmultiphysicsmesh.h"

enum meshTypeE{ createRectangular=1, createTriangular, createZigZag};

STATE ur( const TPZVec<REAL> & x){
  //return ( 2.-imaginary*0.1 );
  return 1.;
}
STATE er( const TPZVec<REAL> & x){
  return 4.26;
}
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL L, REAL theta, REAL lambda, REAL kz, REAL e0, REAL scale);

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int xDiv, const int zDiv);

int main(int argc, char *argv[])
{
  HDivPiola = 1;//use piola mapping
  TPZTimer timer;
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  
  //PARAMETROS FISICOS DO PROBLEMA
  REAL lambda = 1550*1e-9;
  REAL theta = 0.;
  REAL e0 = 1.;
  //PARAMETROS DA GEOMETRIA
  //  REAL L = 4.;
  //  REAL hDomain = L;
  //  REAL wDomain = L;
  //  REAL scale = (5.*lambda)/L;
  
  REAL L = 5*lambda;
  REAL hDomain = L;
  REAL wDomain = L;
  REAL scale = 1.;
  
  REAL w=2.*M_PI*M_C/lambda;
  REAL kZ=w*sqrt(M_UZERO*M_EZERO)/5;
  
  int pOrder = 1; //ordem polinomial de aproximacao
  int dim = 2;
  int xDiv = 1;
  int zDiv = 1;
  
  const int meshType = createTriangular;
  timer.start();
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  CreateGMesh(gmesh, meshType, hDomain, wDomain,  xDiv, zDiv);
  
  TPZCompMesh *cmesh = CMesh(gmesh, pOrder, ur , er , L , theta , lambda , kZ ,  e0 , scale); //funcao para criar a malha computacional
  bool optimizeBandwidth = false;
  TPZAnalysis an(cmesh,optimizeBandwidth);
  
  //configuracoes do objeto de analise
  
  
  TPZSBMatrixLapack<STATE> sbndmtrx( cmesh->NEquations() , cmesh->BandWidth());
  TPZSkylineNSymStructMatrix skylstr(cmesh);
  skylstr.SetNumThreads(0);
  an.SetStructuralMatrix(skylstr);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU); //caso simetrico
  an.SetSolver(step);
  
  TPZStack<std::string> scalnames, vecnames;
  vecnames.Push("absE");//setando para imprimir campoeletrico
  //vecnames.Push("solAnal");//setando para imprimir campoeletrico
  std::string plotfile= "../ValidacaoHCurlFran2EField.vtk";//arquivo de saida que estara na pasta debug
  an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
  int postProcessResolution = 2 ;//define resolucao do pos processamento
  //fim das configuracoes do objeto de analise
  
  //  const REAL k0 = 2*M_PI*freq*sqrt(M_UZERO*M_EZERO);
  //  const REAL L = wDomain ;
#ifdef DEBUG
  {
    std::ofstream file("../cmesh.txt");
    cmesh->Print(file);
  }
#endif
  // Applying the filter
  //  TPZManVector<long> ActiveEquations;
  //  FilterOthersFunctions(ActiveEquations,cmesh);
  //  an.StructMatrix()->EquationFilter().Reset();
  //  an.StructMatrix()->EquationFilter().SetActiveEquations(ActiveEquations);
  
  // Resolvendo o Sistema
  std::cout<<"entrando no assemble"<<std::endl;
  an.Assemble();
  TPZFMatrix< std::complex<double> > stiff;
  stiff = *an.Solver().Matrix().operator->();
#ifdef DEBUG
  stiff.Print("KPZglobalNed = " , std::cout , EMathematicaInput);
  std::cout<<"saindo do assemble"<<std::endl;
  std::cout<<"entrando no solver"<<std::endl;
#endif
  an.Solve();
#ifdef DEBUG
  std::cout<<"saindo do solver"<<std::endl;
#endif
  const TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
  
  solucao.Print("solucao");
  //  TPZFMatrix<STATE> sds(solucao);
  //  solucao.Print(std::cout);
  
  an.PostProcess(postProcessResolution);//realiza pos processamento*)
  
  //  for (int i = 0; i < solucao.Rows(); i++) {
  //
  //    sds.Zero();
  //    sds(i,0) = solucao.GetVal(i,0);
  //    cmesh->Solution() = sds;
  //    an.PostProcess(postProcessResolution);//realiza pos processamento*)
  //  }
  
  timer.stop();
  
  
  
  std::cout <<"Tempo de simulacao total = "<<timer.seconds()<<" s\n";
  
  //computar projecao
  
  
  std::cout << "FINISHED!" << std::endl;
  
  return 0;
}

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int xDiv, const int zDiv)
{
  
  TPZManVector<int,3> nx(3,0);
  TPZManVector<REAL,3> llCoord(3,0.) , ulCoord(3,0.) , urCoord(3,0.) , lrCoord(3,0.);
  llCoord[0] = 0.;
  llCoord[1] = 0.;
  
  ulCoord[0] = 0.;
  ulCoord[1] = hDomain;
  
  urCoord[0] = wDomain;
  urCoord[1] = hDomain;
  
  lrCoord[0] = wDomain;
  lrCoord[1] = 0.;
  
  nx[0]=xDiv;
  nx[1]=zDiv;
  int numl = 1;
  TPZGenGrid *gengrid = NULL;
  switch (meshType) {
    case createRectangular:
    {
      REAL rot = 0.0;
      gengrid = new TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
      gengrid->SetElementType(EQuadrilateral);
#ifdef DEBUG
      std::ofstream outTxt("../gmeshRectangular.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshRectangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
#endif
    }
      break;
    case createTriangular:
    {
      REAL rot = 0.0;
      gengrid = new TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
      gengrid->SetElementType(ETriangle);
#ifdef DEBUG
      std::ofstream outTxt("../gmeshTriangular.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshTriangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
#endif
    }
      break;
    case createZigZag:
    {
      REAL rot = 0.0;
      gengrid = new TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
      gengrid->SetElementType(EQuadrilateral);
      gengrid->SetZigZagPattern();
#ifdef DEBUG
      std::ofstream outTxt("../gmeshZigZag.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshZigZag.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
#endif
    }
      break;
    default:
      DebugStop();
      break;
  }
  gmesh = new TPZGeoMesh();
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno dirichlet)
  gengrid->Read(gmesh , matId);
  
  
  gengrid->SetBC(gmesh, ulCoord, llCoord, bc0);
  gengrid->SetBC(gmesh, urCoord, ulCoord, bc0);
  gengrid->SetBC(gmesh, lrCoord, urCoord, bc0);
  gengrid->SetBC(gmesh, llCoord, lrCoord, bc0);
  
  //gmesh->ResetConnectivities();
  
  gmesh->BuildConnectivity();
#ifdef DEBUG
  std::ofstream outTxt , outVtk;
  switch (meshType) {
    case createRectangular:
    {
      outTxt.open("../gmeshRectangular.txt"); //define arquivo de saida para impressao da malha no
      outVtk.open("../gmeshRectangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      
    }
      break;
    case createTriangular:
    {
      outTxt.open("../gmeshTriangular.txt"); //define arquivo de saida para impressao da malha no
      outVtk.open("../gmeshTriangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      
    }
      break;
    case createZigZag:
    {
      outTxt.open("../gmeshZigZag.txt"); //define arquivo de saida para impressao da malha no
      outVtk.open("../gmeshZigZag.vtk"); //define arquivo de saida para impressao da malha no paraview
    }
      break;
    default:
      DebugStop();
      break;
  }
  gmesh->Print(outTxt);
  outTxt.close();
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVtk, true); //imprime a malha no formato vtk
  outVtk.close();
#endif
  delete gengrid;
}

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL L, REAL theta, REAL lambda, REAL kz, REAL e0, REAL scale)
{
  
  const int dim = 2; //dimensao do problema
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno dirichlet)
  enum{ dirichlet = 0, neumann, mixed}; //tipo da condicao de contorno do problema
  // Criando material

  
  ///criar malha computacional H1
  
  TPZCompMesh * cmeshH1 = new TPZCompMesh(gmesh);
  cmeshH1->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmeshH1->SetDimModel(dim);//seta dimensao do modelo
  // Inserindo material na malha
  REAL kZ = 1.;
  TPZMatComplexH12D *matH1 = new TPZMatComplexH12D( matId, lambda, kZ , ur, er , e0, theta , scale);
  cmeshH1->InsertMaterialObject(matH1);
		
  ///electrical conductor boundary conditions
  TPZFNMatrix<1,STATE> val1(1,1,0.), val2(1,1,0.);
  
  TPZMaterial * BCondH1Dir = matH1->CreateBC(matH1, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
  
  cmeshH1->InsertMaterialObject(BCondH1Dir);//insere material na malha
  //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
  cmeshH1->AutoBuild();
  ///criar malha computacional HCurl
  
  TPZCompMesh * cmeshHCurl = new TPZCompMesh(gmesh);
  cmeshHCurl->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmeshHCurl->SetDimModel(dim);//seta dimensao do modelo
  // Inserindo material na malha
  TPZMatHCurl2D *matHCurl = new TPZMatHCurl2D(matId , lambda , ur , er , e0 , theta, scale);
  cmeshHCurl->InsertMaterialObject(matHCurl);
  
  val1( 0, 0 ) = 0.;
  val2( 0, 0 ) = 0.;
  TPZMaterial * BCondHCurlDir = matH1->CreateBC(matHCurl, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
  
  cmeshHCurl->InsertMaterialObject(BCondHCurlDir);//insere material na malha
  
  cmeshHCurl->SetAllCreateFunctionsHDiv();//define espaco de aproximacao
  cmeshHCurl->AutoBuild();
  
  TPZAdmChunkVector< TPZCompEl* > elVec = cmeshHCurl->ElementVec();

  for (int i = 0; i < cmeshHCurl->NElements(); i++) {
    TPZCompElHDiv < pzshape::TPZShapeQuad > *el = dynamic_cast<TPZCompElHDiv <pzshape::TPZShapeQuad > *>( elVec[i] );
    if ( el == NULL) {
      continue;
    }
    el->SetSideOrient(4,  1);
    el->SetSideOrient(5,  1);
    el->SetSideOrient(6, -1);
    el->SetSideOrient(7, -1);
  }

  for (int i = 0; i < cmeshHCurl->NElements(); i++) {
    TPZCompElHDiv < pzshape::TPZShapeTriang > *el = dynamic_cast<TPZCompElHDiv <pzshape::TPZShapeTriang > *>( elVec[i] );
    if ( el == NULL) {
      continue;
    }
    if ( i % 2 == 1) {
      el->SetSideOrient(3, -1);
      el->SetSideOrient(4, -1);
      el->SetSideOrient(5, -1);
    } 
  }


  if (pOrder == 1) {
    //cmesh->CleanUpUnconnectedNodes();
    TPZCreateApproximationSpace::MakeRaviartThomas(*cmeshHCurl);
    cmeshHCurl->CleanUpUnconnectedNodes();
  }
  else{//for now only lowest order elements are avaliable
    DebugStop();
  }

  
  TPZVec<TPZCompMesh *> meshVec(2);
  meshVec[ 0 ] = cmeshH1;
  meshVec[ 1 ] = cmeshHCurl;
  gmesh->ResetReference();

  TPZCompMesh *cmeshMF = new TPZCompMesh( gmesh );
  cmeshMF->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmeshMF->SetDimModel(dim);
  cmeshMF->SetAllCreateFunctionsMultiphysicElem();
  
  TPZMatMFHCurlFran *matMultiPhysics = new TPZMatMFHCurlFran(matId , lambda , kz , ur , er , e0 , theta, scale)  ;//criando material que implementa a
  //formulacao fraca do problema de validacao
  
  val1( 0, 0 ) = 0.;
  val2( 0, 0 ) = 0.;
  TPZMaterial * BCondMFDir = matMultiPhysics->CreateBC(matMultiPhysics, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
  
  
  cmeshMF->InsertMaterialObject(matMultiPhysics);
  cmeshMF->InsertMaterialObject(BCondMFDir);//insere material na malha
  
  //cmeshMF->AutoBuild();
  //Creating multiphysic elements containing skeletal elements.
  TPZBuildMultiphysicsMesh::AddElements(meshVec, cmeshMF);
  TPZBuildMultiphysicsMesh::AddConnects(meshVec, cmeshMF);
  TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVec, cmeshMF);

  return cmeshMF;

}


//#include "pzquad.h"
//#include "pzmaterialdata.h"
//void FillMatDataSca( const int &nScaFunc ,  TPZVec<REAL> &pos ,  TPZMaterialData &matData);
//
//void FillMatDataVec( const int &nScaFunc , const int &nVecFunc ,  TPZVec<REAL> &pos ,  TPZMaterialData &matData);
//
//
//
//
//
//int main(int argc, char *argv[])
//{
//  //integration rule settings
//  const int order = 10;
//  TPZIntQuad intRule;
//  intRule.SetType(order, order); //gauss-legendre
//  TPZVec<REAL> pos(3 , 0.);
//  REAL w = 0.;
//  
//  
//  //material settings
//  REAL lambda = 1.;//1550 * 1e-6;
//  REAL k0 = 2 * M_PI * 1. / lambda * M_C * sqrt(M_UZERO * M_EZERO);
//  REAL kz = k0 / 4 ;
//  REAL e0 = 1;
//  REAL t = 0;
//  REAL scale = 1;
//  
//  TPZMatMFHCurlFran *mat = new TPZMatMFHCurlFran(1, lambda, kz, ur,er, e0, t, scale);
//  TPZVec<TPZMaterialData> datavec(2);
//  datavec[0].fShapeType = TPZMaterialData::EVecShape;
//  datavec[0].detjac = 1;
//  const int nVecFunc = 4;
//  datavec[1].fShapeType = TPZMaterialData::EScalarShape;
//  datavec[1].detjac = 1;
//  const int nScaFunc = 4;
//  
//  //create matrices
//  TPZFMatrix<STATE> ek , ef;
//  ek.Resize(nVecFunc + nScaFunc , nVecFunc + nScaFunc);
//  ef.Resize(nVecFunc + nScaFunc , 1);
//  
//  for (int i = 0; i < intRule.NPoints(); i++) {
//    intRule.Point(i, pos, w);
//    FillMatDataSca(nScaFunc, pos , datavec[0]);
//    FillMatDataVec(nScaFunc , nVecFunc ,  pos ,  datavec[1]);
//    mat->Contribute(datavec, w, ek, ef);
//  }
//  ek.Print("ek = ", std::cout, EMathematicaInput);
//  return 0;
//}
//
//void FillMatDataSca( const int &nScaFunc ,  TPZVec<REAL> &pos , TPZMaterialData &matData)
//{
//  
//  matData.x = pos;
//  TPZFMatrix<REAL> &phiSca = matData.phi , &dphiSca = matData.dphix;
//  phiSca.Resize(nScaFunc , 1);
//  phiSca.Zero();
//  dphiSca.Resize(3 , nScaFunc);
//  dphiSca.Zero();
//  ///////////
//  TPZVec<REAL> phiXsi(2 , 0.) , phiEta ( 2 , 0. ) , dphiXsi(2 , 0.) , dphiEta ( 2 , 0. );
//  phiXsi[0] = ( 1 - matData.x[0] ) / 2.;
//  dphiXsi[0] = -0.5;
//  phiXsi[1] = ( 1 + matData.x[0] ) / 2.;
//  dphiXsi[1] = 0.5;
//  phiEta[0] = ( 1 - matData.x[1] ) / 2.;
//  dphiEta[0] = -0.5;
//  phiEta[1] = ( 1 + matData.x[1] ) / 2.;
//  dphiEta[1] = 0.5;
//  
//  phiSca(0 , 0) = phiEta[0] * phiXsi[0];
//  dphiSca(0 , 0) = phiEta[0] * dphiXsi[0];
//  dphiSca(1 , 0) = dphiEta[0] * phiXsi[0];
//  /************************/
//  phiSca(1 , 0) = phiEta[0] * phiXsi[1];
//  dphiSca(0 , 1) = phiEta[0] * dphiXsi[1];
//  dphiSca(1 , 1) = dphiEta[0] * phiXsi[1];
//  /************************/
//  phiSca(2 , 0) = phiEta[1] * phiXsi[1];
//  dphiSca(0 , 2) = phiEta[1] * dphiXsi[1];
//  dphiSca(1 , 2) = dphiEta[1] * phiXsi[1];
//  /************************/
//  phiSca(3 , 0) = phiEta[1] * phiXsi[0];
//  dphiSca(0 , 3) = phiEta[1] * dphiXsi[0];
//  dphiSca(1 , 3) = dphiEta[1] * phiXsi[0];
//  
//  return;
//}
//
//void FillMatDataVec( const int &nScaFunc , const int &nVecFunc ,  TPZVec<REAL> &pos , TPZMaterialData &matData)
//{
//  
//  matData.x = pos;
//  TPZFMatrix<REAL> &phiSca = matData.phi , &dphiSca = matData.dphix;
//  phiSca.Resize(nScaFunc , 1);
//  phiSca.Zero();
//  dphiSca.Resize(3 , nScaFunc);
//  dphiSca.Zero();
//  ///////////
//  TPZVec<REAL> phiXsi(2 , 0.) , phiEta ( 2 , 0. ) , dphiXsi(2 , 0.) , dphiEta ( 2 , 0. );
//  phiXsi[0] = ( 1 - matData.x[0] ) / 2.;
//  dphiXsi[0] = -0.5;
//  phiXsi[1] = ( 1 + matData.x[0] ) / 2.;
//  dphiXsi[1] = 0.5;
//  phiEta[0] = ( 1 - matData.x[1] ) / 2.;
//  dphiEta[0] = -0.5;
//  phiEta[1] = ( 1 + matData.x[1] ) / 2.;
//  dphiEta[1] = 0.5;
//  
//  
//  phiSca(0 , 0) = phiEta[0];
//  dphiSca(1 , 0) = dphiEta[0];
//  /************************/
//  phiSca(1 , 0) = phiXsi[1];
//  dphiSca(0 , 1) = dphiXsi[1];
//  /************************/
//  phiSca(2 , 0) = phiEta[1];
//  dphiSca(1 , 2) = dphiEta[1];
//  /************************/
//  phiSca(3 , 0) = phiXsi[0];
//  dphiSca(0 , 3) = dphiXsi[0];
//  
//  matData.fNormalVec.Resize( 3, 4);
//  matData.fNormalVec( 0 , 0) = 1;
//  matData.fNormalVec( 1 , 0) = 0;
//  matData.fNormalVec( 2 , 0) = 0;
//  /*****************************/
//  matData.fNormalVec( 0 , 1) = 0;
//  matData.fNormalVec( 1 , 1) = 1;
//  matData.fNormalVec( 2 , 1) = 0;
//  /*****************************/
////  matData.fNormalVec( 0 , 2) =-1;
//  matData.fNormalVec( 0 , 2) = 1;
//  matData.fNormalVec( 1 , 2) = 0;
//  matData.fNormalVec( 2 , 2) = 0;
//  /*****************************/
//  matData.fNormalVec( 0 , 3) = 0;
////  matData.fNormalVec( 1 , 3) =-1;
//  matData.fNormalVec( 1 , 3) = 1;
//  matData.fNormalVec( 2 , 3) = 0;
//  std::pair<int , long> index;
//  matData.fVecShapeIndex.Resize( nVecFunc , std::make_pair(-1,-1));
//  int scaIndexes[4] = { 0 , 1 , 2 , 3};
//  int vecIndexes[4] = { 0 , 1 , 2 , 3};
//  
//  for (int i = 0 ; i < 4; i ++) {
//    matData.fVecShapeIndex[i] = std::make_pair( scaIndexes[i] , vecIndexes[i] );
//  }
//  
//  //axes_{i,j} =  dx_i / dxchapeu_j
//  matData.axes.Resize(3 , 2);
//  matData.axes( 0 , 0 ) = 1;
//  matData.axes( 0 , 1 ) = 0;
//  
//  matData.axes( 1 , 0 ) = 0;
//  matData.axes( 1 , 1 ) = 1;
//
//  matData.axes( 2 , 0 ) = 0;
//  matData.axes( 2 , 1 ) = 0;
//  return;
//}
