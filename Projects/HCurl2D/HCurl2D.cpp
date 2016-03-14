/**
 * @file
 * @brief A second contact with complex state-variables using HCurl functions. An electromagnetic wave
 * incides in a dielectric slab with a metal backing. Its
 * reflection coefficient is analised as function of the
 * angle of incidence
 * @author Francisco Orlandini
 * @since 2015
 */


#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZMatHCurl2D.h"
#include "pzlog.h"
#include "TPZTimer.h"

#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"
#include "TPZSkylineNSymStructMatrix.h"


#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzgengrid.h"
#include "tpzhierarquicalgrid.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"

enum meshTypeE{ createRectangular=1, createTriangular, createZigZag};

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int xDiv, const int zDiv);

/**
 * @brief Creates cmesh
 * @note All the relevant parameters are arguments of this function (except constants)
 * @param gmesh the geometric mesh
 * @param pOrder polynomial approximation order
 * @param ur relative permeability of the dielectric
 * @param er relative permittivity of the dielectric
 * @param freq frequency of the plane-wave
 */
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL L, REAL theta, REAL lambda,REAL e0, REAL scale);

/**
 * @brief Implements permeability of possibly inhomongeneous media
 * @param x spatial coordinates
 */
static inline STATE urSubs(const TPZVec<REAL> &x);

/**
 * @brief Implements permittivity of possibly inhomongeneous media
 * @param x spatial coordinates
 */
static inline STATE erSubs(const TPZVec<REAL> &x);

void FilterOthersFunctions(TPZManVector<long> &active, TPZCompMesh * fcmesh);

void ExportCSV(const TPZFMatrix<STATE> matriz, const char* name);

void solAnal(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
	REAL x = loc[0];
	REAL y = loc[1];
	REAL theta = 0.0;
	REAL e0 = 1.;
	REAL lambda = 1550*1e-9;
	REAL c = 3e8;
	REAL L = 5*lambda;
	STATE gamma = imaginary *2.*M_PI*c/lambda*sqrt(M_UZERO*M_EZERO*urSubs(loc)*erSubs(loc));
	STATE e_analytic = e0 * exp( -1. * gamma * ( (L - x) * cos(theta) + y * sin(theta) ) ) ;
	result[2] =  e_analytic;
}

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
	
//	REAL w=2.*M_PI*M_C/lambda;
//	REAL kZero=w*sqrt(M_UZERO*M_EZERO);
  
  int pOrder = 1; //ordem polinomial de aproximacao
  int dim = 2;
  int xDiv = 2;
  int zDiv = 1;
  
  const int meshType = createRectangular;
  timer.start();
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
	CreateGMesh(gmesh, meshType, hDomain, wDomain,  xDiv, zDiv);

  TPZCompMesh *cmesh = CMesh(gmesh, pOrder, urSubs , erSubs , L , theta , lambda , e0 , scale); //funcao para criar a malha computacional
	bool optimizeBandwidth = false;
  TPZAnalysis an(cmesh,optimizeBandwidth);

  //configuracoes do objeto de analise
   //CUIDADO
  TPZVec<long> skyVec;
  cmesh->Skyline(skyVec);
  TPZSkylineNSymStructMatrix skylstr(cmesh);
  skylstr.SetNumThreads(0);
  an.SetStructuralMatrix(skylstr);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU); //caso simetrico
  an.SetSolver(step);
  
  TPZStack<std::string> scalnames, vecnames;
  vecnames.Push("absE");//setando para imprimir campoeletrico
	//vecnames.Push("solAnal");//setando para imprimir campoeletrico
  std::string plotfile= "../HCurl2DEField.vtk";//arquivo de saida que estara na pasta debug
  an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
  int postProcessResolution = 2 ;//define resolucao do pos processamento
  //fim das configuracoes do objeto de analise
	
  //  const REAL k0 = 2*M_PI*freq*sqrt(M_UZERO*M_EZERO);
  //  const REAL L = wDomain ;
  {
    std::ofstream file("../cmesh.txt");
    cmesh->Print(file);
  }
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
  stiff.Print("KPZglobalNed = " , std::cout , EMathematicaInput);
  std::cout<<"saindo do assemble"<<std::endl;
  std::cout<<"entrando no solver"<<std::endl;
  an.Solve();
  std::cout<<"saindo do solver"<<std::endl;
  const TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
  std::string name("../../sol");
  name.append( std::to_string(meshType) );
  name.append(".csv");
  ExportCSV(solucao,name.c_str());
  
  
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

void FilterOthersFunctions(TPZManVector<long> &active, TPZCompMesh * fcmesh){
  active.Resize(0, 0);
  int ncon = fcmesh->NConnects();
  
  // DOF related with the Q system
  for(int i = 0; i < ncon; i++)
  {
    TPZConnect &con = fcmesh->ConnectVec()[i];
    int seqnum = con.SequenceNumber();
    int pos = fcmesh->Block().Position(seqnum);
    int blocksize = fcmesh->Block().Size(seqnum);// Just for quads
    
    int vs = active.size();
    active.Resize(vs+blocksize);
    for(int ieq = 0; ieq<blocksize; ieq++)
    {
      active[vs+ieq] = pos+ieq;
    }
  }
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
      std::ofstream outTxt("../gmeshRectangular.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshRectangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
    }
      break;
    case createTriangular:
    {
			REAL rot = 0.0;
			gengrid = new TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
			gengrid->SetElementType(ETriangle);
      std::ofstream outTxt("../gmeshTriangular.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshTriangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
      
    }
      break;
    case createZigZag:
    {
			REAL rot = 0.0;
			gengrid = new TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
			gengrid->SetElementType(EQuadrilateral);
			gengrid->SetZigZagPattern();
      std::ofstream outTxt("../gmeshZigZag.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshZigZag.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
    }
      break;
    default:
      DebugStop();
      break;
  }
	gmesh = new TPZGeoMesh();
	const int matId = 1; //define id para um material(formulacao fraca)
	const int bc0 = -1; //define id para um material(cond contorno dirichlet)
	const int bc1 = -2; //define id para um material(cond contorno mista)
	gengrid->Read(gmesh , matId);
	
	
	gengrid->SetBC(gmesh, ulCoord, llCoord, bc0);
	gengrid->SetBC(gmesh, urCoord, ulCoord, bc0);
	gengrid->SetBC(gmesh, lrCoord, urCoord, bc1);
	gengrid->SetBC(gmesh, llCoord, lrCoord, bc0);
	
	//gmesh->ResetConnectivities();
	//converte malha xy para xz
	for (int i = 0; i < gmesh->NNodes() ; i ++)
	{
		TPZGeoNode node = gmesh->NodeVec()[i];
		node.SetCoord( 2, node.Coord(1) );
		node.SetCoord( 1, 0. );
		gmesh->NodeVec()[i] = node;
	}
	
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

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL L, REAL theta, REAL lambda, REAL e0, REAL scale)
{
  
  const int dim = 2; //dimensao do problema
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno dirichlet)
  const int bc1 = -2; //define id para um material(cond contorno mista)
  enum{ dirichlet = 0, neumann, mixed}; //tipo da condicao de contorno do problema
  // Criando material
  TPZMatHCurl2D *material = new TPZMatHCurl2D(matId , lambda , ur , er , e0 , theta, scale);//criando material que implementa a
  //formulacao fraca do problema de validacao
  
  ///set forcing function for L2 projection of desired solution
//  TPZDummyFunction<STATE> *Ltracer = new TPZDummyFunction<STATE>(solAnal);
//  TPZAutoPointer<TPZFunction<STATE> > fLTracer = Ltracer;
//  material->SetForcingFunctionExact(fLTracer);
  
  
  
  
  ///criar malha computacional
  
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmesh->SetDimModel(dim);//seta dimensao do modelo
  // Inserindo material na malha
  cmesh->InsertMaterialObject(material);
		
  ///electrical conductor boundary conditions
  TPZFNMatrix<1,STATE> val1(1,1,0.), val2(1,1,0.);
  
  TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
  
  
  ///mixed boundary conditions for open wall
	REAL w=2.*M_PI*M_C/lambda;
	REAL k0=w*sqrt(M_UZERO*M_EZERO);
	
	val1(0,0) = -1.*imaginary*k0*cos(theta);
  val2(0,0) = -2.*imaginary*k0*cos(theta)*e0*exp(imaginary*k0*L*cos(theta));
	
  TPZMaterial * BCond1 = material->CreateBC(material, bc1, mixed, val1, val2);//cria material que implementa a condicao de contorno mista
  
  cmesh->InsertMaterialObject(BCond0);//insere material na malha
  cmesh->InsertMaterialObject(BCond1);//insere material na malha
  
  cmesh->SetAllCreateFunctionsHDiv();//define espaco de aproximacao
  
  //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
  cmesh->AutoBuild();
  TPZAdmChunkVector< TPZCompEl* > elVec = cmesh->ElementVec();
  
  for (int i = 0; i < cmesh->NElements(); i++) {
    TPZCompElHDiv < pzshape::TPZShapeQuad > *el = dynamic_cast<TPZCompElHDiv <pzshape::TPZShapeQuad > *>( elVec[i] );
    if ( el == NULL) {
      continue;
    }
    el->SetSideOrient(4,  1);
    el->SetSideOrient(5,  1);
    el->SetSideOrient(6, -1);
    el->SetSideOrient(7, -1);
  }
  
  for (int i = 0; i < cmesh->NElements(); i++) {
    TPZCompElHDiv < pzshape::TPZShapeTriang > *el = dynamic_cast<TPZCompElHDiv <pzshape::TPZShapeTriang > *>( elVec[i] );
    if ( el == NULL) {
      continue;
    }
    el->SetSideOrient(3,  1);
    el->SetSideOrient(4, -1);
    el->SetSideOrient(5,  1);
  }
  
  
  if (pOrder == 1) {
    //cmesh->CleanUpUnconnectedNodes();
    TPZCreateApproximationSpace::MakeRaviartThomas(*cmesh);
    cmesh->CleanUpUnconnectedNodes();
  }
  else{//for now only lowest order elements are avaliable
    DebugStop();
  }
  
  //cmesh->AutoBuild();
  return cmesh;
}

inline STATE urSubs(const TPZVec<REAL> &x)
{
	//return 1.;
	return ( 2.-imaginary*0.1 );
}

inline STATE erSubs(const TPZVec<REAL> &x)
{
	//return 4.26 - imaginary * 0.6;
	return ( 4.+(2.-imaginary*0.1)*(1.-x[0]/(5*1550*1e-9))*(1.-x[0]/(5*1550*1e-9)) );
}
void ExportCSV(const TPZFMatrix<STATE> matriz, const char* name)
{
  std::ofstream file;
  if(name!=NULL) file.open(name);
  
  for (int i = 0 ; i < matriz.Rows(); i++) {
    for(int j = 0 ; j < matriz.Cols() ; j++){
      file<<matriz.GetVal(i, j);
      if (j != matriz.Cols() - 1) {
        file<<" , ";
      }
    }
    file<<std::endl;
  }
  file.close();
}

