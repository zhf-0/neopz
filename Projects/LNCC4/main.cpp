 /***************************************************************************
 *   Copyright (C) 2009 by joao *
 *   joao@joao-laptop   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h" 
#include "pzburger.h"
#include "pzbndcond.h"


#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "tpzsparseblockdiagonalstructmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzbndmat.h" 
#include "pzfstrmatrix.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "tpzautopointer.h"
#include "math.h"
#include "pzlog.h"
#include <sstream>
#include <fstream>
#include "TPZCopySolve.h"

#include "TPZTimer.h"
#include "TPZRefPattern.h"
#include "pzgeoquad.h"
#include "pzespmat.h"

#include "pztransientanalysis.h"
#include "pztransientmat.h"
#include "pzl2projection.h"


// PARABÓLICO
using namespace std;

/*----------MALHAS-----------*/

TPZGeoMesh * MalhaGeo ( const int h );
TPZCompMeshReferred *CreateCompMesh ( TPZGeoMesh &gmesh, int porder );
TPZCompMeshReferred *CreateCompMeshProjection ( TPZGeoMesh &gmesh, int porder );

/*----------SOLVERS----------*/
void SolveIterative(TPZAnalysis &an);
void SolveTransient ( TPZTransientAnalysis<TPZBurger> &an, TPZAnalysis &anProjection );

/*------------GRAFICOS---------------*/
void SolutionGrafic(TPZAnalysis& an_p1,ofstream& AE, int h, int p, ofstream& OUT);
void SolutionError(TPZAnalysis& an ,ofstream& AE,  int h, int p, ofstream& OUT);
 
//vetor de carga e condições de contorno
static REAL pi=4.*atan(1.);
void ForcingInitialCondition( const TPZVec<REAL> &x, TPZVec<REAL> &disp );
void Forcing1( const TPZVec<REAL> &x, TPZVec<REAL> &disp );
void CC1 ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC2 ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC3 ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC4 ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC1Projection ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC2Projection ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC3Projection ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void CC4Projection ( const TPZVec<REAL> &x, TPZVec<REAL> &f );
void ExactSolution( TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv);



int main() 
{
	/* Dados iniciais*/
	int h=0;
	int p=2; 

	std::ofstream AE ( "results.txt" );
	AE.precision(16);
	
	TPZGeoMesh *gmesh = MalhaGeo(h);
	
	TPZCompMeshReferred *cmeshProjection = CreateCompMeshProjection(*gmesh,p);
	cmeshProjection->SetName ( "cmeshProjection" );
	TPZAnalysis anProjection(cmeshProjection);
	gmesh->ResetReference();
	TPZCompMeshReferred *cmesh = CreateCompMesh(*gmesh,p);
	cmesh->SetName ( "cmesh" );
	bool IsLinear=false;
	TPZTransientAnalysis<TPZBurger> an(cmesh, IsLinear);
	
	SolveTransient ( an, anProjection );	
	
	return EXIT_SUCCESS;
}

TPZGeoMesh * MalhaGeo( const int h ){
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}}; 
	TPZGeoEl *elvec[2];
	int nnode = 4;
	int nelem = 2;
	int nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}
	
	int el;
	
	int nodindAll[2][3]={{0,1,2},{2,3,0}};
	//int nodindAll[1][4]={{0,1,2,3}};
	for ( el=0; el<nelem; el++ )
	{
		TPZVec<int> nodind(4);//(3);
		nodind[0]=nodindAll[el][0];
		nodind[1]=nodindAll[el][1];
		nodind[2]=nodindAll[el][2];
		nodind[3]=nodindAll[el][3];		
		int index;
		elvec[el] = gmesh->CreateGeoElement (ETriangle,nodind,1,index );
		//elvec[el] = gmesh->CreateGeoElement (EQuadrilateral,nodind,1,index );		
	}
//	int nodindAll[1][4]={{0,1,2,3}};
//	for ( el=0; el<nelem; el++ )
//	{
//		TPZVec<int> nodind(4);
//		nodind[0]=nodindAll[el][0];
//		nodind[1]=nodindAll[el][1];
//		nodind[2]=nodindAll[el][2];
//		nodind[3]=nodindAll[el][3];
//		int index;
//		elvec[el] = gmesh->CreateGeoElement (EQuadrilateral,nodind,1,index );
//	}
	
	
	gmesh->BuildConnectivity();
	
//	TPZGeoElBC gbc1 ( elvec[0],4,-1);
//	TPZGeoElBC gbc2 ( elvec[0],5,-2);
//	TPZGeoElBC gbc3 ( elvec[0],6,-3);
//	TPZGeoElBC gbc4 ( elvec[0],7,-4);
	
	TPZGeoElBC gbc1 ( elvec[0],3,-1);
	TPZGeoElBC gbc2 ( elvec[0],4,-2);
	TPZGeoElBC gbc3 ( elvec[1],3,-3);
	TPZGeoElBC gbc4 ( elvec[1],4,-4);
	
//	TPZGeoElBC gbc5 ( elvec[0],0,-5);// um ponto (a origem)
	
	ofstream arcg ( "gmesh.txt" );
	gmesh->Print ( arcg );
	
	///Refinamento uniforme
	for ( int ref = 0; ref < h; ref++ )
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ )
		{
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
		}//for i
	}//ref
	return gmesh;
}

TPZCompMeshReferred *CreateCompMesh ( TPZGeoMesh &gmesh, int porder ){
	TPZCompEl::SetgOrder ( porder );
	TPZCompMeshReferred *result = new TPZCompMeshReferred( &gmesh );
	result->SetDimModel ( 2 );
	//result->SetAllCreateFunctionsDiscontinuous(); 
	result->SetAllCreateFunctionsContinuous();
	
	TPZTransientMaterial< TPZBurger> *material ;
	material = new TPZTransientMaterial< TPZBurger > ( 1 , 2, 1.); // id, dim, deltat

	
	TPZVec<REAL> convdir(3);//direção da convecção
	convdir[0]=0.;
	convdir[1]=0.;
	REAL diff= 1.;
	REAL conv= 0.;
	material-> SetParameters(diff, conv, convdir);
	TPZAutoPointer<TPZMaterial> mat ( material );
	//material->SetNonSymmetric();
	material->SetSymmetric();
	material->SetLinearContext(false);
	
	TPZAutoPointer<TPZFunction> LoadVector = new TPZDummyFunction(Forcing1);
	TPZAutoPointer<TPZFunction> BC1 = new TPZDummyFunction(CC1);	
	TPZAutoPointer<TPZFunction> BC2 = new TPZDummyFunction(CC2);	
	TPZAutoPointer<TPZFunction> BC3 = new TPZDummyFunction(CC3);	
	TPZAutoPointer<TPZFunction> BC4 = new TPZDummyFunction(CC4);		
	material->SetForcingFunction ( LoadVector );
	
	result->InsertMaterialObject ( mat );
	
	TPZFMatrix val1 ( 1,1,0. ), val2 ( 1,1,0. );// 0 é Dirichlet, 1 é Neumann, 2 é Robin(implementada apenas no Contínuo)

	TPZAutoPointer<TPZMaterial> bnd1 = material->CreateBC ( mat,-1,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd2 = material->CreateBC ( mat,-2,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd3 = material->CreateBC ( mat,-3,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd4 = material->CreateBC ( mat,-4,0, val1, val2 );
	
	bnd1->SetForcingFunction ( BC1 );
	bnd2->SetForcingFunction ( BC2 );
	bnd3->SetForcingFunction ( BC3 );
	bnd4->SetForcingFunction ( BC4 ); 
	
	result->InsertMaterialObject ( bnd1 ); 
	result->InsertMaterialObject ( bnd2 );
	result->InsertMaterialObject ( bnd3 );
	result->InsertMaterialObject ( bnd4 );
	
	result->AutoBuild();
	//
	ofstream arc ( "cmesh.txt" );
	arc << "NEquation = "<< result->NEquations() << endl;
	result->Print ( arc );
	ofstream arcg ( "gmesh.txt" );
	result->Reference()->Print ( arcg );
	//
	result->SetName("CMesh1");
	//
	return result; 
}

TPZCompMeshReferred *CreateCompMeshProjection ( TPZGeoMesh &gmesh, int porder ){
	TPZCompEl::SetgOrder ( porder );
	TPZCompMeshReferred *result = new TPZCompMeshReferred( &gmesh );
	result->SetDimModel ( 2 );
	//result->SetAllCreateFunctionsDiscontinuous(); 
	result->SetAllCreateFunctionsContinuous();
	
	TPZL2Projection *material ;
	TPZVec<REAL> Sol;
	material = new TPZL2Projection ( 1,2, 1, Sol );
	
//	TPZL2Projection(int id, int dim, int nstate, TPZVec<REAL> &sol,int IntegrationOrder = -1);
	
	TPZAutoPointer<TPZMaterial> mat ( material );
	
	TPZAutoPointer<TPZFunction> LoadVector = new TPZDummyFunction(ForcingInitialCondition);
	TPZAutoPointer<TPZFunction> BC1 = new TPZDummyFunction(CC1Projection);	
	TPZAutoPointer<TPZFunction> BC2 = new TPZDummyFunction(CC2Projection);	
	TPZAutoPointer<TPZFunction> BC3 = new TPZDummyFunction(CC3Projection);	
	TPZAutoPointer<TPZFunction> BC4 = new TPZDummyFunction(CC4Projection);		
	material->SetForcingFunction ( LoadVector );
	
	result->InsertMaterialObject ( mat );
	
	TPZFMatrix val1 ( 1,1,0. ), val2 ( 1,1,0. );// 0 é Dirichlet, 1 é Neumann, 2 é Robin(implementada apenas no Contínuo)
	
	TPZAutoPointer<TPZMaterial> bnd1 = material->CreateBC ( mat,-1,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd2 = material->CreateBC ( mat,-2,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd3 = material->CreateBC ( mat,-3,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd4 = material->CreateBC ( mat,-4,0, val1, val2 );
	
	bnd1->SetForcingFunction ( BC1 );
	bnd2->SetForcingFunction ( BC2 );
	bnd3->SetForcingFunction ( BC3 );
	bnd4->SetForcingFunction ( BC4 ); 
	
	result->InsertMaterialObject ( bnd1 ); 
	result->InsertMaterialObject ( bnd2 );
	result->InsertMaterialObject ( bnd3 );
	result->InsertMaterialObject ( bnd4 );
	
	result->AutoBuild();
	//
	ofstream arc ( "cmesh.txt" );
	arc << "NEquation = "<< result->NEquations() << endl;
	result->Print ( arc );
	ofstream arcg ( "gmesh.txt" );
	result->Reference()->Print ( arcg );
	//
	result->SetName("CMesh1");
	//
	return result;
}

/**------------------------- SOLVERS -------------------------*/

void SolveTransient ( TPZTransientAnalysis<TPZBurger > &an, TPZAnalysis &anProjection)
{
	std::ofstream AE ( "results.txt" );
	AE.precision(16);
	/* 1 - projetar condição incial */
	TPZCompMesh *ProjectionCMesh = anProjection.Mesh();
	TPZBandStructMatrix matProj(ProjectionCMesh);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix Un;
	
	anProjection.SetStructuralMatrix ( matProj );
	TPZStepSolver solvProjection;
	solvProjection.SetDirect ( ELDLt);//ECholesky);// ELU , ELDLt , 
	cout << "ELU " << endl;
	anProjection.SetSolver ( solvProjection );
	cout << endl;
	anProjection.Solution().Redim ( 0,0 );
	cout << "Assemble " << endl;
	anProjection.Assemble();
	anProjection.Solve();
	Un = anProjection.Solution();	
	anProjection.LoadSolution(Un);
//	SolutionGrafic(anProjection, AE, 1, 1, AE);
	
	
	/* 2 - */
	TPZCompMesh *malha = an.Mesh();
	cout << "No equacoes = " << malha->NEquations() << endl;
	TPZBandStructMatrix mat(malha);
	TPZStepSolver solv;
	solv.SetDirect ( ELDLt);//ECholesky);// ELU , ELDLt , 
	cout << "ELDLt " << endl;
	an.SetSolver ( solv );
	an.SetStructuralMatrix ( mat );
	cout << endl;
	cout << "Run " << endl;
	an.LoadSolution(Un);// carrega solução inicial
	an.TimeStep()=0.01;
	an.SetConvergence( 10, 1.0E-6, true);
	an.SetNewtonConvergence(10, 1.0E-5);
	an.SetSaveFrequency( 1, 2);
	bool FromBegining=true;
	bool linesearch =true;//false;
	
	TPZVec<std::string> scalar_names(1);
	TPZVec<std::string> vec_names(0);
	scalar_names[0] = "Solution";
	char buf[256] ;
	sprintf(buf,"iteracao_linearsearchTrue.vtk"); 
	an.DefineGraphMesh(2, scalar_names, vec_names,buf);
	an.RunTransient(std::cout, FromBegining, linesearch);
}


/** ------------------------- GRÁFICOS -------------------------*/
void SolutionGrafic(TPZAnalysis& an ,ofstream& AE,  int h, int p, ofstream& OUT){

	AE<<"------------------------------------"<<endl;
	AE << "hxh, h="<< h<< ", p="<< p<< endl;
	AE << "DOF = " <<  an.Mesh()->NEquations()    << endl;
	AE<<"------------------------------------"<<endl;
	
	
	/// Saida de Dados primal1
	/// 1. Coeficientes da Solucao Numerica 
	//         ofstream filep1("Solution_primal1.out");
	//         TPZFMatrix toprintp1 = an_p1.Solution();
	//         toprintp1.Print("solution", filep1);
	///2. Saida para dx ou vtk
	TPZVec<std::string> scalar_names(1);
	TPZVec<std::string> vec_names(0);

	scalar_names[0] = "Solution";
//	scalar_names[1] = "dx";
//	scalar_names[2] = "dy";
//	scalar_names[3] = "ExactSolution";
//	scalar_names[4] = "Exactdx";
//	scalar_names[5] = "Exactdy";
	
	
	char buf[256] ;
	sprintf(buf,"iteracao_h%d_p%d.vtk",h,p); 
	an.DefineGraphMesh(2, scalar_names, vec_names,buf);

	an.PostProcess(2);// o parametro é a resolução que será acrescida na visualização.
}

void SolutionError(TPZAnalysis& an ,ofstream& AE,  int h, int p, ofstream& OUT){
///Plotar as normas dos erros 
an.SetExact(ExactSolution);
TPZVec<REAL> pos;
an.PostProcess(pos,AE); // Calcula os erros.(pzanalysis.cpp)
}

void ForcingInitialCondition(const TPZVec<REAL> &x, TPZVec<REAL> &disp )
{
	disp.Resize(1);
	disp[0] = sin(pi* x[0])*sin(pi*x[1]);//u_0
}

void Forcing1(const TPZVec<REAL> &x, TPZVec<REAL> &disp )
{
	disp.Resize(1);
	disp[0] = 0.;//u_0
}
void CC1 ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-1
{
	f[0] = 0.;
}

void CC2 ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-2
{
	f[0] = 0.;
}

void CC3 ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-3
{
	f[0] = 0.;
}

void CC4 ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-4
{
	f[0] = 0.;
}
void CC1Projection ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-1
{
	f[0] =  sin(pi* x[0])*sin(pi*x[1]);//u_0
}
void CC2Projection ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-2
{
	f[0] =  sin(pi* x[0])*sin(pi*x[1]);//u_0
}
void CC3Projection ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-3
{
	f[0] =  sin(pi* x[0])*sin(pi*x[1]);//u_0
}
void CC4Projection ( const TPZVec<REAL> &x, TPZVec<REAL> &f )//-4
{
	f[0] =  sin(pi* x[0])*sin(pi*x[1]);//u_0
}
void ExactSolution( TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv)
{
	u.Resize(1);
	u[0] = sin((7.*pi*(1. + x[0]))/3.)*sin((2.*pi*(1. + pow(x[1],2)))/3.);//u		
}

