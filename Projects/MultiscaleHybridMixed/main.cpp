    #ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZCompElMHM.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"
#include "pzelasmat.h" 

#include "pzbuildmultiphysicsmesh.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"

#include "GlobalMHM.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.MHM"));
#endif

using namespace std;

const int matIdLoc = 1;
const int matIdGlob = 1;  // must be different from matInterno.
const int matIdLocGlob = 1;

const int dirichlet = 0;
const int neumann = 1;
const int mixed = 2;

const int bc1 = -1;
const int bc2 = -2;
const int bc3 = -3;
const int bc4 = -4;

TPZGeoMesh* MalhaGeomMHM(int NRefUnif, REAL Lx, REAL Ly);
TPZCompMesh* MalhaCompLoc(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh* MalhaCompGlob(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh* MalhaCompLocGlob(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh*>& meshvec);
void RemoveInterfaces( TPZCompMesh* cmesh );
void SolveMHM( TPZAnalysis &an );

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	InitializePZLOG();
#endif

	int p = 1;
	int NRefUnif = 1;
    REAL Lx = 1.;
    REAL Ly = 1.;

    // ==============
    // Geometric mesh
    // ==============

    TPZGeoMesh* gmesh = MalhaGeomMHM(NRefUnif,Lx,Ly);
	ofstream gm("gmesh.txt");
	gmesh->Print(gm);

    // ========================
    // Local computational mesh: must be the 1st in the meshVec
    // ========================

	TPZCompMesh* cmeshLoc = MalhaCompLoc(gmesh, p);
	ofstream cmLoc("cmeshLoc.txt");
	cmeshLoc->Print(cmLoc);

    // =========================
    // Global computational mesh: must be the 2nd in the meshVec
    // =========================

    TPZCompMesh* cmeshGlob = MalhaCompGlob(gmesh, 0);
	ofstream cmGlob("cmeshGlob.txt");
	cmeshGlob->Print(cmGlob);

    // =================================
    // Local + Global computational mesh
    // =================================

    TPZVec< TPZCompMesh* > meshvec(2);
    meshvec[0] = cmeshLoc;
    meshvec[1] = cmeshGlob;

    TPZCompMesh* cmeshLocGlob = MalhaCompLocGlob(gmesh, meshvec);
    std::ofstream cmLocGlob("cmeshLocGlob.txt");
    cmeshLocGlob->Print( cmLocGlob );

    // ==============================
    // Solution of the global problem
    // ==============================

    TPZAnalysis an(cmeshLocGlob);

    SolveMHM( an );

	return EXIT_SUCCESS;
}

TPZGeoMesh *MalhaGeomMHM(int NRefUnif, REAL Lx, REAL Ly)
{
    int Qnodes = 4;

	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);

	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolLine(2);

	//indice dos nos
	int id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}

	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}

	//indice dos elementos
	id = 0;

    //elementos internos
    TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 2;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id, TopolQuad, matIdLoc, *gmesh);
	id++;

    //elementos de contorno
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id, TopolLine, bc1, *gmesh);
	id++;

	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id, TopolLine, bc2, *gmesh);
	id++;

	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id, TopolLine, bc3, *gmesh);
	id++;

	TopolLine[0] = 3;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id, TopolLine, bc4, *gmesh);
	id++;

    // Define entre quais materiais vou criar interfaces. Terceiro argumento é o tipo de material que quero nessa interface.
    // exemplo: gmesh->AddInterfaceMaterial(lagrangemat, matInterno, matInterno);
    //gmesh->AddInterfaceMaterial(matIdGlob, matIdLoc, matIdLoc);
    //gmesh->AddInterfaceMaterial(matIdGlob, bc1, bc1);
    //gmesh->AddInterfaceMaterial(matIdGlob, bc2, bc2);
    //gmesh->AddInterfaceMaterial(matIdGlob, bc3, bc3);
    //gmesh->AddInterfaceMaterial(matIdGlob, bc4, bc4);

    //construir a malha
	gmesh->BuildConnectivity();

    //Refinamento uniforme
	for( int ref = 0; ref < NRefUnif; ref++ )
    {
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ )
        {
			TPZGeoEl * gel = gmesh->ElementVec()[i];
            gel->Divide (filhos);
		}
	}

	return gmesh;

}

TPZCompMesh* MalhaCompLoc(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
    TPZMatPoisson3d *matLoc = new TPZMatPoisson3d( matIdLoc, 2);
	TPZAutoPointer<TPZMaterial> mat1(matLoc);
	matLoc->NStateVariables();

	TPZCompMesh * cmeshLoc = new TPZCompMesh(gmesh);
    cmeshLoc->SetDefaultOrder(pOrder);
	cmeshLoc->SetDimModel(dim);

	cmeshLoc->InsertMaterialObject(mat1);

	///Inserir condicao de contorno
	TPZFMatrix val1(2,2,0.), val2(2,1,0.);
	REAL uN=0.;
	val2(0,0)=uN;
	TPZMaterial * BCondN1 = matLoc->CreateBC(mat1, bc1, neumann, val1, val2);
	cmeshLoc->InsertMaterialObject(BCondN1);

    TPZMaterial * BCondN2 = matLoc->CreateBC(mat1, bc3, neumann, val1, val2);
    cmeshLoc->InsertMaterialObject(BCondN2);

	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	REAL uD=0.;
	val22(0,0)=uD;
	TPZMaterial * BCondD1 = matLoc->CreateBC(mat1, bc2, dirichlet, val12, val22);
	cmeshLoc->InsertMaterialObject(BCondD1);

	TPZMaterial * BCondD2 = matLoc->CreateBC(mat1, bc4, dirichlet, val12, val22);
	cmeshLoc->InsertMaterialObject(BCondD2);

    cmeshLoc->SetAllCreateFunctionsMHM();

    cmeshLoc->AutoBuild();
    cmeshLoc->AdjustBoundaryElements();
    cmeshLoc->CleanUpUnconnectedNodes();

	return cmeshLoc;
}

TPZCompMesh* MalhaCompGlob(TPZGeoMesh * gmesh,int pOrder)
{
    gmesh->ResetReference();
    TPZCompMesh* cmeshGlob = new TPZCompMesh(gmesh);

    cmeshGlob->SetDefaultOrder(pOrder);
    cmeshGlob->SetDimModel(2);
    cmeshGlob->SetAllCreateFunctionsDiscontinuous();

    TPZMaterial* matGlob = new TPZMatPoisson3d( matIdGlob, 2);  // fake...
    TPZAutoPointer<TPZMaterial> mat1(matGlob);
    cmeshGlob->InsertMaterialObject( matGlob );

	///Inserir condicao de contorno
	TPZFMatrix val1(2,2,0.), val2(2,1,0.);
	REAL uN=0.;
	val2(0,0)=uN;
	TPZMaterial * BCondN1 = matGlob->CreateBC(mat1, bc1, neumann, val1, val2);
	cmeshGlob->InsertMaterialObject(BCondN1);
    
    TPZMaterial * BCondN2 = matGlob->CreateBC(mat1, bc3, neumann, val1, val2);
    cmeshGlob->InsertMaterialObject(BCondN2);
    
	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	REAL uD=0.;
	val22(0,0)=uD;
	TPZMaterial * BCondD1 = matGlob->CreateBC(mat1, bc2, dirichlet, val12, val22);
	cmeshGlob->InsertMaterialObject(BCondD1);
    
	TPZMaterial * BCondD2 = matGlob->CreateBC(mat1, bc4, dirichlet, val12, val22);
	cmeshGlob->InsertMaterialObject(BCondD2);

    cmeshGlob->AutoBuild();
    RemoveInterfaces( cmeshGlob );
    cmeshGlob->AdjustBoundaryElements();
    cmeshGlob->CleanUpUnconnectedNodes();

    return cmeshGlob;
}

TPZCompMesh* MalhaCompLocGlob(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh*>& meshvec)
{
    gmesh->ResetReference();
    TPZCompMesh* mphysics = new TPZCompMesh(gmesh);

    mphysics->SetDimModel(2);
    mphysics->SetAllCreateFunctionsMultiphysicElemMHM();

    GlobalMHM *material = new GlobalMHM(matIdLocGlob, 2);
	TPZAutoPointer<TPZMaterial> mat1(material);
	material->NStateVariables();
	mphysics->InsertMaterialObject(mat1);

	///Inserir condicao de contorno
	TPZFMatrix val1(2,2,0.), val2(2,1,0.);
	REAL uN=0.;
	val2(0,0)=uN;
	TPZMaterial * BCondN1 = material->CreateBC(mat1, bc1, neumann, val1, val2);
	mphysics->InsertMaterialObject(BCondN1);

    TPZMaterial * BCondN2 = material->CreateBC(mat1, bc3, neumann, val1, val2);
    mphysics->InsertMaterialObject(BCondN2);

	TPZFMatrix val12(2,2,0.), val22(2,1,0.);
	REAL uD=0.;
	val22(0,0)=uD;
	TPZMaterial * BCondD1 = material->CreateBC(mat1, bc2, dirichlet, val12, val22);
	mphysics->InsertMaterialObject(BCondD1);

	TPZMaterial * BCondD2 = material->CreateBC(mat1, bc4, dirichlet, val12, val22);
	mphysics->InsertMaterialObject(BCondD2);

    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();

    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);

    return mphysics;
}

void SolveMHM ( TPZAnalysis &an )
{
	TPZCompMesh *malha = an.Mesh();

    //TPZBandStructMatrix mat(malha);
	//TPZSkylineStructMatrix mat(malha);// requer decomposição simétrica, não pode ser LU!
	//TPZBlockDiagonalStructMatrix mat(malha);//ok!
	//TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
	TPZFStructMatrix mat( malha );// ok! matriz estrutural cheia
	//TPZSpStructMatrix mat( malha );//matriz estrutural esparsa (???? NÃO FUNCIONOU !!!!!!!!!!)
	TPZStepSolver solv;
	solv.SetDirect (  ELU );//ECholesky);// ELU , ELDLt ,

    //	cout << "ELDLt " << endl;
	an.SetSolver ( solv );
	an.SetStructuralMatrix ( mat );
	cout << endl;
	an.Solution().Redim ( 0,0 );
	cout << "Assemble " << endl;
	an.Assemble();
    std::ofstream fileout("RigidezGlobal.txt");
    an.Solver().Matrix()->Print("Rigidez", fileout, EMathematicaInput);

	//an.Solve();

//    TPZVec<std::string> scalar_names(1);
//	TPZVec<std::string> vec_names(0);

//	scalar_names[0] = "Solution";

//	char buf[256] ;
//	sprintf(buf, "solution.vtk");
//	an.DefineGraphMesh(2, scalar_names, vec_names, buf);
//
//	an.PostProcess(0);

}

void RemoveInterfaces( TPZCompMesh* cmesh )
{
    int n = cmesh->ElementVec().NElements();
    for(int i = 0; i < n; i++)
    {
        TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc*>( cmesh->ElementVec()[i] );
        if(!disc)
            continue;
        disc->RemoveInterfaces();
    }
}
