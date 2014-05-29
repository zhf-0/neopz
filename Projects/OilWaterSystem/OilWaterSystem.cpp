#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzmatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzconvectionproblem.h"
#include "pzmultiphase.h"
#include "pzl2projection.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzlog.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"

#include "pzequationfilter.h"
#include "pzgradientreconstruction.h"
#include "pzl2projection.h"
#include "pzbfilestream.h"

#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzpoisson3d.h"
#include <time.h>
#include <stdio.h>

// Using Log4cXX as logging tool
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.multiphase.data"));
#endif
//
// End Using Log4cXX as logging tool

TPZCompMesh *ComputationalMeshElasticity(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshPseudopressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshBulkflux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshGravitationalflux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshWaterSaturation(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshMultiphase(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
TPZCompMesh *L2ProjectionQ(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
TPZCompMesh *L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
TPZCompMesh *L2ProjectionS(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle);

void SolveSyst(TPZAnalysis &an, TPZCompMesh *cmesh);
void InitialFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void InitialSaturation(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void PosProcessBulkflux(TPZAnalysis &an, std::string plotfile);
void PosProcessL2(TPZAnalysis &an, std::string plotfile);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void UniformRefinement(TPZGeoMesh *gMesh, int nh);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void SolveSystemTransient(REAL deltaT,REAL maxTime, TPZAnalysis *NonLinearAn, TPZAnalysis *NonLinearAnTan, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);
void GetElSolution(TPZCompEl * cel, TPZCompMesh * mphysics);
void CheckConvergence(TPZFMatrix<STATE> &RUattn,TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);
void ComputeResidual(TPZFMatrix<STATE> &RUattn,REAL &alpha, TPZFMatrix<STATE> &DeltaU, TPZFMatrix<STATE> &ResAlpha, TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);
void CheckElConvergence(TPZFMatrix<STATE> &RUattn,TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);
void GetElSolution(TPZCompEl * cel, TPZCompMesh * mphysics);
void BulkFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &Val);
void CleanGradientSolution(TPZFMatrix<STATE> &Solution, TPZManVector<long> &Gradients);

void FilterPressureFluxEquation(TPZMultiphase *mymaterial, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics, TPZAnalysis &an);
void FilterHigherOrderSaturations(TPZManVector<long> &active,TPZManVector<long> &nonactive, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);

bool ftriang = false;
REAL angle = 1.0*M_PI/16.0;

int main()
{
	
#ifdef LOG4CXX
	InitializePZLOG("OilWaterLog4cxx.cfg");
#endif		
	
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif		
	
	//  	gRefDBase.InitializeAllUniformRefPatterns();
	
	//	Reading mesh
	
	std::string GridFileName, dirname = PZSOURCEDIR;
	GridFileName = dirname + "/Projects/OilWaterSystem/";
	//GridFileName += "OilWaterSystemUnit.dump";
    GridFileName += "Labyrinth.dump";    

//	GridFileName = "Labyrinth.dump";	
	//	GridFileName = "OilWaterSystemUnitTwo.dump";
	//	GridFileName = "OilWaterSystemUnitOneHRef.dump";	
	//	GridFileName = "OilWaterSystemUnitTwoHRef.dump";
	
	
	bool twoMaterial = false;
	TPZReadGIDGrid GeometryInfo;
	GeometryInfo.SetfDimensionlessL(1.0);
	TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
//	RotateGeomesh(gmesh, angle);
	{
		//	Print Geometrical Base Mesh
		std::ofstream argument("GeometicMesh.txt");
		gmesh->Print(argument);
		std::ofstream Dummyfile("GeometricMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
	}
	
	// 	std::vector<REAL> dd(2,0);
	
	
	int Href = 2;
	int div = 0;		
	int POrderBulkFlux = 1;
	int POrderGravitationalFlux = 1;	
	int POrderPseudopressure = 1;
	int POrderWaterSaturation = 1;	
	
	UniformRefinement(gmesh, Href);	
	
	{
		//	Print Geometrical Base Mesh
		std::ofstream argument("RefGeometicMesh.txt");
		gmesh->Print(argument);
		std::ofstream Dummyfile("RefGeometricMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
		//		gmesh->BuildConnectivity();
	}	
	
	
	// Computational meshes	
	
	//	First computational mesh
	TPZCompMesh * CMeshBulkflux = ComputationalMeshBulkflux(gmesh, POrderBulkFlux);
	//	Print Second computational mesh
	std::ofstream ArgumentBulkflux("ComputationalMeshForBulkflux.txt");
	CMeshBulkflux->Print(ArgumentBulkflux);	
	
	//	Auxiliar computational mesh
	TPZCompMesh * CMeshGravitationalflux = ComputationalMeshGravitationalflux(gmesh, POrderGravitationalFlux);
	//	Print Second computational mesh
	std::ofstream ArgumentGravitationalflux("ComputationalMeshFoflux.txt");
	CMeshGravitationalflux->Print(ArgumentGravitationalflux);		
	
	//	Second computational mesh
	TPZCompMesh * CMeshPseudopressure = ComputationalMeshPseudopressure(gmesh, POrderPseudopressure);
	//	Print First computational mesh
	std::ofstream ArgumentPseudopressure("ComputationalMeshForPseudopressure.txt");
	CMeshPseudopressure->Print(ArgumentPseudopressure);
	
	//	Third computational mesh
	TPZCompMesh * CMeshWaterSaturation = ComputationalMeshWaterSaturation(gmesh, POrderWaterSaturation);
	//	Print Second computational mesh
	std::ofstream ArgumentWaterSaturation("ComputationalMeshForWaterSaturation.txt");
	CMeshWaterSaturation->Print(ArgumentWaterSaturation);
	
	TPZAnalysis Anbulkflux(CMeshBulkflux);
    std::string outputfile1;
	outputfile1 = "SolutionBulkflux";
    std::stringstream outputfiletemp1;
    outputfiletemp1 << outputfile1 << ".vtk";
    std::string plotfilebuklflux = outputfiletemp1.str();
	TPZFMatrix<STATE> InitialQSolution = Anbulkflux.Solution();	
 	int rwosQ= InitialQSolution.Rows();
//  	for (int i=0; i < InitialQSolution.Rows(); i++) 
//  	{
//  		InitialQSolution(i)=0.0;
//  	}
	
//	InitialQSolution(0) =  0.1;
//	InitialQSolution(1) =  0.1;	
//	
//	InitialQSolution(6) = -0.1;
//	InitialQSolution(7) = -0.1;	
	
	Anbulkflux.LoadSolution(InitialQSolution);
	int num= InitialQSolution.Rows();
	
	
	TPZVec<STATE> soliniQ(num,0.0);
    TPZCompMesh  * cmeshQL2 = L2ProjectionQ(gmesh, POrderBulkFlux, soliniQ);
    
    TPZAnalysis anQL2(cmeshQL2);
//    SolveSyst(anQL2, cmeshQL2);
    Anbulkflux.LoadSolution(InitialQSolution);	
	PosProcessBulkflux(Anbulkflux,plotfilebuklflux);
	
	
	
	TPZAnalysis AnPressure(CMeshPseudopressure);
	
    std::string outputfile2;
	outputfile2 = "SolutionPressure";
    std::stringstream outputfiletemp2;
    outputfiletemp2 << outputfile2 << ".vtk";
    std::string plotfilePressure = outputfiletemp2.str();
	TPZFMatrix<STATE> InitialPSolution = AnPressure.Solution();	
// 	for (int i=0; i < InitialPSolution.Rows(); i++) {
// 		InitialPSolution(i)=0.0;
// 	}
	
	
	TPZVec<STATE> solini(InitialPSolution.Rows(),0.0);
    TPZCompMesh  * cmeshL2 = L2ProjectionP(gmesh, POrderPseudopressure, solini);
    
    TPZAnalysis anL2(cmeshL2);
//    SolveSyst(anL2, cmeshL2);
	
    AnPressure.LoadSolution(anL2.Solution());
	
    PosProcessL2(AnPressure,plotfilePressure);	
	
	
	TPZAnalysis AnSaturation(CMeshWaterSaturation);
    std::string outputfile3;
    outputfile3 = "SolutionSaturation";
    std::stringstream outputfiletemp3;
    outputfiletemp3 << outputfile3 << ".vtk";
    std::string plotfileSaturation = outputfiletemp3.str();
	TPZFMatrix<STATE> InitialSSolution = AnSaturation.Solution();	
// 	for (int i=0; i < InitialSSolution.Rows(); i++) 
//     {
// 		InitialSSolution(i)=0.0;
// 	}
	
	TPZVec<STATE> solSini(InitialSSolution.Rows(),0.0);
    TPZCompMesh  * cmeshsatL2 = L2ProjectionS(gmesh, POrderWaterSaturation, solSini);
    
    TPZAnalysis ansatL2(cmeshsatL2);
//    SolveSyst(ansatL2, cmeshsatL2);
	
	AnSaturation.LoadSolution(InitialSSolution);
	PosProcessL2(AnSaturation,plotfileSaturation);
	
	// 	//	This is so rare!!
	// 	if (twoMaterial) 
	// 	{
	// 		gmesh->AddInterfaceMaterial(1,2, 1);
	// 		gmesh->AddInterfaceMaterial(2,1, 1);
	// 	}  

    //	Multiphysics Mesh
    TPZVec<TPZCompMesh *> meshvec(3);//4);
    meshvec[0] = CMeshBulkflux;
    meshvec[1] = CMeshPseudopressure;	
    meshvec[2] = CMeshWaterSaturation;
//    meshvec[3] = CMeshGravitationalflux;	
	
    
    TPZCompMesh * MultiphysicsMesh = ComputationalMeshMultiphase(gmesh,meshvec);
    std::ofstream ArgumentMultiphysic("MultiphysicsMesh.txt");
    MultiphysicsMesh->Print(ArgumentMultiphysic);
	
	
    TPZAnalysis *MultiphysicsAn = new TPZAnalysis(MultiphysicsMesh);
    TPZAnalysis *MultiphysicsAnTan = new TPZAnalysis(MultiphysicsMesh);	
    int	Nthreads = 4;
    TPZSkylineNSymStructMatrix matsk(MultiphysicsMesh);
    TPZSkylineNSymStructMatrix matskTan(MultiphysicsMesh);	

	//TPZFrontStructMatrix <TPZFrontNonSym<STATE> > matsk(MultiphysicsMesh);
    
	MultiphysicsAn->SetStructuralMatrix(matsk);
    MultiphysicsAn->StructMatrix()->SetNumThreads(Nthreads);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU); 					
    MultiphysicsAn->SetSolver(step);

	// I can't believe that this is needed!!!
	MultiphysicsAnTan->SetStructuralMatrix(matskTan);
    MultiphysicsAnTan->StructMatrix()->SetNumThreads(Nthreads);
    TPZStepSolver<STATE> steptan;
    steptan.SetDirect(ELU); 					
    MultiphysicsAnTan->SetSolver(steptan);	
    TPZManVector<long> AllConnects(0),NoGradients(0),WithGradients(0);
    FilterHigherOrderSaturations(NoGradients,WithGradients,meshvec,MultiphysicsMesh);
    TPZFMatrix<STATE> SolutiontoMod=MultiphysicsMesh->Solution();
    CleanGradientSolution(SolutiontoMod,WithGradients);
    MultiphysicsAn->LoadSolution(SolutiontoMod);
    std::string outputfile;
    outputfile = "TransientSolutionini";
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    PosProcessMultphysics(meshvec,MultiphysicsMesh,*MultiphysicsAn,plotfile);		

	//     Tima control parameters
    
    REAL hour = 60.0*60.0;
    REAL day = 24.0*hour;
    REAL year = 365.0*day;
	
    REAL deltaT = 0.5;
    REAL maxTime = 10.0;
    SolveSystemTransient(deltaT, maxTime, MultiphysicsAn, MultiphysicsAnTan, meshvec, MultiphysicsMesh);
	
    return 0;
	
}

TPZCompMesh * ComputationalMeshElasticity(TPZGeoMesh *gmesh, int pOrder)
{
//	// Plane strain assumption
//	int planestress = 0;
//	
//	// Getting mesh dimension
//	int dim = 2;
//	
//	TPZCompEl::SetgOrder(pOrder);
//	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
//	cmesh->SetDimModel(dim);
//	
//	
//	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
//	TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
//	TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);		
//	//		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
//	TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
//	//		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,0, val1, val2);
//	TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
//	//		TPZMaterial * BCond8 = material2->CreateBC(mat2,8,0, val1, val2);		
//	
//	cmesh->SetAllCreateFunctionsContinuous(); // H1 approximation space
//	cmesh->InsertMaterialObject(BCond2);
//	cmesh->InsertMaterialObject(BCond3);		
//	//		cmesh->InsertMaterialObject(BCond4);
//	cmesh->InsertMaterialObject(BCond4);
//	//		cmesh->InsertMaterialObject(BCond6);
//	cmesh->InsertMaterialObject(BCond5);
//	//		cmesh->InsertMaterialObject(BCond8);	
//	
//	cmesh->AutoBuild();
//	
//	
//	return cmesh;
	
}

TPZCompMesh *ComputationalMeshBulkflux(TPZGeoMesh *gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	int matId1 = 1;
	int matId2 = 2;		
	
	TPZMatPoisson3d *material1;
	material1 = new TPZMatPoisson3d(matId1,dim);
	TPZMaterial * mat1(material1);
	material1->NStateVariables();
	
	//		TPZMatPoisson3d *material2;
	//		material2 = new TPZMatPoisson3d(matId2,dim);
	//		TPZMaterial * mat2(material2);
	//		material2->NStateVariables();		
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(mat1);
	//		cmesh->InsertMaterialObject(mat2);		
	
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
	TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);		
	TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
	TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
	TPZMaterial * BCond6 = material1->CreateBC(mat1,6,0, val1, val2);
	TPZMaterial * BCond7 = material1->CreateBC(mat1,7,0, val1, val2);
//	TPZMaterial * BCond8 = material2->CreateBC(mat2,8,0, val1, val2);		
	
	cmesh->SetAllCreateFunctionsHDiv();
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);		
	cmesh->InsertMaterialObject(BCond4);
	cmesh->InsertMaterialObject(BCond5);
	cmesh->InsertMaterialObject(BCond6);
	cmesh->InsertMaterialObject(BCond7);
	//		cmesh->InsertMaterialObject(BCond8);
	
	
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;
}

TPZCompMesh *ComputationalMeshGravitationalflux(TPZGeoMesh *gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	int matId1 = 1;
	int matId2 = 2;		
	
	TPZMatPoisson3d *material1;
	material1 = new TPZMatPoisson3d(matId1,dim);
	TPZMaterial * mat1(material1);
	material1->NStateVariables();
	
	//		TPZMatPoisson3d *material2;
	//		material2 = new TPZMatPoisson3d(matId2,dim);
	//		TPZMaterial * mat2(material2);
	//		material2->NStateVariables();		
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(mat1);
	//		cmesh->InsertMaterialObject(mat2);		
	
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
	TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);		
	TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
	TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
	TPZMaterial * BCond6 = material1->CreateBC(mat1,6,0, val1, val2);
	TPZMaterial * BCond7 = material1->CreateBC(mat1,7,0, val1, val2);
	//	TPZMaterial * BCond8 = material2->CreateBC(mat2,8,0, val1, val2);		
	
	cmesh->SetAllCreateFunctionsHDiv();
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);		
	cmesh->InsertMaterialObject(BCond4);
	cmesh->InsertMaterialObject(BCond5);
	cmesh->InsertMaterialObject(BCond6);
	cmesh->InsertMaterialObject(BCond7);
	//		cmesh->InsertMaterialObject(BCond8);
	
	
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;
}

TPZCompMesh *ComputationalMeshPseudopressure(TPZGeoMesh *gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	int matId1 = 1;
	int matId2 = 2;
	
	TPZMatPoisson3d *material1;
	material1 = new TPZMatPoisson3d(matId1,dim);
	TPZMaterial * mat1(material1);
	material1->NStateVariables();
	
	//		TPZMatPoisson3d *material2;
	//		material2 = new TPZMatPoisson3d(matId2,dim);
	//		TPZMaterial * mat2(material2);
	//		material2->NStateVariables();
	
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(mat1);
	//		cmesh->InsertMaterialObject(mat2);		
	
	
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
	TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);
	TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
	TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
	TPZMaterial * BCond6 = material1->CreateBC(mat1,6,0, val1, val2);
	TPZMaterial * BCond7 = material1->CreateBC(mat1,7,0, val1, val2);	
//	TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);		
	//		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
//	TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
	//		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,0, val1, val2);
//	TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
	//		TPZMaterial * BCond8 = material2->CreateBC(mat2,8,0, val1, val2);		
	
	cmesh->SetAllCreateFunctionsDiscontinuous(); // L2 approximation space
	
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);		
	cmesh->InsertMaterialObject(BCond4);
	cmesh->InsertMaterialObject(BCond5);
	cmesh->InsertMaterialObject(BCond6);
	cmesh->InsertMaterialObject(BCond7);	
	
//	cmesh->InsertMaterialObject(BCond2);
//	cmesh->InsertMaterialObject(BCond3);		
//	//		cmesh->InsertMaterialObject(BCond4);
//	cmesh->InsertMaterialObject(BCond4);
//	//		cmesh->InsertMaterialObject(BCond6);
//	cmesh->InsertMaterialObject(BCond5);
//	//		cmesh->InsertMaterialObject(BCond8);
	
	
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	cmesh->AutoBuild();
	
	
	///inserir connect da pressao
	int ncon = cmesh->NConnects();
	for(int i=0; i<ncon; i++)
	{
		TPZConnect &newnod = cmesh->ConnectVec()[i];
		//newnod.SetPressure(true);
		newnod.SetLagrangeMultiplier(1);
	}
	
	///set order total da shape
	int nel = cmesh->NElements();
	for(int i=0; i<nel; i++){
		TPZCompEl *cel = cmesh->ElementVec()[i];
		TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
		celdisc->SetConstC(1.);
		celdisc->SetCenterPoint(0, 0.);
		celdisc->SetCenterPoint(1, 0.);
		celdisc->SetCenterPoint(2, 0.);
		celdisc->SetTrueUseQsiEta();
		if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
		{
			if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
			else celdisc->SetTensorialShape();
		}
	}
	//		
	//		
	//#ifdef DEBUG
	//		int ncel = cmesh->NElements();
	//		for(int i =0; i<ncel; i++){
	//			TPZCompEl * compEl = cmesh->ElementVec()[i];
	//			if(!compEl) continue;
	//			TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
	//			if(facel)DebugStop();
	//			
	//		}
	//#endif		
	
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	
	
	return cmesh;
}


TPZCompMesh *ComputationalMeshWaterSaturation(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	int matId1 = 1;
	int matId2 = 2;	
	
	TPZMatConvectionProblem *material1 = new TPZMatConvectionProblem(matId1,dim);
	TPZMaterial * mat1(material1);
	
	//		TPZMatConvectionProblem *material2 = new TPZMatConvectionProblem(matId2,dim);
	//		TPZMaterial * mat2(material2);		
	
	TPZVec<REAL> convdir(dim,0.);
	convdir[0]=1.;
	REAL flux = 0.;
	REAL rho = 1.;
	
	material1->SetParameters(rho,convdir);
	material1->SetInternalFlux(flux);
	material1->NStateVariables();
	
	//		material2->SetParameters(rho,convdir);
	//		material2->SetInternalFlux(flux);
	//		material2->NStateVariables();		
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->InsertMaterialObject(mat1);
	//		cmesh->InsertMaterialObject(mat2);		
	
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
	TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);
	TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
	TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
	TPZMaterial * BCond6 = material1->CreateBC(mat1,6,0, val1, val2);	
	TPZMaterial * BCond7 = material1->CreateBC(mat1,7,0, val1, val2);
	
//	TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);		
//	//		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
//	TPZMaterial * BCond4 = material1->CreateBC(mat1,4,0, val1, val2);
//	//		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,0, val1, val2);
//	TPZMaterial * BCond5 = material1->CreateBC(mat1,5,0, val1, val2);
//	//		TPZMaterial * BCond8 = material2->CreateBC(mat2,8,0, val1, val2);
	
	cmesh->SetAllCreateFunctionsDiscontinuous();//  L2 approximation space
	
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);		
	cmesh->InsertMaterialObject(BCond4);
	cmesh->InsertMaterialObject(BCond5);
	cmesh->InsertMaterialObject(BCond6);
	cmesh->InsertMaterialObject(BCond7);
    
    // Void material
    int matIdL2Proj = 2;
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,material1->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);     
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
	
//	cmesh->InsertMaterialObject(BCond2);
//	cmesh->InsertMaterialObject(BCond3);		
//	//		cmesh->InsertMaterialObject(BCond4);
//	cmesh->InsertMaterialObject(BCond4);
//	//		cmesh->InsertMaterialObject(BCond6);
//	cmesh->InsertMaterialObject(BCond5);
//	//		cmesh->InsertMaterialObject(BCond8);
	
	///set order total da shape
	int nel = cmesh->NElements();
	for(int i=0; i<nel; i++){
		TPZCompEl *cel = cmesh->ElementVec()[i];
		TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        celdisc->SetConstC(1.);
        celdisc->SetCenterPoint(0, 0.);
        celdisc->SetCenterPoint(1, 0.);
        celdisc->SetCenterPoint(2, 0.);
        celdisc->SetTrueUseQsiEta();
		if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
		{
			if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
			else celdisc->SetTensorialShape();
		}
	}
	
			
	
	
	return cmesh;
}

TPZCompMesh *ComputationalMeshMultiphase(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
	
	
	//Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
	//		mphysics->SetAllCreateFunctionsMultiphysicElem();
	
	bool Dimensionless = true;
	
	int dim =2;
	int matId1 = 1;
	int matId2 = 2;
	// Setting data
	
	TPZMultiphase *material1 = new TPZMultiphase(matId1,dim);
	//		TPZMultiphase *material2 = new TPZMultiphase(matId2,dim);
	
	REAL deltaT = 0.1;
	REAL maxTime = 0.1;
	REAL MPa = 1.0e+6;
	
	if (Dimensionless) 
	{	
		REAL Rhoref,Etaref,Lref,Kref,Pref,Qref;
		Rhoref = 1000.0;
		Etaref = 1.0e-3;
		Lref = 1.0;
		Kref = 1.0e-13;
		Pref = 20.0*MPa;
		Qref = (Rhoref*(Kref/Etaref))*(Pref/Lref);
		
		material1->SetRhoSCreference(Rhoref);
		material1->SetEtaSCreference(Etaref);
		material1->SetLreference(Lref);
		material1->SetKreference(Kref);
		material1->SetPreference(Pref);
	}
	
	std::string GridFileName, dirname = PZSOURCEDIR;
	GridFileName = dirname + "/Projects/OilWaterSystem/";
    GridFileName += "LabyrinthKvaluesH.txt";    
	
	material1->SetTimeStep(1.0);
	material1->fnewWS=true; 
	material1->SetTScheme(1.0, 1.0);
//	material1->LoadKMap("Permeabilities.txt");
	material1->LoadKMap(GridFileName);	
	material1->SetYorN(false);		
	
	TPZMaterial *mat1(material1);
	mphysics->InsertMaterialObject(mat1);
	mphysics->SetDimModel(dim);
	
	
	TPZFMatrix<STATE> val1(4,2,0.), val2(4,1,0.);
	
	val2(0,0)=0.000*cos(angle);// qx
	val2(1,0)=0.000*sin(angle);// qy
	val2(2,0)=1.0*20.0*MPa;// P
	val2(3,0)=1.0;// S	
	TPZMaterial * BCond5 = material1->CreateBC(mat1,6,1, val1, val2);
	
//	val2(0,0)=0.0;// qx
//	val2(1,0)=0.0;// qy
//	val2(2,0)=0.0;// P
//	val2(3,0)=0.0;// S			
//	TPZMaterial * BCond3 = material1->CreateBC(mat1,2,4, val1, val2);		
	//		TPZMaterial * BCond4 = material2->CreateBC(mat2,4,0, val1, val2);
	
	val2(0,0)=0.0;// qx
	val2(1,0)=0.0;// qy
	val2(2,0)=0.0*MPa;// P
	val2(3,0)=0.0;// S			
//	TPZMaterial * BCond2 = material1->CreateBC(mat1,4,4, val1, val2);		
	TPZMaterial * BCond2Nflux = material1->CreateBC(mat1,2,3, val1, val2);
	TPZMaterial * BCond2Nflux2 = material1->CreateBC(mat1,3,3, val1, val2);
	val2(0,0)=0.0;// qx
	val2(1,0)=0.0;// qy
	val2(2,0)=0.0*MPa;// P
	val2(3,0)=0.0;// S		
	TPZMaterial * BCond2Nflux3 = material1->CreateBC(mat1,4,3, val1, val2);
	TPZMaterial * BCond2Nflux4 = material1->CreateBC(mat1,5,3, val1, val2);	
	//		TPZMaterial * BCond4 = material1->CreateBC(mat1,4,2, val1, val2);		
	//		TPZMaterial * BCond6 = material2->CreateBC(mat2,6,3, val1, val2);
	
	val2(0,0)=0.000;// qx
	val2(1,0)=0.0;// qy
	val2(2,0)=1.0*18.0*MPa;// P
	val2(3,0)=0.0;// S			
	TPZMaterial * BCond4 = material1->CreateBC(mat1,7,2, val1, val2);				
		
	
	mphysics->SetAllCreateFunctionsMultiphysicElem();		
//	mphysics->InsertMaterialObject(BCond2);
	mphysics->InsertMaterialObject(BCond2Nflux);
	mphysics->InsertMaterialObject(BCond2Nflux2);
	mphysics->InsertMaterialObject(BCond2Nflux3);
	mphysics->InsertMaterialObject(BCond2Nflux4);	
	//		mphysics->InsertMaterialObject(BCond4);
//	mphysics->InsertMaterialObject(BCond3);
	//		mphysics->InsertMaterialObject(BCond6);
	mphysics->InsertMaterialObject(BCond4);
	//		mphysics->InsertMaterialObject(BCond8);
	mphysics->InsertMaterialObject(BCond5);		
	
	
	mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
	
	// Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);	
	
	mphysics->Reference()->ResetReference();
	mphysics->LoadReferences();		
	
	// Creation of interface elements
	int nel = mphysics->ElementVec().NElements();
	for(int el = 0; el < nel; el++)
	{
		TPZCompEl * compEl = mphysics->ElementVec()[el];
		if(!compEl) continue;
		int index = compEl ->Index();
		if(compEl->Dimension() == mphysics->Dimension())
		{
			TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(mphysics->ElementVec()[index]);
			if(!InterpEl) continue;
			InterpEl->CreateInterfaces();
		}
	}
	
	return mphysics;
}

void PosProcessBulkflux(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(0), vecnames(1);
	vecnames[0]= "FluxL2";
    
	const int dim = 2;
	int div = 5;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malhaflux.txt");
	an.Print("nothing",out);
}

void PosProcessL2(TPZAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(1), vecnames(0);
	scalnames[0]= "Solution";
    
	const int dim = 2;
	int div = 2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malhaflux.txt");
	an.Print("nothing",out);
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
	TPZManVector<std::string,10> scalnames(3), vecnames(1);
	
    scalnames[0] = "WeightedPressure";
    scalnames[1] = "WaterSaturation";
    scalnames[2] = "OilSaturation";
    vecnames[0] = "BulkVelocity";
    
	const int dim = 2;
	int div =2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
    
}

void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		long n = gMesh->NElements();
		for ( long i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
		}//for i
	}//ref
}

void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
{
	
	TPZVec<long > subindex;
	for (long iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		long nel = elvec.NElements(); 
		for(long el=0; el < nel; el++){
			TPZCompEl * compEl = elvec[el];
			if(!compEl) continue;
			long ind = compEl->Index();
			compEl->Divide(ind, subindex, 0);
		}
	}
	
}

void SolveSystemTransient(REAL deltaT,REAL maxTime, TPZAnalysis *NonLinearAn, TPZAnalysis *NonLinearAnTan, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics){
    
	TPZFMatrix<STATE> SolutiontoLoad;
	TPZFMatrix<STATE> SolutiontoSave = meshvec[2]->Solution();
	
//	{
//		TPZBFileStream load;
//		load.OpenRead("MultiphaseSaturationSol.bin");
//		SolutiontoLoad.Read(load,0);
//		meshvec[2]->LoadSolution(SolutiontoLoad);
//		TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);		
//	}


	TPZGradientReconstruction *gradreconst = new TPZGradientReconstruction(false,1.);	
	std::string OutPutFile = "TransientSolutionGRActiveHref2";
    TPZMaterial *mat1 = mphysics->FindMaterial(1);
	//    TPZMaterial *mat2 = mphysics->FindMaterial(2);	
	
    TPZMultiphase * material1 = dynamic_cast<TPZMultiphase *>(mat1);  
	//    TPZMultiphase * material2 = dynamic_cast<TPZMultiphase *>(mat2);  	
    material1->SetTimeStep(deltaT);
    material1->SetTScheme(1.0,1.0);
	bool UsingGradient = true;
	int matIdL2Proj = 2;	
	
	
	//	Starting Newton Iterations
	TPZFMatrix<STATE> DeltaX = mphysics->Solution();
	TPZFMatrix<STATE> Uatn = mphysics->Solution();
	TPZFMatrix<STATE> Uatk = mphysics->Solution();		
	
	
	REAL TimeValue = 0.0;
	REAL Tolerance = 1.0e-7;
	int cent = 0;
	int MaxIterations = 50;
	TimeValue = cent*deltaT;
	REAL NormValue =1.0;
	REAL NormValueT =1.0;
	bool StopCriteria = false;
	TPZFMatrix<STATE> RhsAtn, RhsAtnPlusOne, Residual;
	TPZFMatrix<STATE> RhsAtnT, RhsAtnPlusOneT, ResidualT;	
	
	if (UsingGradient) 
	{	
		meshvec[2]->Reference()->ResetReference();
		meshvec[2]->LoadReferences();		
		TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
		gradreconst->ProjectionL2GradientReconstructed(meshvec[2], matIdL2Proj);
		TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
		NonLinearAn->LoadSolution(mphysics->Solution());	
	}
	
	std::string outputfile;
	outputfile = OutPutFile;
	std::stringstream outputfiletemp;
	outputfiletemp << outputfile << ".vtk";
	std::string plotfile = outputfiletemp.str();
	PosProcessMultphysics(meshvec,mphysics,*NonLinearAn,plotfile);		
	
	
	TPZManVector<long> AllConnects(0),NoGradients(0),WithGradients(0);
	FilterHigherOrderSaturations(NoGradients,WithGradients,meshvec,mphysics);
	AllConnects = NoGradients;	
	AllConnects.Resize(NoGradients.size()+WithGradients.size());
	for (int i=0; i<(WithGradients.size()); i++) 
	{
		AllConnects[i+NoGradients.size()]=WithGradients[i];
	};
	
	NonLinearAn->StructMatrix()->EquationFilter().Reset();
	NonLinearAn->StructMatrix()->EquationFilter().SetActiveEquations(NoGradients);
	
	
	std::cout << " Starting the time computations. " << std::endl;	
	while (TimeValue < maxTime)
	{
		//RotateGeomesh(mphysics->Reference(),angle);
		
		material1->SetLastState();
		NonLinearAn->AssembleResidual();
		RhsAtn = NonLinearAn->Rhs();
		
		material1->SetCurrentState();
		NonLinearAn->Assemble();
		RhsAtnPlusOne = NonLinearAn->Rhs();
		Residual= RhsAtn + RhsAtnPlusOne;		
		NormValue = Norm(Residual);
		
//		material1->SetLastState();
//		NonLinearAnTan->Assemble();
//		RhsAtnT = NonLinearAnTan->Rhs();	
//		
//		material1->SetCurrentState();
//		NonLinearAnTan->Assemble();
//		RhsAtnPlusOneT = NonLinearAnTan->Rhs();
//		ResidualT= RhsAtnT + RhsAtnPlusOneT;		
//		NormValueT = Norm(ResidualT);

		
		int iterations= 0;		
		while (NormValue > Tolerance)
		{		
			
			Residual*=-1.0;
			NonLinearAn->Rhs()=Residual;
			NonLinearAn->Solve();			
			DeltaX = NonLinearAn->Solution();
			Uatk = (Uatn + DeltaX);
			
			
			mphysics->LoadSolution(Uatn + DeltaX);			
			TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
			
#ifdef LOG4CXX
			if(logdata->isDebugEnabled())
			{
				std::stringstream sout;
				sout.precision(20);
                Residual.Print(sout);
				Uatk.Print(sout);		
				LOGPZ_DEBUG(logdata,sout.str());
			}
#endif			
			
//#ifdef LOG4CXX
//			if(logdata->isDebugEnabled())
//			{
//				std::stringstream sout;
//				Uatk.Print("Uatk = ",sout,EMathematicaInput);			
//				LOGPZ_DEBUG(logdata,sout.str());
//			}
//#endif				


			
			if (UsingGradient) 
			{
				gradreconst->ProjectionL2GradientReconstructed(meshvec[2], matIdL2Proj);
				TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
                Uatk=mphysics->Solution();
                CleanGradientSolution(Uatk,WithGradients);
				NonLinearAn->LoadSolution(Uatk);
			}
		
			
//			material1->SetCurrentState();
//			NonLinearAnTan->Assemble();
//			RhsAtnPlusOneT = NonLinearAnTan->Rhs();			
//			ResidualT= RhsAtnT + RhsAtnPlusOneT;
//			NormValueT = Norm(ResidualT);			
			
			material1->SetCurrentState();
			NonLinearAn->Assemble();
			RhsAtnPlusOne = NonLinearAn->Rhs();
			Residual= RhsAtn + RhsAtnPlusOne;
			NormValue = Norm(Residual);	
				
			
#ifdef LOG4CXX
			if(logdata->isDebugEnabled())
			{
				std::stringstream sout;
				sout.precision(15);				
				Uatk.Print(sout);
				Residual.Print("Res = ",sout,EMathematicaInput);
				LOGPZ_DEBUG(logdata,sout.str());
			}
#endif		

			
//#ifdef LOG4CXX
//			if(logdata->isDebugEnabled())
//			{
//				std::stringstream sout;
//				sout.precision(20);
//				Residual.Print(sout);
//				NonLinearAn->Solution().Print(sout);		
//				LOGPZ_DEBUG(logdata,sout.str());
//			}
//#endif			
			
//			NonLinearAn->StructMatrix()->EquationFilter().Reset();
//			NonLinearAn->StructMatrix()->EquationFilter().SetActiveEquations(NoGradients);
//			
//			numofequactive = NonLinearAn->StructMatrix()->EquationFilter().NActiveEquations();			
			
			iterations++;
			std::cout << " Newton's Iteration = : " << iterations  << "     L2 norm = : " << NormValue <<  std::endl;
			if (iterations == MaxIterations) 
			{
				StopCriteria = true;
				std::cout << " Time Step number = : " << iterations  << "\n Exceed max iterations numbers = : " << MaxIterations <<  std::endl;					
				break;
			}

				
			Uatn = Uatk;
			
		}	
			
		
        if (UsingGradient) 
        {
			TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
			gradreconst->ProjectionL2GradientReconstructed(meshvec[2], matIdL2Proj);
			TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
            Uatk=mphysics->Solution();
//             CleanGradientSolution(Uatk,WithGradients);
            NonLinearAn->LoadSolution(Uatk);		
		}	
			
		TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);		
		SolutiontoSave = mphysics->Solution();
		TPZBFileStream save;
        save.OpenWrite("GRSolution.bin");
		SolutiontoSave.Write(save,0);
		
		
		outputfile = OutPutFile;
		std::stringstream outputfiletemp;
		outputfiletemp << outputfile << ".vtk";
		std::string plotfile = outputfiletemp.str();
		PosProcessMultphysics(meshvec,mphysics,*NonLinearAn,plotfile);		
		
		if (StopCriteria) {
			std::cout << " Newton's Iteration = : " << iterations  << "     L2 norm = : " << NormValue <<  std::endl;		
			break;
		}
		
		cent++;
		TimeValue = cent*deltaT;
		
		std::cout << " Time Step :  " << cent  << "  Time :  " << TimeValue <<  std::endl;	
        
/*        if(cent==15)
        {
            UsingGradient=true;
            std::cout << " Now Using Gradient Reconstruction. " << std::endl;                    
        }  */      
		
	}
	
//	CheckConvergence(RhsAtn,NonLinearAn, meshvec, mphysics);
//	CheckElConvergence(RhsAtn,NonLinearAn, meshvec, mphysics);
	
	
}


// Setting up initial conditions

void SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
	
	TPZSkylineStructMatrix full(Cmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//	//Saida de Dados: solucao e  grafico no VT
	//	ofstream file("Solutout");
	//	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

TPZCompMesh *L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialPressure);
    material->SetForcingFunction(forcef);
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	cmesh->AutoBuild();
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }    
    
	return cmesh;
	
}

TPZCompMesh *L2ProjectionQ(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialFlux);
    material->SetForcingFunction(forcef);
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	cmesh->AutoBuild();
	
	return cmesh;
	
}

TPZCompMesh *L2ProjectionS(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialSaturation);
    material->SetForcingFunction(forcef);
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	cmesh->AutoBuild();
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }    
    
	return cmesh;
	
}

// It requires modfify L2 number os state variables
void InitialFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
	REAL x = pt[0];
	REAL y = pt[1];
    disp[0] = 0.0;
	//    disp[1] = 0.0;	
    
}

void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
	REAL x = pt[0];
	REAL y = pt[1];
    disp[0] = 0.0*(1.0 - 0.1 * x);
    
}

void InitialSaturation(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
	REAL x = pt[0];
	REAL y = pt[1];
	
	if (x<=1.0) {
		    disp[0] = 0.0;
	}	else {
		    disp[0] = 0.0;
	}
    
}


void CheckConvergence(TPZFMatrix<STATE> &RUattn,TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics)
{
	
	TPZFMatrix<REAL> U = NonLinearAn->Solution();
	long neq = mphysics->NEquations();
	
	int nsteps = 10;
	REAL alpha;
	TPZFMatrix<REAL> alphas(nsteps,1,0.0),ResNorm(nsteps,1,0.0),ConvergenceOrder(nsteps-1,1,0.0);
	
	TPZFMatrix<REAL> DeltaX(neq,1,0.00001),ResAlpha(neq,0.0);
	
    for(int i = 0; i < nsteps; i++)
    {
        alpha = (1.0*i+1.0)/10.0;
        alphas(i,0) = log(alpha);
	    ComputeResidual(RUattn,alpha, DeltaX, ResAlpha, NonLinearAn, meshvec, mphysics);	    
	    ResNorm(i,0)=log(Norm(ResAlpha));
		
	    NonLinearAn->LoadSolution(U);
	    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);		
		
	}
	
    for(int i = 1; i < nsteps; i++){ ConvergenceOrder(i-1,0) =  (ResNorm(i,0)-ResNorm(i-1,0))/(alphas(i,0)-alphas(i-1,0));}	
	
	std::ofstream outfile("CheckConvergence.txt");
	
	ConvergenceOrder.Print("CheckConv =",outfile,EMathematicaInput);	
	NonLinearAn->LoadSolution(U);
	TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);	
	
}

void ComputeResidual(TPZFMatrix<STATE> &RUattn, REAL &alpha, TPZFMatrix<STATE> &DeltaU, TPZFMatrix<STATE> &ResAlpha, TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics)
{
	
	TPZFMatrix<STATE> TangentRes;
	TPZMaterial *mat1 = mphysics->FindMaterial(1);
	TPZMultiphase * material1 = dynamic_cast<TPZMultiphase *>(mat1); 	
	
	
	// Computing the first part of the residual expresion.		  
	
	material1->SetCurrentState();
	NonLinearAn->Assemble();	
	TPZFMatrix<STATE> RhsAtnPlusOne = NonLinearAn->Rhs();
	
	TPZFMatrix<STATE> ResidualAtU = RUattn + RhsAtnPlusOne;		
	NonLinearAn->Solver().Matrix()->Multiply((1.0)*alpha*DeltaU,TangentRes);
	
	NonLinearAn->LoadSolution(NonLinearAn->Solution()+alpha*DeltaU);			
	TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);			
	
	material1->SetCurrentState();
	NonLinearAn->Assemble();
	RhsAtnPlusOne = NonLinearAn->Rhs();
	TPZFMatrix<STATE> ResidualAtUplusAlphaDeltaU = RUattn + RhsAtnPlusOne;
	
	
	ResAlpha = ((1.0)*ResidualAtUplusAlphaDeltaU - ((1.0)*ResidualAtU + (1.0) * TangentRes));
    
	
}

void CheckElConvergence(TPZFMatrix<STATE> &RUattn,TPZAnalysis *NonLinearAn, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics)
{
	
	TPZMaterial *mat = mphysics->FindMaterial(1);
	TPZMultiphase * material = dynamic_cast<TPZMultiphase *>(mat);
	
	TPZFMatrix<STATE> Usol = NonLinearAn->Solution();
	std::ofstream outsol("Solution.txt");
	Usol.Print("Sol =",outsol,EMathematicaInput);
	outsol.flush();
	
	int NumberofEl = mphysics->ElementVec().NElements();      
	long neq = mphysics->NEquations();
	TPZElementMatrix elk(mphysics, TPZElementMatrix::EK),elf(mphysics, TPZElementMatrix::EF);      
	
	int nsteps = 9;
	STATE du=0.0001;
	
	
	TPZFMatrix<REAL> DeltaU(neq,1,du);
	TPZFNMatrix<4,REAL> alphas(nsteps,1,0.0),ElConvergenceOrder(nsteps-1,1,0.0);
	TPZFNMatrix<9,REAL> res(nsteps,1,0.0);
	
	std::ofstream outfile("CheckConvergencebyElements.txt");
	
	for(long i = 0; i < NumberofEl; i++ )
	{
		
		NonLinearAn->LoadSolution(RUattn);			
		TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);	
	    
		TPZCompEl * iel = mphysics->ElementVec()[i];  
		material->SetLastState();
		
		iel->Print(outfile);
		
		iel->CalcStiff(elk,elf);
		TPZFNMatrix<9,REAL> elResidualUn = elf.fMat;	
		
		NonLinearAn->LoadSolution(Usol);			
		TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
		
		material->SetCurrentState();	  
		iel->CalcStiff(elk,elf);	  
		TPZFNMatrix<9,REAL> elTangentU = elk.fMat;
		TPZFNMatrix<9,REAL> elResidualU = elf.fMat;
		
		int SizeOfElMat = elTangentU.Rows();
		TPZFNMatrix<9,STATE> deltaUel(SizeOfElMat,1,du),ResTangUel(SizeOfElMat,1,0.0);
		TPZFNMatrix<9,STATE> ResU(SizeOfElMat,1,0.0),ResUalphadu(SizeOfElMat,1,0.0);
		
		ResU= elResidualU+elResidualUn;
		
		std::ofstream outek("TangentEl.txt");
		std::ofstream outef("ResidualEl.txt");
		std::ofstream outefn("ResidualEln.txt");	  
		elTangentU.Print("Tangent = ",outek,EMathematicaInput);
		outek.flush();
		
		elResidualU.Print("ResU = ",outef,EMathematicaInput);	  
		outef.flush();
		
		elResidualUn.Print("Resn = ",outefn,EMathematicaInput);	  
		outefn.flush();
		
		REAL alpha = 0;
		for(int j = 0; j < nsteps; j++)
		{	
			
			alpha = (1.0*j+1.0)/10.0;
			elTangentU.Multiply(alpha*deltaUel,ResTangUel);		
			alphas(j,0) = log(alpha);
			
			NonLinearAn->LoadSolution(Usol+alpha*DeltaU);			
			TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
			
			iel->CalcStiff(elk,elf);
			std::ofstream outefa("ResUa.txt");
			TPZFNMatrix<9,REAL> elResidualUalphadx = elf.fMat;
			elResidualUalphadx.Print("ResUa =",outefa,EMathematicaInput);
			ResUalphadu = elResidualUalphadx + elResidualUn;
			STATE NormValue = Norm(ResUalphadu-(ResU+ResTangUel));
			res(j) = log(NormValue);
			
		}      
		
		for(int j = 1; j < nsteps ; j++){ElConvergenceOrder(j-1,0)=(res(j,0)-res(j-1,0))/(alphas(j,0)-alphas(j-1,0));}  
		// 	  AllOrders[i]= ElConvergenceOrder;
		ElConvergenceOrder.Print("CheckConv = ",outfile,EMathematicaInput);
		outfile.flush();
		
	}	  
	
	
	
}

void GetElSolution(TPZCompEl * cel, TPZCompMesh * mphysics)
{
	if(!cel) {return;}
	
	TPZBlock<STATE> &Block = mphysics->Block(); 
	int NumberOfEquations = cel->NEquations();  
	int NumberOfConnects = cel->NConnects();
	TPZFMatrix<STATE> elSolution(NumberOfEquations,1,0.0);
	long DestinationIndex = 0L;
	
	for(int iconnect = 0; iconnect < NumberOfConnects; iconnect++)
	{
		TPZConnect Connect = cel->Connect(iconnect);
		int seq = Connect.SequenceNumber();
		int SizeOfBlockAtseq = Block.Size(seq);
		int BlockGlobalPosition = Block.Position(seq);
		
		for(int iblock = 0; iblock	 < SizeOfBlockAtseq; iblock++)
		{
			
			elSolution(DestinationIndex++,0) = mphysics->Solution()[BlockGlobalPosition+iblock];  
			
		}
	}
	
}




void FilterHigherOrderSaturations(TPZManVector<long> &active, TPZManVector<long> &nonactive, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics)
{
    int ncon_saturation = meshvec[2]->NConnects();
    int ncon = mphysics->NConnects();
    for(int i = 0; i < ncon-ncon_saturation; i++)
    {
        TPZConnect &con = mphysics->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = mphysics->Block().Position(seqnum);
        int blocksize = mphysics->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+blocksize);
        for(int ieq = 0; ieq<blocksize; ieq++){
            active[vs+ieq] = pos+ieq;
        }
    }	
	
    for(int i = ncon-ncon_saturation; i<ncon; i++)
    {
        TPZConnect &con = mphysics->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = mphysics->Block().Position(seqnum);
        int blocksize = mphysics->Block().Size(seqnum);
        int vs = active.size();
        active.Resize(vs+1);
		
        int ieq = blocksize-1;
        active[vs] = pos+ieq;
    }
	
    for(int i = ncon-ncon_saturation; i<ncon; i++)
    {
        TPZConnect &con = mphysics->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = mphysics->Block().Position(seqnum);
        int blocksize = mphysics->Block().Size(seqnum);
        int vs = nonactive.size();
        nonactive.Resize(vs+blocksize-1);
        for(int ieq = 0; ieq<blocksize-1; ieq++)
		{
            nonactive[vs+ieq] = pos+ieq;
        }		
		
    }
	
}

void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle)
{
	REAL theta = CounterClockwiseAngle; 
	// It represents a 3D rotation around the z axis.
	TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
	RotationMatrix(0,0) =	+cos(theta);
	RotationMatrix(0,1) =	-sin(theta);
	RotationMatrix(1,0) =	+sin(theta);
	RotationMatrix(1,1) =	+cos(theta);
	RotationMatrix(2,2) = 1.0;
	TPZVec<STATE> iCoords(3,0.0);
	TPZVec<STATE> iCoordsRotated(3,0.0);
	
	RotationMatrix.Print("Rotation = ");
	
	int NumberofGeoNodes = gmesh->NNodes();
	for (int inode = 0; inode < NumberofGeoNodes; inode++) 
	{
		TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
		GeoNode.GetCoordinates(iCoords);
		// Apply rotation
		iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
		iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
		iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
//		std::cout << "iCoords[0] " << iCoords[0] << std::endl;
//		std::cout << "iCoords[1] " << iCoords[1] << std::endl;
//		std::cout << "iCoords[2] " << iCoords[2] << std::endl;			
//		std::cout << "iCoordsRotated[0] " << iCoordsRotated[0] << std::endl;
//		std::cout << "iCoordsRotated[1] " << iCoordsRotated[1] << std::endl;
//		std::cout << "iCoordsRotated[2] " << iCoordsRotated[2] << std::endl;		
		GeoNode.SetCoord(iCoordsRotated);
		gmesh->NodeVec()[inode] = GeoNode;
	}
	
}

void CleanGradientSolution(TPZFMatrix<STATE> &Solution, TPZManVector<long> &Gradients)
{
    for(int i=0; i < Gradients.size(); i++ )
    {
        Solution(Gradients[i],0)=0.0;
    }   
}