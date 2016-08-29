
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <math.h>

// Geometry
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "tpzhierarquicalgrid.h"

// Computational mesh
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzlog.h"

// Materials
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TPZPoroPermCoupling.h"

// Analysis
#include "pzanalysis.h"
#include "TPZPoroPermAnalysis.h"

// Matrix
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"

// Simulation data structure
#include "TPZSimulationData.h"

// Methods declarations


// Rectangular geometry
TPZGeoMesh * RockBox(TPZVec<REAL> dx_dy, TPZVec<int> n);
void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);

void UniformRefinement(TPZGeoMesh *gmesh, int nh);
void UniformRefinement(TPZGeoMesh * gmesh, int nh, int mat_id);


// Create a computational mesh for deformation
TPZCompMesh * CMesh_Deformation(TPZGeoMesh * gmesh, int order);

// Create a computational mesh for pore pressure excess
TPZCompMesh * CMesh_PorePressure(TPZGeoMesh * gmesh, int order);

// Create a computational mesh for pore pressure excess
TPZCompMesh * CMesh_PorePermeabilityCoupling(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data);


#ifdef LOG4CXX
static LoggerPtr log_data(Logger::getLogger("pz.permeabilityc"));
#endif

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.0e14;
    
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + "/Projects/PermeabilityCoupling/";
    FileName += "permeability_c_log.cfg";
    InitializePZLOG(FileName);
#endif
    
#ifdef LOG4CXX
    
    if(log_data->isInfoEnabled())
    {
        std::stringstream sout;
        sout << " Defining the case ... " << std::endl;
        LOGPZ_DEBUG(log_data,sout.str())
    }
#endif
    
    TPZSimulationData * sim_data = new TPZSimulationData;
    
    REAL dt = 0.1;
    int n_steps = 10;
    REAL epsilon_res = 0.01;
    REAL epsilon_corr = 0.0001;
    int n_corrections = 2;
    
    /** @brief Definition gravity field */
    TPZVec<REAL> g(2,0.0);
    g[1] = -0.0*9.81;
    
    sim_data->SetGravity(g);
    sim_data->SetTimeControls(n_steps, dt);
    sim_data->SetNumericControls(n_corrections, epsilon_res, epsilon_corr);
    
#ifdef LOG4CXX
    
    if(log_data->isInfoEnabled())
    {
        std::stringstream sout;
        sout << " Computing Geometry ... " << std::endl;
        LOGPZ_DEBUG(log_data,sout.str())
    }
#endif
    
    TPZVec<REAL> dx_dy(2);
    TPZVec<int> n(2);
    
    dx_dy[0] = 0.5; // x - direction
    dx_dy[1] = 1.0; // y - direction
    n[0] = 10; // x - direction
    n[1] = 10; // y - direction
    
    TPZGeoMesh * gmesh = RockBox(dx_dy,n);
    
#ifdef LOG4CXX
    
    if(log_data->isInfoEnabled())
    {
        std::stringstream sout;
        sout << " Computing Geometry accomplished... " << std::endl;
        LOGPZ_DEBUG(log_data,sout.str())
    }
#endif
    
    // Create the approximation space
    int deformation_order = 1;
    int pore_pressure_order = 1;
    
    // Create multiphysisc mesh
    TPZManVector<TPZCompMesh * , 2 > mesh_vector(2);
    
    mesh_vector[0] = CMesh_Deformation(gmesh, deformation_order);
    mesh_vector[1] = CMesh_PorePressure(gmesh, pore_pressure_order);
    
    TPZCompMesh * cmesh_poro_perm_coupling = CMesh_PorePermeabilityCoupling(gmesh, mesh_vector,sim_data);

    // Create the static analysis
    
    // Run Static analysis
    // @omar:: the initial condition is set up to zero for displacement and pore pressure excess
    
    // Create the Transient analysis
    
    bool mustOptimizeBandwidth = true;
    TPZPoroPermAnalysis * time_analysis = new TPZPoroPermAnalysis;
    time_analysis->SetCompMesh(cmesh_poro_perm_coupling,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
    TPZSkylineNSymStructMatrix skyl(cmesh_poro_perm_coupling);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(skyl);
    
    // Run Transient analysis
    time_analysis->Run_Evolution();
    
    std::cout << " Execution finished" << std::endl;
	return EXIT_SUCCESS;
}

void UniformRefinement(TPZGeoMesh *gmesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        long n = gmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (sons);
        }//for i
    }//ref
}

void UniformRefinement(TPZGeoMesh * gmesh, int nh, int mat_id)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        long n = gmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1){
                if (gel->MaterialId()== mat_id){
                    gel->Divide (sons);
                }
            }
        }//for i
    }//ref
}

TPZCompMesh * CMesh_PorePermeabilityCoupling(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
    // Plane strain assumption
    int planestress = 0;
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    REAL l = 6.009e6;
    REAL mu = 3.486e6;
    REAL l_u = 6.009e6;
    
    // Aproximation Space of order -> pOrder
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    TPZPoroPermCoupling * material = new TPZPoroPermCoupling(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetPlaneProblem(planestress);
    material->SetParameters(l, mu, l_u);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet_x_y   = 0;
    int dirichlet_x     = 1;
    int dirichlet_y     = 2;
    int neumann_y     = 5;
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    

    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet_x_y, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet_x, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = -1000.0;
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, neumann_y, val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet_x, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    // Setando Multifisico
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, cmesh);
    
#ifdef PZDEBUG
    std::ofstream out("CMeshMultiPhysics.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

// Create a computational mesh for deformation;
TPZCompMesh * CMesh_PorePressure(TPZGeoMesh * gmesh, int order){
    
    // Plane strain assumption
    //    int planestress = 0;
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    // Aproximation Space of order -> pOrder
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    // Creating a material object
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet = 0;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshPorePressure.txt");
    cmesh->Print(out);
#endif

    return cmesh;
}

// Create a computational mesh for deformation;
TPZCompMesh * CMesh_Deformation(TPZGeoMesh * gmesh, int order){
        
    // Plane strain assumption
//    int planestress = 0;
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    // Aproximation Space of order -> pOrder
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);


    // Creating a material object
    int nstate = 2;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet = 0;

    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshDeformation.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}



TPZGeoMesh * RockBox(TPZVec<REAL> dx_dy, TPZVec<int> n){
    
    REAL t=0.0,dx,dy;
    int n_elements;
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
    GeoMesh1->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh1->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    int matid=1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
    GeoMesh1->BuildConnectivity();
    GeoMesh1->SetDimension(0);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew1.txt");
        GeoMesh1->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh1,Dummyfile, true);
    }
    
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh1);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc = new TPZDummyFunction<STATE>(ParametricfunctionX);
    CreateGridFrom.SetParametricFunction(ParFunc);
    CreateGridFrom.SetFrontBackMatId(bc_left,bc_right);
    dx = dx_dy[0];
    n_elements = n[0];
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh2 = CreateGridFrom.ComputeExtrusion(t, dx, n_elements);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew2.txt");
        GeoMesh2->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew2.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh2,Dummyfile, true);
    }
    
    
    
    TPZHierarquicalGrid CreateGridFrom2(GeoMesh2);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc2 = new TPZDummyFunction<STATE>(ParametricfunctionY);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    CreateGridFrom2.SetFrontBackMatId(bc_bottom,bc_top);
    dy = dx_dy[1];
    n_elements = n[1];
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dy, n_elements);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMeshNew3.txt");
        GeoMesh3->Print(argument);
        std::ofstream Dummyfile("GeometricMeshNew3.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh3,Dummyfile, true);
    }
    
    long last_node = GeoMesh3->NNodes() - 1;
    long last_element = GeoMesh3->NElements() - 1;
    long node_id = GeoMesh3->NodeVec()[last_node].Id();
    long element_id = GeoMesh3->Element(last_element)->Id();
    const std::string name("Geomechanic Reservoir box");
    GeoMesh3->SetName(name);
    GeoMesh3->SetMaxNodeId(node_id);
    GeoMesh3->SetMaxElementId(element_id);
    
    return GeoMesh3;
    
}

void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0]; // x
    X[1] = 0.0; // y
    X[2] = 0.0; // z
}

void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0; // x
    X[1] = par[0]; // y
    X[2] = 0.0; // z
}

