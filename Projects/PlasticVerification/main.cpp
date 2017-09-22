
#include <iostream>
#include <cstdlib>

// Utils
#include "pzstack.h"

// Geometry desciption
#include "pzgmesh.h"
#include "TPZGmshReader.h"
// Geometry Utils
#include "TPZVTKGeoMesh.h"

// FE Analysis
#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "pzelastoplasticanalysis.h"

// Global Matrix
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"

// Solver
#include "pzstepsolver.h"

// Materials
#include "pzbndcond.h"
#include "TPZElasticCriteria.h"
#include "pzelastoplastic2D.h"
#include "TPZPlasticStepPV.h"

// PostProcessing Data
#include "pzpostprocanalysis.h"


/**
 Structure wich define the simulation controls
 */
struct SimulationControl {
    
    bool            nonlinear_model_Q;
    int             elemen_type;
    int             n_h_levels;
    int             n_p_levels;
    int             p_order;
    int             int_order;
    int             n_threads;
    std::string     domain_type;
    std::string     execution_summary;
    std::string     dump_folder;
    TPZStack<int>   Omega_ids;
    TPZStack<int>   Gamma_ids;
    
    SimulationControl() : nonlinear_model_Q(false),elemen_type(0), n_h_levels(0), n_p_levels(1), p_order(1), int_order(1), n_threads(0), domain_type(""),execution_summary(""),dump_folder(""),Omega_ids(),Gamma_ids()
    {
        
    }
    
    SimulationControl(const SimulationControl &copy) :
        nonlinear_model_Q(copy.nonlinear_model_Q),
        elemen_type(copy.elemen_type),
        n_h_levels(copy.n_h_levels),
        n_p_levels(copy.n_p_levels),
        p_order(copy.p_order),
        int_order(copy.int_order),
        n_threads(copy.n_threads),
        domain_type(copy.domain_type),
        execution_summary(copy.execution_summary),
        dump_folder(copy.dump_folder),
        Omega_ids(copy.Omega_ids),
        Gamma_ids(copy.Gamma_ids)
    {
        
    }
    
    SimulationControl & operator=(const SimulationControl &copy)
    {
        nonlinear_model_Q = copy.nonlinear_model_Q;
        elemen_type = copy.elemen_type;
        n_h_levels = copy.n_h_levels;
        n_p_levels = copy.n_p_levels;
        p_order     = copy.p_order;
        int_order = copy.int_order;
        n_threads = copy.n_threads;
        domain_type = copy.domain_type;
        execution_summary = copy.execution_summary;
        dump_folder = copy.dump_folder;
        Omega_ids = copy.Omega_ids;
        Gamma_ids = copy.Gamma_ids;
        return *this;
    }
    
};


/**
 Computes an approximation based on SimulationControl object

 @param sim_ctrl control object that contains all the required information for compute the approximation
 */
void ComputeApproximation(SimulationControl & sim_ctrl);


/**
 Create the TPZGeoMesh object based on a Gmsh *.msh file

 @param sim_ctrl control object for geometry functionalities
 @return the TPZGeoMesh object
 */
TPZGeoMesh * GeomtricMesh(SimulationControl  & sim_ctrl);


/**
 Create the TPZCompMesh object based on a TPZGeoMesh and the control object
 
 @param gmesh geometry mesh description
 @param sim_ctrl control object for computational mesh functionalities
 @return the TPZCompMesh object
 */
TPZCompMesh * ComputationalMesh(TPZGeoMesh * gmesh, SimulationControl & sim_ctrl);


/**
 BC type = {dirichlet = 0, neumann = 1}

 @param bc_index boundary condition physical tag
 @return bc type and correpongind boundary quantity
 */
std::pair<int, TPZVec<STATE> > BcType(int bc_index);

/**
 Define the material model being used00

 @param sim_ctrl control object for the material model functionalities
 @return the Material object
 */
TPZMaterial * MaterialModel(int index, SimulationControl  & sim_ctrl);

/**
 Create the FE analysis based on a TPZCompMesh

 @param cmesh computational mesh description
 @param sim_ctrl control object for the FE analysis mesh functionalities
 @return the FE analysis
 */
TPZElastoPlasticAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationControl & sim_ctrl);



/**
 Create a FE postprocessing over the integrations points

 @param cmesh computational mesh description
 @param sim_ctrl control object for the FE analysis mesh functionalities
 @return the FE postprocessing
 */
TPZPostProcAnalysis * CreatePostProAnalysis(TPZCompMesh * cmesh, SimulationControl & sim_ctrl);


/**
 Execute the postprocessing step

 @param an the FE analysis
 @param file name for the vtk file
 @param sim_ctrl sim_ctrl control object for postprocessing functionalities
 */
void PosProcess(TPZPostProcAnalysis  * an, SimulationControl & sim_ctrl);

int main()
{
    
    // Common parameters
    SimulationControl common;
    common.Omega_ids.Push(1);
    common.Gamma_ids.Push(2);
    common.Gamma_ids.Push(3);
    common.Gamma_ids.Push(4);
    common.Gamma_ids.Push(5);
    
    // Linear case
    SimulationControl linear_elastic(common);
    linear_elastic.dump_folder  = "linear_model";
    linear_elastic.n_h_levels   = 0;
    linear_elastic.p_order      = 2;
    ComputeApproximation(linear_elastic);
    
    
	return 0;
}

void ComputeApproximation(SimulationControl & sim_ctrl){
 
    // Creating the output directory
    std::string command = "mkdir " + sim_ctrl.dump_folder;
    system(command.c_str());
    
    TPZGeoMesh                  * gmesh = GeomtricMesh(sim_ctrl);
    TPZCompMesh                 * cmesh = ComputationalMesh(gmesh, sim_ctrl);
    TPZElastoPlasticAnalysis    * analysis = CreateAnalysis(cmesh, sim_ctrl);
    
    REAL epsilon = 1.0e-4;
    int n_iter  = 10;
    bool linesearchQ = false;
    bool checkconvQ = false;
//    bool consistenMatrixQ = false;
    
    std::stringstream summary;
    summary << sim_ctrl.dump_folder << "/" "conv" << "_" << sim_ctrl.execution_summary << ".txt";
    std::ofstream report(summary.str().c_str(),std::ios::app);
    analysis->IterativeProcess(report, epsilon, n_iter, linesearchQ, checkconvQ);
    
    TPZPostProcAnalysis * postpro = CreatePostProAnalysis(cmesh, sim_ctrl);
    PosProcess(postpro, sim_ctrl);
    
    return;
    
}

TPZGeoMesh * GeomtricMesh(SimulationControl  & sim_ctrl){
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    std::string dirname = PZSOURCEDIR;
    std::string grid;
    grid = dirname + "/Projects/PlasticVerification/Geometry/geometry.msh";
    
    TPZGmshReader Geometry;
    geomesh = Geometry.GeometricGmshMesh(grid);
    const std::string name("Omega description from gmsh script");
    geomesh->SetName(name);
        
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name   << sim_ctrl.dump_folder << "/" "geometry" << "_" << sim_ctrl.domain_type << ".txt";
    vtk_name    << sim_ctrl.dump_folder << "/" "geometry" << "_" << sim_ctrl.domain_type << ".vtk";
    ofstream textfile(text_name.str().c_str());
    geomesh->Print(textfile);

    int n_div = sim_ctrl.n_h_levels;
    for ( int ref = 0; ref < n_div; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = geomesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = geomesh->ElementVec() [i];
            if (gel->Dimension() != 0) gel->Divide (filhos);
        }//for i
    }//ref
    geomesh->BuildConnectivity();
    
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, vtkfile, true);
    
    return geomesh;
}

TPZCompMesh * ComputationalMesh(TPZGeoMesh * gmesh, SimulationControl  & sim_ctrl){
    
    TPZCompMesh * cmesh = new TPZCompMesh;
    
    int nvolumes = sim_ctrl.Omega_ids.size();
    int nboundaries = sim_ctrl.Gamma_ids.size();
    
#ifdef PZDEBUG
    if (nvolumes != 1) {
        std::cout << "Error:: unable to compute the given case = " << & sim_ctrl << std::endl;
        DebugStop();
    }
#endif
    
    int dimension = gmesh->Dimension();
    int p         = sim_ctrl.p_order;
    std::pair<int, TPZVec<STATE> > bc_data;
    
    TPZFNMatrix<10,STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    for (int iv = 0; iv < nvolumes ; iv++) {

        TPZMaterial * volume = MaterialModel(iv,sim_ctrl);
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
            
            int bc_index = sim_ctrl.Gamma_ids[ib];
            bc_data = BcType(bc_index);
            val2(0,0) = bc_data.second[0];
            val2(1,0) = bc_data.second[1];
            TPZBndCond * face = volume->CreateBC(volume,bc_index,bc_data.first,val1,val2);
            cmesh->InsertMaterialObject(face);
        }
        
    }
    cmesh->SetReference(gmesh);
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(p);
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    
    cmesh->AutoBuild();
    cmesh->ExpandSolution();
    
    
//    TPZCompMeshTools::GroupElements(cmesh);
//    TPZCompMeshTools::CreatedCondensedElements(cmesh, true);
    
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_ctrl.dump_folder << "/" << "cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
}

std::pair<int, TPZVec<STATE> > BcType(int bc_index){
    
    TPZVec<STATE> bc_quantity(2);
    std::map<int, std::pair<int, TPZVec<STATE> > > bc_definition;
    
    bc_quantity[0] = 0.0;
    bc_quantity[1] = 0.0;
    bc_definition.insert(std::make_pair(2, std::make_pair(0, bc_quantity)));
    
    bc_quantity[0] = 0.0;
    bc_quantity[1] = 0.0;
    bc_definition.insert(std::make_pair(3, std::make_pair(1, bc_quantity)));
    
    bc_quantity[0] = 0.0;
    bc_quantity[1] = -0.5;
    bc_definition.insert(std::make_pair(4, std::make_pair(0, bc_quantity)));
    
    bc_quantity[0] = 0.0;
    bc_quantity[1] = 0.0;
    bc_definition.insert(std::make_pair(5, std::make_pair(1, bc_quantity)));
    
    std::map<int, std::pair<int, TPZVec<STATE> > >::iterator chunk = bc_definition.find(bc_index);
    if(chunk == bc_definition.end()){
        std::cout<<"Error:: Material not found"<<std::endl;
        DebugStop();
    }
    return chunk->second;
}

TPZMaterial * MaterialModel(int index, SimulationControl  & sim_ctrl){
    
    bool planestrain = true;
    
    /////////////////////////////////////////////////////////////////////////
    // Elastic part
    REAL Eyoung = 4.3614e9;
    REAL nu     = 0.2;
    TPZMatElastoPlastic2D<TPZElasticCriteria> * elastic = new TPZMatElastoPlastic2D< TPZElasticCriteria >(sim_ctrl.Omega_ids[index], planestrain);
    
    TPZElasticCriteria MatEla;
    TPZElasticResponse ER;
    ER.SetUp(Eyoung, nu);
    MatEla.SetElasticResponse(ER);
    elastic->SetPlasticity(MatEla);
    /////////////////////////////////////////////////////////////////////////
    
    
    
    
    TPZMaterial * plastic(elastic);
    return plastic;
}

TPZElastoPlasticAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationControl & sim_ctrl){
    
    std::ofstream out;
    TPZElastoPlasticAnalysis * analysis = new TPZElastoPlasticAnalysis(cmesh, out);
    TPZSkylineStructMatrix matrix(cmesh);
    matrix.SetNumThreads(sim_ctrl.n_threads);
    analysis->SetStructuralMatrix(matrix);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    analysis->SetSolver(step);
    return analysis;
    
}

TPZPostProcAnalysis * CreatePostProAnalysis(TPZCompMesh * cmesh, SimulationControl & sim_ctrl){
    
    TPZPostProcAnalysis * analysis = new TPZPostProcAnalysis;
    
    int n_vols = sim_ctrl.Omega_ids.size();
    TPZVec<int> mat_ids(n_vols);
    for (int i = 0; i < n_vols; i++) {
        mat_ids[i] = sim_ctrl.Omega_ids[i];
    }
    
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Displacement");
    
    TPZStack<std::string> PostProcVars;
    for (int i = 0; i < scalnames.size(); i++) {
        PostProcVars.Push(scalnames[i]);
    }
    for (int i = 0; i < vecnames.size(); i++) {
        PostProcVars.Push(vecnames[i]);
    }

    analysis->SetCompMesh(cmesh);
    analysis->SetPostProcessVariables(mat_ids, PostProcVars);
    TPZFStructMatrix structmatrix(analysis->Mesh());
    analysis->SetStructuralMatrix(structmatrix);
    analysis->TransferSolution();

    return analysis;
}

void PosProcess(TPZPostProcAnalysis  * an, SimulationControl & sim_ctrl)
{
    
    std::stringstream sol_vtk_name;
    sol_vtk_name << sim_ctrl.dump_folder << "/" "sol" << "_" << sim_ctrl.execution_summary << ".vtk";
    std::string file(sol_vtk_name.str());
    
    int dim = 2;
    TPZStack<std::string,10> scalnames, vecnames;
    int div = 1;
    
    vecnames.Push("Displacement");
    
    an->DefineGraphMesh(dim,scalnames,vecnames,file);
    an->PostProcess(div);
}
