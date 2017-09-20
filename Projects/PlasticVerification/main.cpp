
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

// Global Matrix
#include "pzskylstrmatrix.h"

// Solver
#include "pzstepsolver.h"

// Materials
#include "TPZElasticCriteria.h"
#include "pzelastoplastic2D.h"
#include "TPZPlasticStepPV.h"


/**
 Structure wich define the simulation controls
 */
struct SimulationControl {
    
    bool            nonlinear_model_Q;
    int             elemen_type;
    int             n_h_levels;
    int             n_p_levels;
    int             int_order;
    int             n_threads;
    std::string     domain_type;
    std::string     execution_summary;
    std::string     dump_folder;
    TPZStack<int>   Omega_ids;
    TPZStack<int>   Gamma_ids;
    
    SimulationControl() : nonlinear_model_Q(false),elemen_type(0), n_h_levels(0), n_p_levels(1), int_order(1), n_threads(0), domain_type(""),execution_summary(""),dump_folder(""),Omega_ids(),Gamma_ids()
    {
        
    }
    
    SimulationControl(const SimulationControl &copy) :
        nonlinear_model_Q(copy.nonlinear_model_Q),
        elemen_type(copy.elemen_type),
        n_h_levels(copy.n_h_levels),
        n_p_levels(copy.n_p_levels),
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
TPZNonLinearAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationControl & sim_ctrl);



/**
 Execute the postprocessing step

 @param an the FE analysis
 @param file name for the vtk file
 @param sim_ctrl sim_ctrl control object for postprocessing functionalities
 */
void PosProcess(TPZAnalysis  * an, SimulationControl & sim_ctrl);

int main()
{
    
    SimulationControl linear_elastic;
    linear_elastic.dump_folder  = "linear_model";
    linear_elastic.n_h_levels   = 0;
    ComputeApproximation(linear_elastic);
    
    
	return 0;
}

void ComputeApproximation(SimulationControl & sim_ctrl){
 
    // Creating the output directory
    std::string command = "mkdir " + sim_ctrl.dump_folder;
    system(command.c_str());
    
    TPZGeoMesh              * gmesh = GeomtricMesh(sim_ctrl);
    TPZCompMesh             * cmesh = ComputationalMesh(gmesh, sim_ctrl);
    TPZNonLinearAnalysis    * analysis = CreateAnalysis(cmesh, sim_ctrl);
    
    REAL epsilon = 0.0001;
    int n_iter  = 20;
    
    std::stringstream summary;
    summary << sim_ctrl.dump_folder << "/" "conv" << "_" << sim_ctrl.execution_summary << ".txt";
    std::ofstream report(summary.str().c_str(),std::ios::app);
    analysis->IterativeProcess(report, epsilon, n_iter);
    
    
    DebugStop();
    PosProcess(analysis, sim_ctrl);
    
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
    int dirichlet = 0;
    int neumman   = 1;
    int p         = sim_ctrl.n_p_levels;
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    for (int iv = 0; iv < nvolumes ; iv++) {
        
        
        TPZMaterial * volume = MaterialModel(iv,sim_ctrl);
        
//        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f);
//        rhs_exact->SetPolynomialOrder(sim_data.int_order);
//        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
//        volume->SetForcingFunction(rhs);
        
        
//        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic);
//        analytic->SetPolynomialOrder(sim_data.int_order);
//        TPZAutoPointer<TPZFunction<STATE> > analytic_full = analytic;
//        volume->SetForcingFunctionExact(analytic_full);
        
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
//            TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution);
//            analytic_bc->SetPolynomialOrder(sim_data.int_order);
//            TPZAutoPointer< TPZFunction<STATE> > solution = analytic_bc;
            
//            TPZMaterial * face = volume->CreateBC(volume,sim_ctrl.Gamma_ids[ib],dirichlet,val1,val2);
//            face->SetForcingFunction(solution);
//            cmesh->InsertMaterialObject(face);
        }
        
    }
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(p);
    cmesh->SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    cmesh->ExpandSolution();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ndof = cmesh->NEquations();
    
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

TPZMaterial * MaterialModel(int index, SimulationControl  & sim_ctrl){
    
    bool planestrain = true;

    TPZMatElastoPlastic2D<TPZElasticCriteria> * elastic = new TPZMatElastoPlastic2D< TPZElasticCriteria >(sim_ctrl.Omega_ids[index], planestrain);
    TPZElasticCriteria MatEla;
    TPZTensor<REAL> initstress(0.), finalstress(0.);
    elastic->SetPlasticity(MatEla);
    TPZMaterial * plastic(elastic);
    
    return plastic;
}

TPZNonLinearAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationControl & sim_ctrl){
    
    std::ofstream out;
    TPZNonLinearAnalysis * analysis = new TPZNonLinearAnalysis(cmesh, out);
    TPZSkylineStructMatrix matrix(cmesh);
    matrix.SetNumThreads(sim_ctrl.n_threads);
    analysis->SetStructuralMatrix(matrix);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    analysis->SetSolver(step);
    return analysis;
    
}

void PosProcess(TPZAnalysis  * an, SimulationControl & sim_ctrl)
{
    
    std::stringstream sol_vtk_name;
    sol_vtk_name << sim_ctrl.dump_folder << "/" "sol" << "_" << sim_ctrl.execution_summary << ".vtk";
    std::string file(sol_vtk_name.str());
    
    int dim = 2;
    TPZStack<std::string,10> scalnames, vecnames;
    int div = 2;
    
    vecnames.Push("u");
    
    an->DefineGraphMesh(dim,scalnames,vecnames,file);
    an->PostProcess(div,dim);
}
