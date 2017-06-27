//
//  GPU_poisson.cpp
//  Experiemental gpu implementation for Poisson equation
//
//
//

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZRefPatternTools.h"

#include "tpzhierarquicalgrid.h"
#include "TPZReadGIDGrid.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"
#include "tpzquadraticquad.h"
#include "tpzarc3d.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzfunction.h"
#include "tpzchangeel.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "TPZPrimalPoisson.h"
#include "TPZDualPoisson.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZCompMeshTools.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"

#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckmesh.h"
#include "TPZGmshReader.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


/**
 *  Structure containing the project directives and controls.
 */
struct SimulationCase {
    
    bool            UsePardisoQ;
    bool            UseFrontalQ;
    int             n_h_levels;
    int             n_p_levels;
    int             int_order;
    int             n_threads;
    std::string     domain_type;
    std::string     conv_summary;
    std::string     dump_folder;
    TPZStack<int>   omega_ids;
    TPZStack<int>   gamma_ids;
    
    SimulationCase() : UsePardisoQ(true), UseFrontalQ(false), n_h_levels(0), n_p_levels(1), int_order(1), n_threads(0),
    domain_type(""),conv_summary(""),dump_folder(""),omega_ids(),gamma_ids()
    {
        
    }
    
    SimulationCase(const SimulationCase &copy) : UsePardisoQ(copy.UsePardisoQ), UseFrontalQ(copy.UseFrontalQ),
        n_h_levels(copy.n_h_levels), n_p_levels(copy.n_p_levels), int_order(copy.int_order),
        n_threads(copy.n_threads), domain_type(copy.domain_type), conv_summary(copy.conv_summary),
        dump_folder(copy.dump_folder), omega_ids(copy.omega_ids), gamma_ids(copy.gamma_ids)
    {
        
    }
    
    SimulationCase &operator=(const SimulationCase &copy)
    {

        UsePardisoQ = copy.UsePardisoQ;
        UseFrontalQ = copy.UseFrontalQ;
        n_h_levels = copy.n_h_levels;
        n_p_levels = copy.n_p_levels;
        int_order = copy.int_order;
        n_threads = copy.n_threads;
        domain_type = copy.domain_type;
        conv_summary = copy.conv_summary;
        dump_folder = copy.dump_folder;
        omega_ids = copy.omega_ids;
        gamma_ids = copy.gamma_ids;
        return *this;
    }
};

//#define Solution1
#define Solution2


static void Analytic(const TPZVec<REAL> &x, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu);
static void f(const TPZVec<REAL> &p, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf);

TPZGeoMesh * GeomtricMesh(int ndiv, SimulationCase  & sim_data);
void PrintGeometry(TPZGeoMesh * gmesh, SimulationCase & sim_data);
void UniformRefinement(TPZGeoMesh * gmesh, int n_ref);


TPZGeoMesh * MakeCartesianGeometry(int ndiv, SimulationCase  & sim_data);
void Parametricfunction_x(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void Parametricfunction_y(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void Parametricfunction_z(const TPZVec<REAL> &par, TPZVec<REAL> &X);


void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int Axis);


TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, long &ndof, TPZVec<TPZCompMesh *> &meshvec);

TPZCompMesh * PrimalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, long &ndof); //  Primal approximation


TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
void PosProcess(TPZAnalysis * an, std::string file, SimulationCase & sim_data);

// Compute simulation cases
void ComputeCases(TPZStack<SimulationCase> cases);

// Compute simulation cases
void ComputeApproximation(SimulationCase & sim_data);

// Approximation ratess
void ComputeConvergenceRates(TPZVec<REAL> &error, TPZVec<REAL> &convergence);

// Error calculation
void ErrorH1(TPZCompMesh *cmesh, REAL &error_primal , REAL & error_dual, REAL & error_h1);

int main()
{
   
    // load uniform refinement patterns
    gRefDBase.InitializeAllUniformRefPatterns();
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    TPZStack<SimulationCase> simulations;
    
    // Formulations over the sphere 
    struct SimulationCase common;
    common.UsePardisoQ = false;
    common.UseFrontalQ = false;
    common.n_h_levels = 4;
    common.n_p_levels = 1;
    common.int_order  = 10; // needed for rhs function
    common.n_threads  = 0;
    common.conv_summary = "convergence_summary";
    common.omega_ids.Push(1);     // Domain
    common.gamma_ids.Push(-1);    // Gamma_D Dirichlet data boundary
    
    // Primal Formulation over the solid sphere
    struct SimulationCase H1Case_1 = common;
    H1Case_1.domain_type = "cube"; // domain_type = {line, plane, otherwise -> cube}
    H1Case_1.dump_folder = "H1_" + H1Case_1.domain_type;
    simulations.Push(H1Case_1);
    
    ComputeCases(simulations);
    return 0;
}

void ComputeCases(TPZStack<SimulationCase> cases){
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    int n_cases = cases.size();
    for (int i = 0; i < n_cases; i++) {
        ComputeApproximation(cases[i]);
    }
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    std::cout << "End:: Overal time = " << int_t2-int_t1 << std::endl;
#endif
    
}

void ComputeApproximation(SimulationCase & sim_data){
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_case_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    // Creating the directory
    std::string command = "mkdir " + sim_data.dump_folder;
    system(command.c_str());
    
    std::stringstream summary;
    summary   << sim_data.dump_folder << "/" "conv" << "_" << sim_data.domain_type << ".txt";
    std::ofstream convergence(summary.str().c_str(),std::ios::app);
    
    TPZManVector<REAL,10> p_error(sim_data.n_h_levels+1,1.0);
    TPZManVector<REAL,10> d_error(sim_data.n_h_levels+1,1.0);
    TPZManVector<REAL,10> h_error(sim_data.n_h_levels+1,1.0);
    
    TPZManVector<REAL,10> p_conv(sim_data.n_h_levels,0.0);
    TPZManVector<REAL,10> d_conv(sim_data.n_h_levels,0.0);
    TPZManVector<REAL,10> h_conv(sim_data.n_h_levels,0.0);
    
    int n_h_levels = sim_data.n_h_levels;
    int n_p_levels = sim_data.n_p_levels;

    using namespace std;
    
    for (int p = 1; p <= n_p_levels; p++) {
        
        convergence << std::endl;        
        convergence << " Polynomial order  =  " << p << std::endl;
        convergence << setw(5)  << " h" << setw(25) << " ndof" << setw(25) << " ndof_cond" << setw(25) << " assemble_time (msec)" << setw(25) << " solving_time (msec)" << setw(25) << " error_time (msec)" << setw(25) << " Primal l2 error" << setw(25) << " Dual l2 error "  << setw(25) << " H1 error " << endl;
        
        int h_base = 1; // start with 8 elements
        for (int h = 0; h <= n_h_levels; h++) {
            
            // Compute the geometry
            TPZGeoMesh * gmesh;
            gmesh = GeomtricMesh(h_base, sim_data);
    
#ifdef PZDEBUG
            TPZCheckGeom check(gmesh);
            int checkQ = check.PerformCheck();
            if (checkQ) {
                DebugStop();
            }
#endif
            
            UniformRefinement(gmesh, h);

            std::stringstream text_name;
            std::stringstream vtk_name;
            text_name   << sim_data.dump_folder << "/" "geo" << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".txt";
            vtk_name    << sim_data.dump_folder << "/" "geo" << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".vtk";
            ofstream textfile(text_name.str().c_str());
            gmesh->Print(textfile);
            
            std::ofstream vtkfile(vtk_name.str().c_str());
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
            
            // Compute the geometry
            long ndof, ndof_cond;
            TPZManVector<TPZCompMesh *,5> meshvec;
            TPZCompMesh * cmesh = ComputationalMesh(gmesh, p, sim_data, ndof, meshvec);

            // Create Analysis
            TPZAnalysis * analysis = CreateAnalysis(cmesh,sim_data);
            
            REAL assemble_time = 0.0;
            REAL solving_time  = 0.0;
            REAL error_time    = 0.0;
            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            
            analysis->Assemble();
            
            ndof_cond = analysis->Rhs().Rows();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            
            analysis->Solve();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
            assemble_time = boost::numeric_cast<double>((t2-t1).total_milliseconds());
            solving_time  = boost::numeric_cast<double>((t3-t2).total_milliseconds());
#endif
            
            analysis->LoadSolution();
            cmesh->Solution() *= -1.0; /* @omar::consequence of newton correction */
            analysis->LoadSolution(cmesh->Solution());
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmesh);

            // PostProccessing
            std::stringstream sol_vtk_name;
            sol_vtk_name    << sim_data.dump_folder << "/" "sol" << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".vtk";
            std::string file(sol_vtk_name.str());
            PosProcess(analysis, file, sim_data);
            
            
            // compute the error
#ifdef USING_BOOST
            boost::posix_time::ptime error_t1 = boost::posix_time::microsec_clock::local_time();
#endif
            
            ErrorH1(cmesh, p_error[h], d_error[h], h_error[h]);
            
#ifdef USING_BOOST
            boost::posix_time::ptime error_t2 = boost::posix_time::microsec_clock::local_time();
#endif
            
#ifdef USING_BOOST
            error_time = boost::numeric_cast<double>((error_t2 - error_t1).total_milliseconds());
            
#endif
            // current summary
            convergence << setw(5) << h << setw(25) << ndof << setw(25) << ndof_cond << setw(25) << assemble_time << setw(25) << solving_time << setw(25) << error_time << setw(25) << p_error[h] << setw(25) << d_error[h]  << setw(25) << h_error[h] << endl;
            
            delete cmesh;
            for (int i = 0; i < meshvec.size(); i++) {
                delete meshvec[i];
            }
            delete gmesh;
        }
        
        // compute rates
        ComputeConvergenceRates(p_error,p_conv);
        ComputeConvergenceRates(d_error,d_conv);
        ComputeConvergenceRates(h_error,h_conv);
        
        
        // print convergence summary
        convergence << std::endl;
        convergence << " Convergence rates summary " << std::endl;
        convergence << " Polynomial order  =  " << p << std::endl;
        convergence << " Primal convergence rates = " << setw(5) << p_conv << std::endl;
        convergence << " Dual convergence rates = " << setw(5) << d_conv << std::endl;
        convergence << " H1 convergence rates = " << setw(5) << h_conv << std::endl;
        convergence << std::endl;
        convergence << " ------------------------------------------------------------------ " << std::endl;
        convergence.flush();
        
    }
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_case_t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    REAL case_solving_time = boost::numeric_cast<double>((int_case_t2-int_case_t1).total_milliseconds());
#endif
    
#ifdef USING_BOOST
    convergence <<  "Case closed in :" << setw(10) <<  case_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << "Case closed in :" << setw(10) <<  case_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << std::endl;
#endif
    
}

void ComputeConvergenceRates(TPZVec<REAL> &error, TPZVec<REAL> &convergence){
    
    int ndata = error.size();
    STATE log0p5 = log(0.5);
    for (int i = 1; i < ndata; i++) {
        STATE logerror = log(error[i-1]);
        STATE logerrori = log(error[i]);
        convergence[i-1] = (logerrori - logerror)/log0p5;
    }
}


void Analytic(const TPZVec<REAL> &p, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu){
    

    gradu.Resize(4,1);
    gradu.Zero();
    u.resize(1);
    
    STATE x,y,z;
    x = p[0];
    y = p[1];
    z = p[2];
   
#ifdef Solution1
    
    u[0] = x*x + y*y + z*z;
    
    STATE dudx  = 2.0*x;
    STATE dudy  = 2.0*y;
    STATE dudz  = 2.0*z;
    
    gradu(0,0) = -1.0*(dudx);
    gradu(1,0) = -1.0*(dudy);
    gradu(2,0) = -1.0*(dudz);
    
    gradu(3,0) = -6.0;
    
#endif
    
#ifdef Solution2
    
    REAL l = 1.0;
    
    u[0] = sin(M_PI*x*l)*sin(M_PI*x*l) * sin(M_PI*y*l)*sin(M_PI*y*l) * sin(M_PI*z*l)*sin(M_PI*z*l);
    
    STATE dudx  = M_PI*l*(sin(2.0*M_PI*x*l) * sin(M_PI*y*l)*sin(M_PI*y*l) * sin(M_PI*z*l)*sin(M_PI*z*l));
    STATE dudy  = M_PI*l*(sin(M_PI*x*l)*sin(M_PI*x*l) * sin(2.0*M_PI*y*l) * sin(M_PI*z*l)*sin(M_PI*z*l));
    STATE dudz  = M_PI*l*(sin(M_PI*x*l)*sin(M_PI*x*l) * sin(M_PI*y*l)*sin(M_PI*y*l) * sin(2.0*M_PI*z*l));
    
    gradu(0,0) = -1.0*(dudx);
    gradu(1,0) = -1.0*(dudy);
    gradu(2,0) = -1.0*(dudz);
    
    REAL f1 = cos(2.0*M_PI*x*l) * sin(M_PI*y*l)*sin(M_PI*y*l) * sin(M_PI*z*l)*sin(M_PI*z*l);
    REAL f2 = sin(M_PI*x*l)*sin(M_PI*x*l) * cos(2.0*M_PI*y*l) * sin(M_PI*z*l)*sin(M_PI*z*l);
    REAL f3 = sin(M_PI*x*l)*sin(M_PI*x*l) * sin(M_PI*y*l)*sin(M_PI*y*l) * cos(2.0*M_PI*z*l);
    
    gradu(3,0) = - 2.0 * M_PI * M_PI * l * l * (f1 + f2 + f3);
    
#endif
    
    
}

void f(const TPZVec<REAL> &p, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf){
    
    f.resize(1);
    TPZVec<STATE> u;
    TPZFMatrix<STATE> gradu;
    Analytic(p, u, gradu);
    f[0] = gradu(3,0);
    
}

TPZGeoMesh * GeomtricMesh(int ndiv, SimulationCase  & sim_data){
    
    TPZGeoMesh * geometry = NULL;
    
    geometry = MakeCartesianGeometry(ndiv, sim_data);
    return geometry;
}

TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, long &ndof, TPZVec<TPZCompMesh *> &meshvec){
    
    TPZCompMesh * mesh = PrimalMesh(geometry, p, sim_data, ndof);
    return mesh;

}

TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    
    if (sim_data.UsePardisoQ) {
        
        TPZSymetricSpStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    
    if (sim_data.UseFrontalQ) {

        TPZParFrontStructMatrix<TPZFrontSym<STATE> > matrix(cmesh);
        matrix.SetDecomposeType(ELDLt);
        matrix.SetNumThreads(sim_data.n_threads);
        
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    else{
        
        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);        
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
    }
    
   return analysis;
    
}

void PosProcess(TPZAnalysis  * an, std::string file, SimulationCase & sim_data)
{
    int dim = an->Mesh()->Reference()->Dimension();
    TPZStack<std::string,10> scalnames, vecnames;
    
    vecnames.Push("q");
    vecnames.Push("q_exact");
    scalnames.Push("p");
    scalnames.Push("p_exact");
    scalnames.Push("f_exact");

    int div = 1;
    an->DefineGraphMesh(dim,scalnames,vecnames,file);
    an->PostProcess(div,dim);
    
}

TPZCompMesh * PrimalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, long &ndof){
    
    int dimension = geometry->Dimension();
    int dirichlet = 0;
    int nvolumes = sim_data.omega_ids.size();
    int nboundaries = sim_data.gamma_ids.size();
    
#ifdef PZDEBUG
    if (nvolumes != 1) {
    std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
        DebugStop();
    }
#endif
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
        
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    for (int iv = 0; iv < nvolumes ; iv++) {
        
        TPZMaterial * volume = new TPZPrimalPoisson(sim_data.omega_ids[iv], dimension);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic);
        analytic->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > analytic_full = analytic;
        volume->SetForcingFunctionExact(analytic_full);
        
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
            TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ib],dirichlet,val1,val2);
            cmesh->InsertMaterialObject(face);
        }
        
    }
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(p);
    cmesh->SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    cmesh->ExpandSolution();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    ndof = cmesh->NEquations();
    
//    TPZCompMeshTools::GroupElements(cmesh);
//    TPZCompMeshTools::CreatedCondensedElements(cmesh, true);
//    cmesh->CleanUpUnconnectedNodes();
//    cmesh->ExpandSolution();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Primal_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}

TPZGeoMesh * MakeCartesianGeometry(int ndiv, SimulationCase  & sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 1) {
        std::cout << "Cube:: Please pass materials ids for volume (1) and boundary domain (-1) }" << std::endl;
        DebugStop();
    }
#endif
    
    REAL t=0.0;
    int n = pow(2,ndiv);
    REAL dt = 1.0/REAL(n);

    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh_point = new TPZGeoMesh;
    GeoMesh_point->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh_point->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    int matid=1;
    int front = -1;
    int back  = -1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh_point);
    GeoMesh_point->BuildConnectivity();
    GeoMesh_point->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh_point);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(Parametricfunction_x);
    CreateGridFrom.SetParametricFunction(ParFunc);
    CreateGridFrom.SetFrontBackMatId(front, front);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh_line = CreateGridFrom.ComputeExtrusion(t, dt, n);

    if (sim_data.domain_type == "line") {
        return GeoMesh_line;
    }
    
    TPZHierarquicalGrid CreateGridFrom2(GeoMesh_line);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(Parametricfunction_y);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    CreateGridFrom2.SetFrontBackMatId(front, front);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh_surface = CreateGridFrom2.ComputeExtrusion(t, dt, n);

    if (sim_data.domain_type == "plane") {
        return GeoMesh_surface;
    }
    
    TPZHierarquicalGrid CreateGridFrom3(GeoMesh_surface);
    GeoMesh_surface->SetDimension(2);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc3 = new TPZDummyFunction<REAL>(Parametricfunction_z);
    CreateGridFrom3.SetParametricFunction(ParFunc3);
    CreateGridFrom3.SetFrontBackMatId(back, back);

    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh_cube = CreateGridFrom3.ComputeExtrusion(t, dt, n);

    return GeoMesh_cube;
}

void Parametricfunction_x(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = par[0];
    X[1] = 0.0;
    X[2] = 0.0;
}

void Parametricfunction_y(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void Parametricfunction_z(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}


void PrintGeometry(TPZGeoMesh * gmesh, SimulationCase & sim_data){
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name   << sim_data.dump_folder << "/" "geo" << "_" << sim_data.domain_type << "_" << "mhm" << "_l" << ".txt";
    vtk_name    << sim_data.dump_folder << "/" "geo" << "_" << sim_data.domain_type << "_" << "mhm" << "_l" << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    gmesh->Print(textfile);
    
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}

/** @brief Apply uniform refinement on the Geometric mesh */
void UniformRefinement(TPZGeoMesh * gmesh, int n_ref){
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = gmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            if (gel->Dimension() != 0) gel->Divide (filhos);
        }//for i
    }//ref
//    gmesh->BuildConnectivity();
}

void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the axis -> i.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    
    switch (Axis) {
        case 1:
        {
            RotationMatrix(0,0) = 1.0;
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(1,2) =   -sin(theta);
            RotationMatrix(2,1) =   +sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 2:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,2) =   +sin(theta);
            RotationMatrix(1,1) = 1.0;
            RotationMatrix(2,0) =   -sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 3:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
        default:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
    }
    
    TPZVec<REAL> iCoords(3,0.0);
    TPZVec<REAL> iCoordsRotated(3,0.0);
    
    //RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

void ErrorH1(TPZCompMesh *cmesh, REAL &error_primal , REAL & error_dual, REAL & error_h1)
{
    
    long nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    TPZManVector<REAL,10> globalerror(3,0.   );
    for (long iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];

        if (!cel) {
            continue;
        }
        
        if(cel->Reference()->Dimension()!=dim) {
            continue;
        }
        
        TPZManVector<REAL,10> elerror(3,0.);
        elerror.Fill(0.);
        cel->EvaluateError(Analytic, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerror[i] += elerror[i]*elerror[i];
            
        }
    }
    
    error_primal    = sqrt(globalerror[0]);
    error_dual      = sqrt(globalerror[1]);
    error_h1        = sqrt(globalerror[2]);
    
    
}
