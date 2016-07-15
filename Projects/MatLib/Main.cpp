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
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"
#include "pzgeoelbc.h"
#include "TPZRefPatternTools.h"

#include "pzgengrid.h"
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


#include "TPZSpStructMatrix.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"

#include "pznonlinanalysis.h"

#include "pzpostprocanalysis.h"
#include "tpzcompmeshreferred.h"

#ifdef USING_MATLIB // NS: Shouldnt there be a define in here? To choose using the lib matlib?
#include "MatLib.h"
#include "TPZMatLibMaterial.h"
USING_MATLIB_NAMESPACE

#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.main.matlib"));
#endif

// Dummy Boundary Conditions
//const int dirichlet = 0;
//const int neumann = 1;

/// survival test
void TestMatLib(TPZMatLibMaterial &object);

/// generate a square geometric mesh on a unit square
TPZGeoMesh *GenerateGeoMesh(int nelx, int nely);

/// refine the geometry towards a boundary condition
void RefineTowards(TPZGeoMesh *gmesh, int bcmatid, int numref);

/// generate a computation H1 mesh
TPZCompMesh *GenerateCompMesh(TPZGeoMesh *gmesh, int porder);

/// Solve the nonlinear problem
void SolveNonLinear(TPZCompMesh *cmesh);

/// Create the post processing file
void PostProcess(TPZCompMesh *cmesh, const std::string &filename, int resolution);
#endif


int main(int argc, char *argv[])
{

//    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
//    std::string FileName = dirname;
//    FileName = dirname + "/Projects/Elasticity2D/";
//    FileName += "Elasticity2DLog.cfg";
    InitializePZLOG();
#endif

    gRefDBase.InitializeRefPatterns();
#ifdef USING_MATLIB
    std::cout << "Printing something\n";
    ModelDictionary::list(std:: cout);
    
    /**
     *  Create a constitutive model for MatLib
     *
     *  @param "ISOTROPIC_ELASTICITY" Example model, in this case isotropic elastic
     *
     *  @return A dynamically created object
     */
    ConstitutiveModel* model = ModelDictionary::build("ISOTROPIC_ELASTICITY");
    MaterialProperties prop("steel");
    prop.setProperty("YOUNG_MODULUS",2.1e11);
    prop.setProperty("POISSON_COEFFICIENT",0.3);
    MaterialModel mater(*model,prop);
    mater.initialize();
    

    MaterialState locState0,locState1;
    mater.initState(locState0);
    mater.initState(locState1);
    locState1.grad[0] += 1.e-3;
    ParameterSet extPar;
    MatLibMatrix K(model->nExtVar());
    std::cout << "eps=" << locState1.grad << std::endl;
    bool shouldcomputetangent = true;
    mater.updateState(extPar,locState0,locState1,1.0,K,shouldcomputetangent);
    std::cout << "sig=" << locState1.flux << std::endl;
    std::cout << "Internal=" << locState1.internal << std::endl;
    std::cout << "K="<< K << std::endl;
    
    TPZMatLibMaterial mat(1,3,mater);
    TestMatLib(mat);
    int nelx = 4, nely = 4;
    TPZGeoMesh *gmesh = GenerateGeoMesh(nelx, nely);
    
    int numref = 10;
    RefineTowards(gmesh, -8, numref);
    
    int porder = 2;
    TPZCompMesh *cmesh = GenerateCompMesh(gmesh,porder);
    SolveNonLinear(cmesh);
    int resolution = 2;
    PostProcess(cmesh, "../matlib.vtk", resolution);
    
#endif
  
  return 0;
}

void TestMatLib(TPZMatLibMaterial &object)
{
    long memoryindex = object.PushMemItem();
    
    // test in 3d
    int dim = 3;
    TPZMaterialData data;
    TPZFNMatrix<9,STATE> ek(dim,dim,0.), ef(dim,1,0.);
    
    dim = 3;
    object.SetDimension(dim);
    memoryindex = object.PushMemItem();
    data.dphi.Redim(dim,1);
    data.dphix.Redim(dim,1);
    data.dphix(0,0) = 1.;
    data.axes.Redim(dim, 3);
    for (int i=0; i<dim; i++) {
        data.axes(i,i) = 1.;
    }
    data.dsol[0].Redim(dim, dim);
    data.intGlobPtIndex  = memoryindex;
    ek.Redim(dim,dim);
    ef.Redim(dim,1);
    
    object.Contribute(data, 1., ek, ef);
    
    dim = 2;
    object.SetDimension(dim);
    memoryindex = object.PushMemItem();
    data.dphi.Redim(dim,1);
    data.dphix.Redim(dim,1);
    data.dphix(0,0) = 1.;
    data.axes.Redim(dim, 3);
    for (int i=0; i<dim; i++) {
        data.axes(i,i) = 1.;
    }
    data.dsol[0].Redim(dim, dim);
    data.intGlobPtIndex  = memoryindex;
    ek.Redim(dim,dim);
    ef.Redim(dim,1);
    
    object.Contribute(data, 1., ek, ef);
    
    dim = 1;
    object.SetDimension(dim);
    memoryindex = object.PushMemItem();
    data.dphi.Redim(dim,1);
    data.dphix.Redim(dim,1);
    data.dphix(0,0) = 1.;
    data.axes.Redim(dim, 3);
    for (int i=0; i<dim; i++) {
        data.axes(i,i) = 1.;
    }
    data.dsol[0].Redim(dim, dim);
    data.intGlobPtIndex  = memoryindex;
    ek.Redim(dim,dim);
    ef.Redim(dim,1);
    
    object.Contribute(data, 1., ek, ef);
    
}

/// generate a square geometric mesh on a unit square
TPZGeoMesh *GenerateGeoMesh(int nelx, int nely)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    TPZManVector<int,2> nel(2);
    nel[0] = nelx+nelx%2;
    nel[1] = nely+nely%2;
    TPZManVector<REAL,3> x0(3,0.), x1(3,1.);
    x1[2] = 0.;
    TPZGenGrid gengrid(nel, x0, x1);
    gengrid.SetRefpatternElements(true);
    gengrid.Read(gmesh);
    gengrid.SetBC(gmesh, 0, -1);
    gengrid.SetBC(gmesh, 1, -1);
    gengrid.SetBC(gmesh, 4, -4);
    gengrid.SetBC(gmesh, 5, -5);
    gengrid.SetBC(gmesh, 6, -6);
    gengrid.SetBC(gmesh, 7, -7);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZManVector<REAL,3> co(3,0.),qsi(2,0.);
    co[0] = 0.5;
    co[1] = 0.5;
    TPZManVector<long> nodeindexes(1,-1);
    long initindex = 0;
    TPZGeoEl *gelpoint = gmesh->FindElement(co, qsi, initindex, 2);
    int side = gelpoint->WhichSide(qsi);
    TPZGeoElBC(gelpoint,side,-8);
    return gmesh;
}

/// refine the geometry towards a boundary condition
void RefineTowards(TPZGeoMesh *gmesh, int bcmatid, int numref)
{
    std::set<int> matids;
    matids.insert(bcmatid);
    for (int ref = 0; ref<numref; ref++)
    {
        long nel = gmesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            /** @brief Refines the element if it touches an element with a material id included in matids */
            TPZRefPatternTools::RefineDirectional(gel, matids);
        }
    }
}


/// generate a computation H1 mesh
TPZCompMesh *GenerateCompMesh(TPZGeoMesh *gmesh, int porder)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder);
    ConstitutiveModel* model = ModelDictionary::build("ISOTROPIC_ELASTICITY");
    MaterialProperties prop("steel");
    prop.setProperty("YOUNG_MODULUS",2.1e11);
    prop.setProperty("POISSON_COEFFICIENT",0.3);
    MaterialModel mater(*model,prop);
    mater.initialize();
    int matid = 1;
    int dimension = 2;
    TPZMatLibMaterial *material = new TPZMatLibMaterial(matid,dimension,mater);
    
    TPZFNMatrix<4,REAL> val1(2,2,0.),val2(2,1,0.);
    TPZBndCond *bnd0 = material->CreateBC(material, -1, 0, val1, val2);
    val2(0) = -1.;
    TPZBndCond *bnd7 = material->CreateBC(material, -7, 1, val1, val2);
    val2.Zero();
    val2(0) = 1.;
    TPZBndCond *bnd5 = material->CreateBC(material, -5, 1, val1, val2);
    val2.Zero();
    
    val2(1) = 1;
    TPZBndCond *bnd4 = material->CreateBC(material, -4, 3, val1, val2);
    val2.Zero();
    
    val2(1) = 1000.;
    TPZBndCond *bnd8 = material->CreateBC(material, -8, 1, val1, val2);
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(bnd0);
    cmesh->InsertMaterialObject(bnd7);
    cmesh->InsertMaterialObject(bnd5);
    cmesh->InsertMaterialObject(bnd4);
    cmesh->InsertMaterialObject(bnd8);
    
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    return cmesh;
}

/// Solve the nonlinear problem
void SolveNonLinear(TPZCompMesh *cmesh)
{
    TPZNonLinearAnalysis analysis(cmesh, std::cout);
    TPZSkylineStructMatrix strmat(cmesh);
    strmat.SetNumThreads(0);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    
    analysis.SetStructuralMatrix(strmat);
    analysis.SetSolver(step);
    
    analysis.IterativeProcess(std::cout, 1.e-12, 3);
    TPZMatLibMaterial *matlib = dynamic_cast<TPZMatLibMaterial *>(cmesh->FindMaterial(1));
    if (!matlib) {
        DebugStop();
    }
    matlib->LoadState();
    analysis.IterativeProcess(std::cout, 1.e-12, 3);
    
    
}

/// Create the post processing file
void PostProcess(TPZCompMesh *cmesh, const std::string &filename, int resolution)
{
//    TPZCompMeshReferred *refcomp = new TPZCompMeshReferred(cmesh->Reference());
//    refcomp->LoadReferred(cmesh);
    TPZPostProcAnalysis ppanalysis;
    ppanalysis.SetCompMesh(cmesh);
    TPZManVector<int,3> matids(1,1);
    TPZStack<std::string> varnames,scalnames,vecnames,tensornames;
    varnames.Push("stress");
    varnames.Push("deformation");
    varnames.Push("elastically stored energy");
    tensornames.Push("stress");
    tensornames.Push("deformation");
    scalnames.Push("elastically stored energy");
    varnames.Push("state");
    vecnames.Push("state");
    ppanalysis.SetPostProcessVariables(matids, varnames);
    TPZSkylineStructMatrix strmat(ppanalysis.Mesh());
    ppanalysis.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    ppanalysis.SetSolver(step);
    ppanalysis.TransferSolution();
    
    ppanalysis.DefineGraphMesh(2, scalnames, vecnames, tensornames, filename);
    ppanalysis.PostProcess(resolution);
//    delete refcomp;
}

