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

#include "pzgengrid.h"
#include "TPZMatElasticity2D.h"
#include "TPZMatLaplacian.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "TPZLagrangeMultiplier.h"
#include "pzmat1dlin.h"

#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"

#include "pzelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "pzmultiphysicselement.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiPhysicsinterfaceEl.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.hpp"

#include <cmath>
#include <set>
//#include <Accelerate/Accelerate.h>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif


//	This Solve Different analysis
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh);

//	These are tools for spatial and polynomial refinement and Postprocess of solutions
void PostProcessElasticity(TPZAnalysis &an, std::string plotfile);

/// create the 1d skeleton geometric elements
void AddSkeletonElements(TPZAutoPointer<TPZGeoMesh> gmesh);

/// create the collapsed volume elements
/// groupval : geometric element index to which a collapsed element belongs
void AddCollapsedVolumeElements(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<long> &groupval);

/// create the 1d skeleton elements
TPZCompMesh * ComputationalSkeletonMesh(TPZAutoPointer<TPZGeoMesh> & gmesh,int pOrder);

// Condense the internal degrees of freedom
void CondenseInternalEquations(TPZCompMesh *mphysics);

/// Define a rotation tensor which will be applied to the boundary conditions and mesh coordinates
TPZFNMatrix<9,REAL> gRotate(3,3,0.);

enum MMATID {Enomat, Emat1, Emat2, Ebc1, Ebc2, Ebc3, Ebc4, Ewrap, ESkeleton, EInterfaceMat1, EInterfaceMat2, EGroup};

void OneTriangle(TPZAutoPointer<TPZGeoMesh> gmesh)
{
    REAL coord[3][3] = {
        {-1,-1,0},
        {1,-1},
        {0,0}
    };
    gmesh->NodeVec().Resize(3);
    for (int i=0; i<3; i++) {
        TPZManVector<REAL> co(3);
        for(int c=0; c<3; c++) co[c] = coord[i][c];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    TPZManVector<long> indices(3);
    for(int i=0; i<3; i++) indices[i] = (i+0)%3;
    int matid = 1;
    long index;
    TPZGeoEl *gel = gmesh->CreateGeoElement(ETriangle, indices, matid, index);
    gmesh->BuildConnectivity();
    gel->CreateBCGeoEl(3, 2);
    gel->CreateBCGeoEl(4, 3);
    gel->CreateBCGeoEl(5, 4);
}

void OneQuad(TPZAutoPointer<TPZGeoMesh> gmesh)
{
    REAL coord[4][3] = {
        {-1,-1,0},
        {1,-1},
        {0,0},
        {0,0}
    };
    gmesh->NodeVec().Resize(4);
    for (int i=0; i<4; i++) {
        TPZManVector<REAL> co(3);
        for(int c=0; c<3; c++) co[c] = coord[i][c];
        gmesh->NodeVec()[i].Initialize(co, gmesh);
    }
    TPZManVector<long> indices(4);
    for(int i=0; i<4; i++) indices[i] = (i+0)%4;
    int matid = 1;
    long index;
    TPZGeoEl *gel = gmesh->CreateGeoElement(EQuadrilateral, indices, matid, index);
    gmesh->BuildConnectivity();
    gel->CreateBCGeoEl(4, 2);
    gel->CreateBCGeoEl(5, 3);
    gel->CreateBCGeoEl(7, 4);
}


void Displacement(const TPZVec<REAL> &xv, TPZVec<STATE> &result)
{
    //1    result[0] = 100.+x[0]*x[0]+3.*x[0]*x[1]+4.*x[1]*x[1];
    //1    result[1] = 200.+5.*x[0]+6.*x[1]+2.*x[0]*x[0]+4.*x[0]*x[1]+6.*x[1]*x[1];
    // apply the rotation transposed
    REAL x = gRotate(0,0)*xv[0]+gRotate(1,0)*xv[1];
    REAL y = gRotate(0,1)*xv[0]+gRotate(1,1)*xv[1];
//2    result[0] = exp(x-y)*(1.-x)*x*(1.-y)*y;
//2    result[1] = sin(M_PI*x)*sin(M_PI*y);
    result[0] = gRotate(0,0)*x;
    result[1] = gRotate(1,0)*x;
    
}

void Force(const TPZVec<REAL> &xv, TPZVec<STATE> &result)
{
    REAL x = xv[0];
    REAL y = xv[1];
    //1    result[0] = 14.;
    //1    result[1] = 61./2.;
    REAL G = 0.5;
    REAL Lambda = 1.;
    
    
//2    result[0]=exp(x-y)*x*(-1.+y)*(G*(4-4.*x+5*y+3*x*y)+(3+x)*y*Lambda)+M_PI*M_PI*(G+Lambda)*cos(M_PI*x)*cos(M_PI*y);
//2    result[1] = -exp(x-y)*(-1.+x+x*x)*(1.+(-3.+y)*y)*(G+Lambda)-M_PI*M_PI*(3*G+Lambda)*sin(M_PI*x)*sin(M_PI*y);
    result[0] = 0;
    result[1] = 0;
}



void Exact_Sol(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL x = xv[0];
    REAL y = xv[1];
    val.resize(6);
    deriv.Resize(2, 6);
    
    REAL G = 0.5;
    REAL Lambda = 1.;
    
    deriv.Zero();
    val[0] = x;
    val[1] = 0;

    if(y<0.5)
    {
        G = 0.5;
        Lambda = 0.;
    }
    else
    {
        G = 5.;
        Lambda = 0.;
    }
    val[2] = 2.*G+Lambda;
    val[3] = 0.;
    val[4] = 0.;
    val[5] = Lambda;
    
    /* 1111111111111 */
    //    val[0] = 100.+x[0]*x[0]+3.*x[0]*x[1]+4.*x[1]*x[1];
    //    val[1] = 200.+5.*x[0]+6.*x[1]+2.*x[0]*x[0]+4.*x[0]*x[1]+6.*x[1]*x[1];
    //    // stresses
    //    val[2] = 6.+8.*x[0]+18.*x[1];
    //    val[3] = 2.5+3.5*x[0]+6.*x[1];
    //    val[4] = val[3];
    //    val[5] = 12.+10.*x[0]+27.*x[1];
    //
    //    deriv(0,0) = 2.*x[0]+3.*x[1];
    //    deriv(1,0) = 3.*x[0]+8.*x[1];
    //
    //    deriv(0,1) = 5.+4.*x[0]+4.*x[1];
    //    deriv(1,1) = 6.+4.*x[0]+12.*x[1];
    //
    //    deriv(0,2) = 8.;
    //    deriv(1,2) = 18.;
    //    deriv(0,3) = 3.5;
    //    deriv(1,3) = 6.;
    //    deriv(0,4) = deriv(0,3);
    //    deriv(1,4) = deriv(1,3);
    //    deriv(0,5) = 10.;
    //    deriv(1,5) = 27.;
    /* 1111111111111 */
    /* 2222222222222
    val[0] = exp(x-y)*(1.-x)*x*(1.-y)*y;
    val[1] = sin(M_PI*x)*sin(M_PI*y);
    
    val[2] = (2.*G+Lambda)*exp(x-y)*(-1.+x+x*x)*(-1.+y)*y+M_PI*cos(M_PI*y)*sin(M_PI*x);
    val[3] = G*(-exp(x-y)*(-1.+x)*x*(1.+(-3.+y)*y)+M_PI*cos(M_PI*x)*sin(M_PI*y));
    val[4] = val[3];
    val[5] = exp(x-y)*(-1.+x+x*x)*(-1.+y)*y*Lambda+(2.*G+Lambda)*M_PI*cos(M_PI*y)*sin(M_PI*x);
    
    deriv(0,0) = exp(x-y)*(-1.+x+x*x)*(-1.+y)*y;
    deriv(1,0) = -exp(x-y)*(-1+x)*x*(1.-3.*y+y*y);
    
    deriv(0,1) = M_PI*cos(M_PI*x)*sin(M_PI*y);
    deriv(1,1) = M_PI*cos(M_PI*y)*sin(M_PI*x);
    
    deriv(0,2) = exp(x-y)*x*(3.+x)*(-1.+y)*y*(2.*G+Lambda)+M_PI*M_PI*Lambda*cos(M_PI*x)*cos(M_PI*y);
    deriv(0,5) = exp(x-y)*x*(3.+x)*(-1.+y)*y*Lambda + M_PI*M_PI*(2.*G+Lambda)*cos(M_PI*x)*cos(M_PI*y);
    deriv(0,3) = -exp(x-y)*G*(-1.+x+x*x)*(1.+(-3.+y)*y)-G*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    deriv(0,4) = deriv(0,3);
    
    deriv(1,2) = -exp(x-y)*(-1.+x+x*x)*(1.+(-3.+y)*y)*(2.*G+Lambda)-M_PI*M_PI*Lambda*sin(M_PI*x)*sin(M_PI*y);
    deriv(1,5) = -exp(x-y)*(-1.+x+x*x)*(1.+(-3.+y)*y)*Lambda-M_PI*M_PI*(2.*G+Lambda)*sin(M_PI*x)*sin(M_PI*y);
    deriv(1,3) = exp(x-y)*G*(-1.+x)*x*(-4.+y)*(-1.+y)+G*M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y);
    deriv(1,4) = deriv(1,3);
    22222222222222 */
}


int main(int argc, char *argv[])
{
    
    gRotate.Identity();
    
    
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    std::string FileName = dirname;
    FileName = dirname + "/Projects/Elasticity2D/";
    FileName += "Elasticity2DLog.cfg";
    InitializePZLOG();
#endif
    int PElasticity = 2;
    int maxnelx = 2;
    for ( int nelx = 1; nelx < maxnelx; nelx *= 2)
    {
        
        std::cout << "nelx = " << nelx << std::endl;
        TPZManVector<REAL,4> x0(3,-1.),x1(3,1.);
        x0[2] = 0.;
        x1[2] = 0.;
        TPZManVector<int,4> nx(2,nelx);
        TPZGenGrid gengrid(nx,x0,x1);
        gengrid.SetElementType(EQuadrilateral);
        TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
        
//        OneQuad(gmesh);
        gengrid.Read(gmesh,EGroup);
        gengrid.SetBC(gmesh, 4, Ebc1);
        gengrid.SetBC(gmesh, 5, Ebc2);
        gengrid.SetBC(gmesh, 6, Ebc3);
        gengrid.SetBC(gmesh, 7, Ebc4);
        
        AddSkeletonElements(gmesh);
        /// generate the SBFem elementgroups
        
        /// put sbfem pyramids into the element groups
        
        TPZCompMesh *mphysics = ComputationalSkeletonMesh(gmesh,1);
        
        std::map<long,long> geltogroup;
        long nel = gmesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel) continue;
            if (gel->MaterialId() == EGroup) {
                long index;
                new TPZSBFemElementGroup(*mphysics,index);
                geltogroup[el] = index;
            }
        }
//        mphysics->ApproxSpace().CreateCompEl(gel2D, *mphysics, index);
/// Generate the SBFem geometric elements
        TPZManVector<long> groupval;
        AddCollapsedVolumeElements(gmesh, groupval);
        
        gmesh->Print();
        
        std::cout << groupval << std::endl;
        
        mphysics->ApproxSpace().SetAllCreateFunctionsSBFem(2);
        
        mphysics->AutoBuild();
        
        nel = mphysics->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = mphysics->Element(el);
            if (!cel) {
                continue;
            }
            TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
            if (sbfem) {
                TPZGeoEl *gel = sbfem->Reference();
                long gelindex = gel->Index();
                TPZGeoElSide gelside(gel,4);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside && (neighbour.Element()->Dimension() != 1 || !neighbour.Element()->Reference())) {
                    neighbour = neighbour.Neighbour();
                }
                if (neighbour == gelside) {
                    DebugStop();
                }
                long skelindex = neighbour.Element()->Reference()->Index();
                sbfem->SetSkeleton(skelindex);
                
                long gelgroup = groupval[gelindex];
                long celgroupindex = geltogroup[gelgroup];
                
                TPZCompEl *celgr = mphysics->Element(celgroupindex);
                TPZSBFemElementGroup *sbfemgr = dynamic_cast<TPZSBFemElementGroup *>(celgr);
                if (!sbfemgr) {
                    DebugStop();
                }
                sbfemgr->AddElement(sbfem);
                sbfem->SetElementGroupIndex(celgroupindex);
                TPZElementMatrix E0,E1,E2,M0;

                sbfem->ComputeKMatrices(E0, E1, E2,M0);

                E0.fMat.Print("E0");
                E1.fMat.Print("E1");
                E2.fMat.Print("E2");
            }
        }
        
        int nrow = mphysics->NEquations();
        mphysics->Solution()(1,0) = 1.;
        mphysics->Solution()(2,0) = 1.;
//        for (int i=2; i<nrow-2; i+=2) {
//            mphysics->Solution()(i,0) = 1.;
//        }
        std::map<long,long>::iterator it;
        for(it = geltogroup.begin(); it != geltogroup.end(); it++)
        {
            long celindex = it->second;
            TPZCompEl *cel = mphysics->Element(celindex);
            TPZElementMatrix ek,ef;
            cel->CalcStiff(ek, ef);
            cel->LoadSolution();
        }

//        pyr->SetElasticity(100., 0.);
//        pyr->SetReferenceGeometric(gel2D->Index());
//
        // return the state variable
        int var = 0;
        TPZManVector<STATE,3> solout;
        for(it = geltogroup.begin(); it != geltogroup.end() && false; it++)
        {
            long celindex = it->second;
            TPZCompEl *cel = mphysics->Element(celindex);
            TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
            if (elgr)
            {
                TPZStack<TPZCompEl *, 5> els = elgr->GetElGroup();
                int nel = els.size();
                for (int el=0; el<nel; el++) {
                    TPZCompEl *cel = els[el];
                    TPZManVector<REAL,3> qsi(2,0.);
                    cel->Solution(qsi, var, solout);
                    TPZManVector<REAL,3> x(3);
                    cel->Reference()->X(qsi, x);
                    std::cout << "Solution at x " << x << " is " << solout << std::endl;
                }
            }
        }

        if(1)
        {
            std::ofstream out("../CompMultiphysicsbefore.txt");
            mphysics->Print(out);
        }
        
        
        
        // Visualization of computational meshes
        bool mustOptimizeBandwidth = false;
        TPZAnalysis * ElasticAnalysis = new TPZAnalysis(mphysics,mustOptimizeBandwidth);
        
        TPZStack<std::string> vecnames,scalname;
        scalname.Push("State");
//        scalname.Push("SigmaX");
//        scalname.Push("SigmaY");
//        scalname.Push("TauXY");

        ElasticAnalysis->DefineGraphMesh(2, scalname, vecnames, "../SBFem.vtk");
        int resolution = 0;
        ElasticAnalysis->PostProcess(resolution, 2);
        
        std::cout << "neq = " << mphysics->NEquations() << std::endl;
        SolveSist(ElasticAnalysis, mphysics);
        
        std::cout << "Post processing\n";
//        ElasticAnalysis->Solution().Print("Solution");
//        mphysics->Solution().Print("expandec");
        
        ElasticAnalysis->SetExact(Exact_Sol);
        
        TPZManVector<STATE> errors(3,0.);
        
        long neq_condense = ElasticAnalysis->Solution().Rows();
        long neq = mphysics->Solution().Rows();
        
        ElasticAnalysis->PostProcessError(errors);
        
        std::stringstream sout;
        sout << "../Errors" << PElasticity << ".txt";
        std::ofstream results(sout.str(),std::ios::app);
        results << nx[0] << " " << neq << " " << neq_condense << " " << errors << std::endl;
        
        gmesh->ResetReference();
        delete ElasticAnalysis;
        delete mphysics;
    }
    
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}


TPZCompMesh * ComputationalTensorMesh(TPZAutoPointer<TPZGeoMesh> & gmesh,int pOrder, int pInternal, bool hybrid)
{
    // Plane strain assumption
    int planestress = 0;
    
    // Getting mesh dimension
    int dim = 2;
    int matId1 = 1;
    
    TPZMatElasticity2D *material;
    material = new TPZMatElasticity2D(matId1);
    
    
    // Setting up paremeters
    material->SetfPlaneProblem(planestress);
    REAL lamelambda = 0.0e9,lamemu = 0.5e9, fx= 0, fy = 0;
    material->SetParameters(lamelambda,lamemu, fx, fy);
    //material->SetElasticParameters(40.0,0.0);
    REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
    material->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
    REAL Alpha = 1.0;
    //material->SetBiotAlpha(Alpha);cade o metodo?
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    if (hybrid) {
        TPZFNMatrix<4> xk(2,2,0),xc(2,2,0),xb(2,2,0),xf(2,1,0);
        TPZMat1dLin *matskel = new TPZMat1dLin(ESkeleton);
        matskel->SetMaterial(xk, xc, xb, xf);
        cmesh->InsertMaterialObject(matskel);
        
        TPZLagrangeMultiplier *lagrange1 = new TPZLagrangeMultiplier(EInterfaceMat1,1,2);
        cmesh->InsertMaterialObject(lagrange1);
        TPZLagrangeMultiplier *lagrange2 = new TPZLagrangeMultiplier(EInterfaceMat2,1,2);
        lagrange2->SetMultiplier(-1.);
        cmesh->InsertMaterialObject(lagrange2);
        
    }
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsSymTensor(dim);
    if (hybrid) {
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = material->CreateBC(material,Ebc1,1, val1, val2);
    
    val2(0,0) = 1.0*1000.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond3 = material->CreateBC(material,Ebc2,1, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond4 = material->CreateBC(material,Ebc3,1, val1, val2);
    
    val2(0,0) = -1.0*1000.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond5 = material->CreateBC(material,Ebc4,1, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BWrap = material->CreateBC(material,Ewrap,1, val1, val2);
    
    
    cmesh->InsertMaterialObject(material);
    // the other material
    material = dynamic_cast<TPZMatElasticity2D *>( material->NewMaterial());
    material->SetId(Emat2);
    material->SetElasticity(10., 0.);
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    cmesh->InsertMaterialObject(BWrap);
    
    std::set<int> matid;
    matid.insert(Emat1);
    matid.insert(Emat2);
    cmesh->AutoBuild(matid);
    
    
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    long nel = cmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != 2) {
            continue;
        }
        int side = gel->NSides()-1;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            DebugStop();
        }
        intel->SetSideOrder(side, pInternal);
    }
    
    cmesh->ExpandSolution();
    
//    cmesh->Print();
    
    cmesh->LoadReferences();
    
    if (hybrid) {
        TPZVec<TPZManVector<long,3> > elwrap(cmesh->NElements());
        // create a layer of 1D elements around each triangle
        // the computational element has to have the same connect indices
        long nel = cmesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if (!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() < 2) {
                continue;
            }
            elwrap[el].Resize(3, -1);
            int ncorner = gel->NCornerNodes();
            for (int side = ncorner; side < gel->NSides()-1; side++) {
                TPZStack<TPZCompElSide> celstack;
                TPZGeoElSide gelside(gel,side);
                elwrap[el][side-ncorner] = -1;
                gelside.EqualLevelCompElementList(celstack, false, false);
                for (int ic=0; ic<celstack.size(); ic++) {
                    TPZCompElSide celside = celstack[ic];
                    TPZGeoEl *gels = celside.Element()->Reference();
                    if (gels->Dimension() == 2) {
                        continue;
                    }
                    int matid = gels->MaterialId();
                    if (matid >= Ebc1 && matid <= Ebc4) {
                        elwrap[el][side-ncorner] = celside.Element()->Index();
                    }
                }
                if (elwrap[el][side-ncorner] == -1) {
                    // create the boundary element
                    TPZGeoEl *gwrap = gel->CreateBCGeoEl(side, Ewrap);
                    long index;
                    TPZCompEl *wrap = cmesh->CreateCompEl(gwrap, index);
                    int nsc = intel->NSideConnects(side);
                    for (int ic=0; ic<nsc; ic++) {
                        long cindex = intel->SideConnectIndex(ic, side);
                        wrap->SetConnectIndex(ic, cindex);
                    }
                    elwrap[el][side-ncorner] = wrap->Index();
                }
            }
//            std::cout << "el " << el << " elwrap " << elwrap[el] << std::endl;
        }
        
        cmesh->LoadReferences();
        
        TPZVec<TPZManVector<long,3> > pressureelements(cmesh->NElements());
        TPZVec<TPZManVector<long,3> > elwrapsign(cmesh->NElements());
        // create a pressure polynomial element element between the elements
        // set which element contributes positive and negative
        nel = cmesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if (!cel) continue;
            pressureelements[el].Resize(3, -1);
            elwrapsign[el].Resize(3, 1);
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            TPZGeoEl *gel = cel->Reference();
            int ncorner = gel->NCornerNodes();
            for (int side = ncorner; side < gel->NSides()-1; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                bool foundpressure = false;
                bool bcfound = false;
                while(neighbour != gelside) {
                    TPZGeoEl *gels = neighbour.Element();
                    if (gels->Dimension() == 2) {
                        neighbour = neighbour.Neighbour();
                        continue;
                    }
                    int matid = gels->MaterialId();
                    if (matid == ESkeleton) {
                        pressureelements[el][side-ncorner] = neighbour.Element()->Index();
                        elwrapsign[el][side-ncorner] = -1;
                        foundpressure = true;
                        break;
                    }
                    if (matid >= Ebc1 && matid <= Ebc4) {
                        bcfound = true;
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if (foundpressure == false && bcfound == false) {
                    // create the boundary element
                    TPZGeoEl *gdispl = gel->CreateBCGeoEl(side, ESkeleton);
//                    long index;
//                    TPZCompEl *displ = cmesh->CreateCompEl(gdispl, index);
//                    int nconn = displ->NConnects();
//                    for (int ic = 0; ic<nconn; ic++) {
//                        TPZConnect &c = displ->Connect(ic);
//                        c.SetLagrangeMultiplier(2);
//                    }
                    pressureelements[el][side-ncorner] = gdispl->Index();
                    elwrapsign[el][side-ncorner] = 1;
                }
            }
//            std::cout << "el " << el << " pressure elements " << pressureelements[el] << std::endl;

        }

        gmesh->ResetReference();
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);

        // transform the geometry indices into computational element indices
        std::set<int> matids;
        matids.insert(ESkeleton);
        cmesh->AutoBuild(matids);
        cmesh->LoadReferences();
        nel = pressureelements.size();
        for (long el=0; el<nel; el++) {
//            std::cout << "el " << el << " pressure indexes " << pressureelements[el] << std::endl;
            for (int is = 0; is<3; is++) {
                long elindex = pressureelements[el][is];
                if (elindex == -1) {
                    continue;
                }
                TPZGeoEl *gel = gmesh->Element(elindex);
                if (!gel) {
                    DebugStop();
                }
                TPZCompEl *cel = gel->Reference();
                if (!cel) {
                    DebugStop();
                }
                int nc = cel->NConnects();
                for (int ic=0; ic<nc; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    c.SetLagrangeMultiplier(2);
                }
                pressureelements[el][is] = cel->Index();
            }
//            std::cout << "el " << el << " pressure indexes " << pressureelements[el] << std::endl;
        }
        
//        cmesh->Print();
        
        TPZVec<TPZManVector<long,3> > interfaceelements(cmesh->NElements());
        // create an interface element between the wrap and the pressure element
        nel = cmesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if (!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() != 2) {
                continue;
            }
            interfaceelements[el].Resize(3, -1);
            int ncorner = gel->NCornerNodes();
            for (int side = ncorner; side < gel->NSides()-1; side++) {
                if (pressureelements[el][side-ncorner] == -1) {
                    continue;
                }
                long leftelindex = elwrap[el][side-ncorner];
                long dispindex = pressureelements[el][side-ncorner];
                TPZCompElSide compleft(cmesh->Element(leftelindex),2);
                TPZCompElSide compright(cmesh->Element(dispindex),2);
                long intindex;
                int matid = EInterfaceMat1;
                if (elwrapsign[el][side-ncorner] == -1) {
                    matid = EInterfaceMat2;
                }
                TPZGeoEl *ginterf = gel->CreateBCGeoEl(side, matid);
                new TPZInterfaceElement(*cmesh,ginterf,intindex,compleft,compright);
                interfaceelements[el][side-ncorner] = intindex;
            }
//            std::cout << "el " << el << " interface elements " << interfaceelements[el] << std::endl;
        }
        
        // create element groups of the elements, elwrap and pressure elements
        
        // create condensed elements
    }
    
    return cmesh;
    
}

TPZCompMesh * ComputationalSkeletonMesh(TPZAutoPointer<TPZGeoMesh> & gmesh,int pOrder)
{
    // Plane strain assumption
    int planestress = 0;
    
    // Getting mesh dimension
    int dim = 2;
    int matId1 = 1;
    
    TPZMaterial *material;
    int nstate = 1;
    bool elasticity = false;
    if (elasticity)
    {
        TPZMatElasticity2D *matloc = new TPZMatElasticity2D(matId1);
        material = matloc;
        nstate = 2;
        // Setting up paremeters
        matloc->SetfPlaneProblem(planestress);
        REAL lamelambda = 1.0e9,lamemu = 0.5e3, fx= 0, fy = 0;
        matloc->SetParameters(lamelambda,lamemu, fx, fy);
        //material->SetElasticParameters(40.0,0.0);
        REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
        matloc->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
        REAL Alpha = 1.0;
    }
    else
    {
        TPZMatLaplacian *matloc = new TPZMatLaplacian(matId1);
        matloc->SetDimension(2);
        matloc->SetSymmetric();
        material = matloc;
        nstate = 1;
    }
    //material->SetBiotAlpha(Alpha);cade o metodo?
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    TPZFMatrix<STATE> val1(nstate,nstate,0.), val2(nstate,1,0.);
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = material->CreateBC(material,Ebc1,1, val1, val2);
    
    val2(0,0) = 1.0*1000.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond3 = material->CreateBC(material,Ebc2,1, val1, val2);
    
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond4 = material->CreateBC(material,Ebc3,1, val1, val2);
    
    val2(0,0) = -1.0*1000.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BCond5 = material->CreateBC(material,Ebc4,1, val1, val2);
    
    val2(0,0) = 0.0;
    //    val2(1,0) = 0.0;
    TPZMaterial * BSkeleton = material->CreateBC(material,ESkeleton,1, val1, val2);
    
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    cmesh->InsertMaterialObject(BSkeleton);
    // the other material
    material =  material->NewMaterial();
    material->SetId(Emat2);
    //    material->SetElasticity(10., 0.);
    cmesh->InsertMaterialObject(material);
    cmesh->AutoBuild();
    return cmesh;
    
}

#define VTK
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh)
{
    TPZSkylineStructMatrix skymat(fCmesh);
    std::set<int> matids;
    matids.insert(Emat1);
    matids.insert(Emat2);
    matids.insert(Ebc1);
    matids.insert(Ebc2);
    matids.insert(Ebc3);
    matids.insert(Ebc4);
    matids.insert(EInterfaceMat1);
    matids.insert(EInterfaceMat2);
//    skymat.SetMaterialIds(matids);
    
    an->SetStructuralMatrix(skymat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    an->Assemble();
    an->Solve();
}

void PostProcessElasticity(TPZAnalysis &an, std::string plotfile)
{
    TPZManVector<std::string,10> scalnames(0), vecnames(0);
    
    
    scalnames.Resize(2);
    vecnames.Resize(1);
    scalnames[0] = "SigmaX";
    scalnames[1] = "SigmaY";
    vecnames[0]= "Displacement";
    
    const int dim = 2;
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
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

void MultiPhysicsMesh(TPZVec<TPZCompMesh *> &meshvec, TPZCompMesh *cmesh, bool hybrid)
{
    
    cmesh->SetDimModel(2);
    
    
    // Plane strain assumption
    int planestress = 0;
    
    // Getting mesh dimension
    int dim = 2;
    
    TPZMatElasticity2D *material;
    material = new TPZMatElasticity2D(Emat1);
    
//    if(y<0.5)
//    {
//        G = 0.5; material 1
//        Lambda = 0.;
//    }
//    else
//    {
//        G = 5.; material 2
//        Lambda = 0.;
//    }
    TPZGeoMesh *gmesh = cmesh->Reference();
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        int matid = gel->MaterialId();
        if (matid == Emat1) {
            TPZManVector<REAL,3> centerxi(2), centerx(3);
            gel->CenterPoint(gel->NSides()-1, centerxi);
            gel->X(centerxi, centerx);
            if (centerx[1] > 0.5) {
                gel->SetMaterialId(Emat2);
            }
        }
    }
    
    // Setting up paremeters
    material->SetfPlaneProblem(planestress);
    REAL lamelambda = 1.0e0,lamemu = 0.5e0, fx= 3., fy = 1.5;
    material->SetParameters(lamelambda,lamemu, fx, fy);
    material->SetElasticity(1., 0);
    TPZDummyFunction<STATE> *dummy2 = new TPZDummyFunction<STATE>(Force);
    dummy2->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > Force(dummy2);
    material->SetForcingFunction(Force);
    //material->SetElasticParameters(40.0,0.0);
    REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
    material->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
    REAL Alpha = 1.0;
    //material->SetBiotAlpha(Alpha);cade o metodo?
    
    cmesh->SetDimModel(dim);
    {
        TPZFNMatrix<4> xk(2,2,0),xc(2,2,0),xb(2,2,0),xf(2,1,0);
        TPZMat1dLin *matskel = new TPZMat1dLin(ESkeleton);
        matskel->SetMaterial(xk, xc, xb, xf);
        cmesh->InsertMaterialObject(matskel);
        
        int dimension = 2;
        TPZLagrangeMultiplier *lagrange1 = new TPZLagrangeMultiplier(EInterfaceMat1,dimension,2);
        cmesh->InsertMaterialObject(lagrange1);
        TPZLagrangeMultiplier *lagrange2 = new TPZLagrangeMultiplier(EInterfaceMat2,dimension,2);
        lagrange2->SetMultiplier(1.);
        cmesh->InsertMaterialObject(lagrange2);
        
    }
    

    
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Displacement);
    dummy->SetPolynomialOrder(1);
    TPZAutoPointer<TPZFunction<STATE> > localforce(dummy);
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = material->CreateBC(material,Ebc1,TPZMatElasticity2D::EDirichletXY, val1, val2);
    BCond2->SetForcingFunction(localforce);
    
    val2(0,0) = 1.0*1000.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond3 = material->CreateBC(material,Ebc2,TPZMatElasticity2D::EDirichletXY, val1, val2);
    BCond3->SetForcingFunction(localforce);
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond4 = material->CreateBC(material,Ebc3,TPZMatElasticity2D::EDirichletXY, val1, val2);
    BCond4->SetForcingFunction(localforce);
    val2(0,0) = -1.0*1000.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond5 = material->CreateBC(material,Ebc4,TPZMatElasticity2D::EDirichletXY, val1, val2);
    BCond5->SetForcingFunction(localforce);
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;

    TPZMaterial * BWrap = material->CreateBC(material,Ewrap,TPZMatElasticity2D::EDirichletXY, val1, val2);
    
    
    cmesh->InsertMaterialObject(material);
    // the other material
    material = dynamic_cast<TPZMatElasticity2D *>( material->NewMaterial());
    material->SetId(Emat2);
    material->SetElasticity(10., 0.);
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    cmesh->InsertMaterialObject(BWrap);
    // Transferindo para a multifisica
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmesh);
    
    cmesh->LoadReferences();
    
    if (hybrid) {
        long nel = cmesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) DebugStop();
            int matid = gel->MaterialId();
            // look for all elements of type diplacement (1D)
            if (matid != ESkeleton) {
                continue;
            }
            TPZStack<TPZGeoElSide> gelstack;
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZCompElSide celside(gelside.Reference());
            if(!celside) DebugStop();
            // loop over the neighbours, look for the wrap elements (there should be exactly 2)
            TPZGeoElSide neighbour(gelside.Neighbour());
            TPZGeoElSide gelinterface1,gelinterface2;
            while (neighbour != gelside) {
                if (neighbour.Element()->MaterialId() == Ewrap) {
                    gelstack.push_back(neighbour);
                }
                if (neighbour.Element()->MaterialId() == EInterfaceMat1) {
                    gelinterface1 = neighbour;
                }
                if (neighbour.Element()->MaterialId() == EInterfaceMat2) {
                    gelinterface2 = neighbour;
                }
                neighbour = neighbour.Neighbour();
            }
            if (gelstack.size() != 2 || !gelinterface1 || !gelinterface2 || !gelstack[0].Reference() || !gelstack[1].Reference()) {
                DebugStop();
            }
//            /** @brief Constructor */
//            TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, long &index, TPZCompElSide left, TPZCompElSide right);

            // create the two interface elements
            long index;
            new TPZMultiphysicsInterfaceElement(*cmesh,gelinterface1.Element(), index, gelstack[0].Reference(), celside);
            new TPZMultiphysicsInterfaceElement(*cmesh,gelinterface2.Element(), index, gelstack[1].Reference(), celside);
        }
    }
    
}

// Condense the internal degrees of freedom
void CondenseInternalEquations(TPZCompMesh *mphysics)
{
    int numskeleton = 0;
    // loop over all elements
    mphysics->LoadReferences();
    // group the elements with the boundary elements
    long nel = mphysics->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = mphysics->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() == 1) {
            continue;
        }
        std::set<long> celconnects;
        cel->BuildConnectList(celconnects);
        int ncorner = gel->NCornerNodes();
        int nsides = gel->NSides();
        // Add the connects of all neighbouring skeleton elements
        for (int is=ncorner; is<nsides-1; is++) {
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gel,is);
            gelside.ConnectedCompElementList(celstack, 0, 0);
            for (int ist=0; ist<celstack.size(); ist++) {
                TPZCompEl *cel = celstack[ist].Element();
                int matid = cel->Reference()->MaterialId();
                if (matid == ESkeleton) {
                    numskeleton++;
                    cel->BuildConnectList(celconnects);
                }
            }
        }
        
//        std::cout << "connects associated with element " << el << " | ";
//        for(int const &val : celconnects)
//        {
//            std::cout << val << " ";
//        }
//        std::cout << std::endl;
        
        TPZStack<TPZCompEl *> celstack;
        for (int is=gel->NCornerNodes(); is<nsides-1; is++) {
            if (gel->SideDimension(is) != 1) {
                DebugStop();
            }
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> celstackinternal;
            gelside.ConnectedCompElementList(celstackinternal, 0, 0);
            int nneigh = celstackinternal.size();
            for (int in=0; in<nneigh; in++) {
                std::set<long> neighconnect;
                TPZCompEl *neighcompel = celstackinternal[in].Element();
                int matid = neighcompel->Reference()->MaterialId();
                if (matid == ESkeleton) {
                    continue;
                }
                neighcompel->BuildConnectList(neighconnect);
                TPZManVector<long> intersect(neighconnect.size(),-1);
                if(neighconnect.size()==0) continue;
                long *output;
                output = std::set_intersection(celconnects.begin(), celconnects.end(), neighconnect.begin(), neighconnect.end(), &intersect[0]);
                if (output - &intersect[0] == neighconnect.size()) {
                    celstack.push_back(celstackinternal[in].Element());
                }
            }
        }
        if (celstack.size()) {
            long index;
            TPZElementGroup *elgr = new TPZElementGroup(*mphysics,index);
            elgr->AddElement(cel);
            long nstack = celstack.size();
            for (int ist=0; ist<nstack; ist++) {
                TPZCompEl *celst = celstack[ist];
                elgr->AddElement(celst);
            }
//            elgr->Print();
        }
    }
    // compute the number of elements connected to each node
    mphysics->ComputeNodElCon();
    
    if (numskeleton == 0)
    {
        // increment the number of connected elements for two of the elasticity nodes
        nel = mphysics->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = mphysics->Element(el);
            if (!cel) {
                continue;
            }
            int nc = cel->NConnects();
            int count = 0;
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                int lagrange = c.LagrangeMultiplier();
                if (lagrange == 0 || c.NElConnected() != 1) {
                    continue;
                }
                count++;
                c.IncrementElConnected();
                if (count == 2) {
                    break;
                }
            }
        }
    }
    // condense the internal nodes
    nel = mphysics->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = mphysics->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel && gel->MaterialId() == ESkeleton) {
            continue;
        }
        TPZCondensedCompEl *celcond = new TPZCondensedCompEl(cel);
    }
    // just to be sure
    mphysics->ComputeNodElCon();
}

/// create the 1d skeleton geometric elements
void AddSkeletonElements(TPZAutoPointer<TPZGeoMesh> gmesh)
{
    // create a lower dimension element on each boundary
    int dim = gmesh->Dimension();
    
    long nel = gmesh->NElements();
    
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide thisside(gel,is);
            TPZGeoElSide neighbour = thisside.Neighbour();
            while (neighbour != thisside && neighbour.Element()->Dimension() != dim-1) {
                neighbour = neighbour.Neighbour();
            }
            // if we didnt find a lower dimension element
            if (thisside == neighbour) {
                gel->CreateBCGeoEl(is,ESkeleton);
            }
        }
    }
}

/// create the collapsed volume elements
/// groupval : geometric element index to which a collapsed element belongs
void AddCollapsedVolumeElements(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<long> &groupval)
{
    int dim = gmesh->Dimension();
    long nel = gmesh->NElements();
    long ncreated = 0;
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->MaterialId() != EGroup || gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is=gel->NCornerNodes(); is<nsides-1; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> connected;
            gelside.EqualorHigherCompElementList2(connected,0,0);
            ncreated += connected.size();
        }
    }
    groupval.resize(nel+ncreated);
    groupval.Fill(-1);
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->MaterialId() != EGroup || gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        TPZManVector<REAL,3> xicenter(dim),xcenter(3);
        gel->CenterPoint(nsides-1,xicenter);
        gel->X(xicenter,xcenter);
        long middlenode = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[middlenode].Initialize(xcenter,gmesh);
        for (int is=gel->NCornerNodes(); is<nsides-1; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> connected;
            TPZManVector<long,20> nodeindexes;
            gelside.EqualorHigherCompElementList2(connected,0,0);
            for (int sbel=0; sbel<connected.size(); sbel++) {
                TPZGeoEl *sbgel = connected[sbel].Element()->Reference();
                long ncorner = sbgel->NCornerNodes();
                nodeindexes.resize(2*ncorner);
                for (int in=0; in < ncorner; in++) {
                    nodeindexes[in] = sbgel->NodeIndex(in);
                }
                for (int in=ncorner; in<2*ncorner; in++) {
                    nodeindexes[in] = middlenode;
                }
                switch (sbgel->Type()) {
                    case EOned:
                        long index;
                        gmesh->CreateGeoElement(EQuadrilateral, nodeindexes, Emat1, index);
                        groupval[index] = el;
                        break;
                        
                    default:
                        DebugStop();
                        break;
                }
            }
            ncreated += connected.size();
        }
    }
    gmesh->BuildConnectivity();
}

