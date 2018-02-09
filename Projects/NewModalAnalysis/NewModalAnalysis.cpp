/**
 * @file
 * @brief Performs modal analysis in an electromagnetic waveguide.
 * @details This project is aimed at electromagnetic analysis
 * of waveguides, but it is not restricted to it. More generally,
 * it could be a starting point for problems with longitudinal axis symmetry.
 * The problem is stated as a generalised eigenvalue problem, as seen in
 * Jin, Jian-Ming. The finite element method in electromagnetics.(contd.)
 * John Wiley & Sons, 2015 (Chapter 8).
 * The transverse field components are approximated using H(curl,\Omega)
 * functions,
 * and the longitudinal field components are approximated using H^1(\Omega)
 * functions.
 *
 * @author Francisco Orlandini
 * @since 2015
 */

#include <stddef.h>               // for NULL
#include <fstream>
#include <TPZEigenAnalysis.h>
#include <TPZSpStructMatrix.h>
#include <tpzgeoelmapped.h>
#include <tpzarc3d.h>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "pzstrmatrix.h"
#include "pzl2projection.h"
#include "pzfstrmatrix.h"
#include "pzbuildmultiphysicsmesh.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif
#ifdef USING_SLEPC
#include <TPZSlepcEPSHandler.h>
#include <TPZSlepcSTHandler.h>
#endif
#include "parameter_handler.h"
#include "TPZMatWaveguideCutOffAnalysis.h"
#include "TPZMatModalAnalysis.h"
#include "TPZGmshReader.h"
#include "SPZModalAnalysisDataReader.h"
#include "SPZModalAnalysisData.h"

void RunSimulation(SPZModalAnalysisData &simData);

void ReadGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec);

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh,
                 int pOrder, TPZVec<int> &matIdVec, TPZVec<STATE> &urVec,
                 TPZVec<STATE> &er, REAL f0, bool isCutOff);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> cmeshMF,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif



    ParameterHandler prm;
    SPZModalAnalysisDataReader reader(prm,argc,argv);
    SPZModalAnalysisData simData;
    reader.ReadParameters(simData);

    RunSimulation(simData);

    return 0;
}

void RunSimulation(SPZModalAnalysisData &simData) {
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    std::cout<<"Creating GMesh...";
    #ifdef USING_BOOST
    boost::posix_time::ptime t1_g =
        boost::posix_time::microsec_clock::local_time();
    #endif
    TPZVec<int> matIdVec;
    ReadGMesh(gmesh, simData.physicalOpts.meshFile, matIdVec);
    #ifdef USING_BOOST
    boost::posix_time::ptime t2_g =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created!  "<<t2_g-t1_g<<std::endl;
    #endif
    TPZVec<TPZCompMesh *> meshVec(1);
    std::cout<<"Creating CMesh...";
    #ifdef USING_BOOST
    boost::posix_time::ptime t1_c =
        boost::posix_time::microsec_clock::local_time();
    #endif

    CreateCMesh(meshVec, gmesh, simData.pzOpts.pOrder, matIdVec,
                simData.physicalOpts.urVec, simData.physicalOpts.erVec,
                simData.physicalOpts.fOp,simData.physicalOpts.isCutOff); // funcao para criar a malha computacional

    #ifdef USING_BOOST
    boost::posix_time::ptime t2_c =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_c-t1_c<<std::endl;
    #endif
    TPZCompMesh *cmesh = meshVec[0];

    TPZEigenAnalysis an(cmesh, true);

    TPZManVector<long, 1000> activeEquations;
    int neq = 0;
    int neqOriginal = 0;

    TPZAutoPointer<TPZStructMatrix> strmtrx;
    strmtrx = new TPZSpStructMatrix(cmesh);
    strmtrx->SetNumThreads(simData.pzOpts.nThreads);
    FilterBoundaryEquations(meshVec, activeEquations, neq, neqOriginal);
    strmtrx->EquationFilter().SetActiveEquations(activeEquations);
    an.SetStructuralMatrix(strmtrx);

    //const int nSolutions = neq >= simData.solverOpts.eps_nev ? simData.solverOpts.eps_nev : neq;
    TPZSlepcEPSHandler<STATE> solver;
    TPZSlepcSTHandler stHandler;
    solver.SetTolerances(simData.solverOpts.eps_tol,simData.solverOpts.eps_max_its);
    solver.SetConvergenceTest(simData.solverOpts.eps_conv_test);
    solver.SetWhichEigenpairs(simData.solverOpts.eps_which_eig);
    solver.SetTargetEigenvalue(simData.solverOpts.target);

    stHandler.SetPrecond(simData.solverOpts.st_precond);
    stHandler.SetSolver(simData.solverOpts.st_solver);
    stHandler.SetSolverTol(simData.solverOpts.ksp_rtol,simData.solverOpts.ksp_atol,
                           simData.solverOpts.ksp_dtol,simData.solverOpts.ksp_max_its);
    stHandler.SetType(simData.solverOpts.st_type,simData.solverOpts.target);
    solver.SetST(stHandler);

    solver.SetTrueResidual(simData.solverOpts.eps_true_res);
    solver.SetProblemType(simData.solverOpts.eps_prob_type);
    solver.SetType(simData.solverOpts.eps_type);
    solver.SetKrylovOptions(simData.solverOpts.eps_krylov_locking,simData.solverOpts.eps_krylov_restart);
    solver.SetEPSDimensions(simData.solverOpts.eps_nev, simData.solverOpts.eps_ncv, simData.solverOpts.eps_mpd);
    solver.SetVerbose(simData.solverOpts.eps_verbose);
    an.SetSolver(solver);

    std::cout << "Assembling..." << std::endl;
#ifdef USING_BOOST
    boost::posix_time::ptime t1 =
        boost::posix_time::microsec_clock::local_time();
#endif
    an.Assemble();
#ifdef USING_BOOST
    boost::posix_time::ptime t2 =
        boost::posix_time::microsec_clock::local_time();
#endif
    std::cout << "Finished assembly." << std::endl;

    std::cout << "Solving..." << std::endl;
#ifdef USING_BOOST
    boost::posix_time::ptime t3 =
        boost::posix_time::microsec_clock::local_time();
#endif

    an.Solve();
#ifdef USING_BOOST
    boost::posix_time::ptime t4 =
        boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2 - t1 << " Time for solving "
              << t4 - t3 << std::endl;
#endif
  std::cout << "FINISHED!" << std::endl;
  return;
}

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> meshVec,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal) {
    TPZCompMesh *cmesh = meshVec[0];

    TPZManVector<long, 1000> allConnects;
    std::set<long> boundConnects;

    for (int iel = 0; iel < cmesh->NElements(); iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if (cel == NULL) {
            continue;
        }
        if (cel->Reference() == NULL) {

            continue;
        }
        if (cel->Reference()->MaterialId() == -1) {
            std::set<long> boundConnectsEl;
            std::set<long> depBoundConnectsEl;
            std::set<long> indepBoundConnectsEl;
            cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
            cel->BuildConnectList(boundConnectsEl);
            for (std::set<long>::iterator iT = boundConnectsEl.begin();
                 iT != boundConnectsEl.end(); iT++) {
                const long val = *iT;
                if (boundConnects.find(val) == boundConnects.end()) {
                    boundConnects.insert(val);
                }
            }
        }
    }
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        if (boundConnects.find(iCon) == boundConnects.end()) {
            TPZConnect &con = cmesh->ConnectVec()[iCon];
            if(con.HasDependency())continue;
            int seqnum = con.SequenceNumber();
            int pos = cmesh->Block().Position(seqnum);
            int blocksize = cmesh->Block().Size(seqnum);
            if (blocksize == 0)
                continue;

            int vs = activeEquations.size();
            activeEquations.Resize(vs + blocksize);
            for (int ieq = 0; ieq < blocksize; ieq++) {
                activeEquations[vs + ieq] = pos + ieq;
            }
        }
    }

    neqOriginal = cmesh->NEquations();
    neq = 0;
    TPZCompMesh *cmeshHCurl = meshVec[1];
    TPZCompMesh *cmeshH1 = meshVec[2];
    int nHCurlEquations = 0, nH1Equations = 0;
    TPZMatModalAnalysis *mat =
        dynamic_cast<TPZMatModalAnalysis *>(cmesh->FindMaterial(1));
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        bool isH1;
        if (boundConnects.find(iCon) == boundConnects.end()) {
            if (cmesh->ConnectVec()[iCon].HasDependency())
                continue;
            int seqnum = cmesh->ConnectVec()[iCon].SequenceNumber();
            int blocksize = cmesh->Block().Size(seqnum);
            if (mat->H1Index() == 0 && iCon < cmeshH1->NConnects()) {
                isH1 = true;
            } else if (mat->H1Index() == 1 && iCon >= cmeshHCurl->NConnects()) {
                isH1 = true;
            } else {
                isH1 = false;
            }
            for (int ieq = 0; ieq < blocksize; ieq++) {
                neq++;
                isH1 == true ? nH1Equations++ : nHCurlEquations++;
            }
        }
    }
    std::cout << "------\tactive eqs\t-------" << std::endl;
    std::cout << "# H1 equations: " << nH1Equations << std::endl;
    std::cout << "# HCurl equations: " << nHCurlEquations << std::endl;
    std::cout << "# equations: " << neq << std::endl;
    std::cout << "------\t----------\t-------" << std::endl;
    return;
}



void ReadGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec) {
    TPZGmshReader meshReader;
    gmesh = meshReader.GeometricGmshMesh(mshFileName);
#ifdef PZDEBUG
//	TPZCheckGeom * Geometrytest = new TPZCheckGeom(gmesh);
//	int isBadMeshQ = Geometrytest->PerformCheck();
//
//	if (isBadMeshQ) {
//		DebugStop();
//	}
#endif
    auto matIds = meshReader.fMaterialDataVec.fMatID;
    matIdVec.Resize(matIds.size());
    int i = 0;
    for(auto id = matIds.begin(); id != matIds.end(); id++,i++ )    {
      matIdVec[i] = *id;
    }

    std::string meshFileName(mshFileName);
    const size_t strlen = meshFileName.length();
    meshFileName.append(".vtk");
    std::ofstream outVTK(meshFileName.c_str());
    meshFileName.replace(strlen, 4, ".txt");
    std::ofstream outTXT(meshFileName.c_str());

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
    gmesh->Print(outTXT);
    outTXT.close();
    outVTK.close();

    return;
}

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh,
                 int pOrder, TPZVec<int> &matIdVec, TPZVec<STATE> &urVec,
                 TPZVec<STATE> &erVec, REAL f0, bool isCutOff) {
    enum {
      dirichlet = 0
    }; // tipo da condicao de contorno do problema
    const int dim = 2;   // dimensao do problema

    TPZVec<int> volMatId(matIdVec.size()-1);
    for(int i = 0; i< volMatId.size(); i++){
      volMatId[i] = matIdVec[i];
    }
    const int boundMatId = matIdVec[matIdVec.size()-1];
    if(volMatId.size()!=urVec.size()) DebugStop();

    /// criar malha computacional H1

    TPZCompMesh *cmeshH1 = new TPZCompMesh(gmesh);
    cmeshH1->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmeshH1->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha
    const int nState = 1;
    TPZVec<STATE> sol; // only for creating material. this material will not be
    // used in reality

    TPZL2Projection *matH1;
    for (int i = 0; i < volMatId.size(); ++i) {
        matH1 = new TPZL2Projection(volMatId[i], dim, nState, sol);
        cmeshH1->InsertMaterialObject(matH1);
    }

    /// electrical conductor boundary conditions
    TPZFNMatrix<1, STATE> val1(1, 1, 0.), val2(1, 1, 0.);

    TPZMaterial *BCondH1Dir = matH1->CreateBC(
        matH1, boundMatId, dirichlet, val1, val2); // cria material que implementa a
    // condicao de contorno de dirichlet

    cmeshH1->InsertMaterialObject(BCondH1Dir); // insere material na malha
    // Cria elementos computacionais que gerenciarao o espaco de aproximacao da
    // malha
    cmeshH1->AutoBuild();
    cmeshH1->CleanUpUnconnectedNodes();

    /// criar malha computacional HCurl

    TPZCompMesh *cmeshHCurl = new TPZCompMesh(gmesh);
    cmeshHCurl->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmeshHCurl->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha

    TPZVecL2 *matHCurl;
    for (int i = 0; i < volMatId.size(); ++i) {
        matHCurl = new TPZVecL2(volMatId[i]);
        cmeshHCurl->InsertMaterialObject(matHCurl);
    }

    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *BCondHCurlDir =
        matHCurl->CreateBC(matHCurl, boundMatId, dirichlet, val1,
                           val2); // cria material que implementa a condicao de
    // contorno de dirichlet

    cmeshHCurl->InsertMaterialObject(BCondHCurlDir); // insere material na malha

    cmeshHCurl->SetAllCreateFunctionsHCurl(); // define espaco de aproximacao

    cmeshHCurl->AutoBuild();
    cmeshHCurl->CleanUpUnconnectedNodes();

    TPZCompMesh *cmeshMF = new TPZCompMesh(gmesh);
    TPZVec<TPZMatModalAnalysis *>matMultiPhysics(volMatId.size());
    if (isCutOff) {
        for (int i = 0; i < volMatId.size(); ++i) {
            matMultiPhysics[i] =
                new TPZMatWaveguideCutOffAnalysis(volMatId[i], f0, urVec[i], erVec[i]);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    } else {
        for (int i = 0; i < volMatId.size(); ++i) {
            matMultiPhysics[i] =
                new TPZMatModalAnalysis(volMatId[i], f0, urVec[i], erVec[i]);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    }

    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *BCondMFDir =
        matMultiPhysics[0]->CreateBC(matMultiPhysics[0], boundMatId, dirichlet, val1,
                                     val2); // any material is fine

    cmeshMF->InsertMaterialObject(BCondMFDir); // insere material na malha
    std::set<int> set;
    for (int i = 0; i < volMatId.size(); ++i) {
        set.insert(volMatId[i]);
    }
    set.insert(boundMatId);

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild(set);
    cmeshMF->CleanUpUnconnectedNodes();

    TPZVec<TPZCompMesh *> meshVec(2);
    meshVec[matMultiPhysics[0]->H1Index()] = cmeshH1;
    meshVec[matMultiPhysics[0]->HCurlIndex()] = cmeshHCurl;

    TPZBuildMultiphysicsMesh::AddElements(meshVec, cmeshMF);
    TPZBuildMultiphysicsMesh::AddConnects(meshVec, cmeshMF);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVec, cmeshMF);
    cmeshMF->CleanUpUnconnectedNodes();

    cmeshMF->ExpandSolution();
    cmeshMF->CleanUpUnconnectedNodes();
    cmeshMF->ComputeNodElCon();
    cmeshMF->CleanUpUnconnectedNodes();
    std::ofstream fileH1("cmeshH1.txt");
    cmeshH1->Print(fileH1);
    std::ofstream fileHCurl("cmeshHCurl.txt");
    cmeshHCurl->Print(fileHCurl);
    std::ofstream fileMF("cmeshMFHCurl.txt");
    cmeshMF->Print(fileMF);

    meshVecOut.resize(3);

    meshVecOut[0] = cmeshMF;
    meshVecOut[1 + matMultiPhysics[0]->H1Index()] = cmeshH1;
    meshVecOut[1 + matMultiPhysics[0]->HCurlIndex()] = cmeshHCurl;
    return;
}
