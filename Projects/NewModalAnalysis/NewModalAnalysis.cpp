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
#include "boost/date_time/posix_time/posix_time.hpp"
#ifdef USING_SLEPC
#include <TPZSlepcEPSHandler.h>
#include <TPZSlepcSTHandler.h>
#include <Complex/TPZMatWaveguidePml.h>

#endif
#include "parameter_handler.h"
#include "TPZMatWaveguideCutOffAnalysis.h"
#include "TPZMatModalAnalysis.h"
#include "pzintel.h"
#include "TPZGmshReader.h"
#include "SPZModalAnalysisDataReader.h"
#include "SPZModalAnalysisData.h"

void RunSimulation(SPZModalAnalysisData &simData,std::ostringstream &eigeninfo);

void
ReadGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix, const bool &print, const REAL &scale = 1.);

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh,
                 int pOrder, TPZVec<int> &matIdVec, TPZVec<STATE> &urVec,
                 TPZVec<STATE> &erVec, REAL f0, bool isCutOff, const std::string &prefix,
                 const bool &print, const REAL &scale, const bool &hasPML, const REAL &alphaMax);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> cmeshMF,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);


void CreateGmshMesh(const std::string &meshName, const std::string &newName,
                    const REAL &factor, const int &nThreads,
                    const REAL &scale, const int &meshOrder);

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif



    ParameterHandler prm;
    SPZModalAnalysisDataReader reader(prm,argc,argv);
    SPZModalAnalysisData simData;
    reader.ReadParameters(simData);
    //freq sweep
    int nSteps = 1;
    REAL firstLambda = simData.physicalOpts.lambda;
    REAL lastLambda = firstLambda;
    REAL stepSize = 0;

    if(simData.pzOpts.freqSweep){
        firstLambda = simData.pzOpts.lambdaMin;
        lastLambda = simData.pzOpts.lambdaMax;
        stepSize = (lastLambda-firstLambda)/nSteps;
    }
    else{
        stepSize = 0;
        simData.pzOpts.freqSteps = 1;
    }
    std::string meshOriginal = simData.pzOpts.meshFile;
    const int pOrderOrig = simData.pzOpts.pOrder;
    std::ostringstream eigeninfo;

    boost::posix_time::ptime t1_total =
            boost::posix_time::microsec_clock::local_time();

    for (int i = 0; i < simData.pzOpts.hSteps; ++i) {
        std::cout<<"Beginning step "<<i+1<<" out of "<<simData.pzOpts.hSteps<<"h steps."<<std::endl;
        const REAL factorVal = simData.pzOpts.factorVec[i];
        simData.pzOpts.meshFile = meshOriginal.substr(0,meshOriginal.size()-4)
                                  +"ord"+std::to_string(simData.pzOpts.meshOrder)
                                  +"h"+std::to_string(i)
                                  +".msh";
        CreateGmshMesh(meshOriginal, simData.pzOpts.meshFile,
                       factorVal, simData.pzOpts.nThreads,
                       simData.pzOpts.scaleFactor,simData.pzOpts.meshOrder);
        simData.pzOpts.pOrder = pOrderOrig;
        for (int j = 0; j < simData.pzOpts.pSteps; ++j) {
            std::cout<<"h step: "<< i+1<<". Beginning step "<<j+1<<" out of "<<simData.pzOpts.hSteps<<"p steps."<<std::endl;
            for (int k = 0; k < simData.pzOpts.freqSteps; ++k) {
                if(simData.pzOpts.freqSteps > 1){
                    std::cout<<"h step: "<< i+1<<". p step: "<<j+1;
                    std::cout<<". Beginning step "<<k+1<<" out of "<<simData.pzOpts.freqSteps<<"freq steps."<<std::endl;
                }
                simData.physicalOpts.lambda = firstLambda + k * stepSize;
                RunSimulation(simData,eigeninfo);
            }
            simData.pzOpts.pOrder++;
        }
    }
    boost::posix_time::ptime t2_total =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Total time: "<<t2_total-t1_total<<std::endl;
    if(simData.pzOpts.exportEigen){
        std::cout<<"Exporting results..."<<std::endl;
        std::string eigenFileName = simData.pzOpts.prefix;
        eigenFileName +="mapord_"+std::to_string(simData.pzOpts.meshOrder)+"_";
        if(simData.pzOpts.freqSweep){
            eigenFileName+="from_l_"+std::to_string(firstLambda)+"_to_"+std::to_string(lastLambda)+"_";
        }else{
            eigenFileName+="l_"+std::to_string(firstLambda)+"_";
        }
        if(simData.pzOpts.pSteps>1){
            eigenFileName+="from_p_"+std::to_string(pOrderOrig)+"_to_";
            eigenFileName+=std::to_string(pOrderOrig+simData.pzOpts.pSteps-1)+"_";
        }else{
            eigenFileName+="p_"+std::to_string(pOrderOrig);
        }
        eigenFileName+="_n_hsteps_"+std::to_string(simData.pzOpts.hSteps)+".csv";

        std::ofstream file(eigenFileName.c_str());
        file << eigeninfo.str();
    }
    return 0;
}

void CreateGmshMesh(const std::string &meshName, const std::string &newName,
                    const REAL &factor, const int &nThreads,
                    const REAL &scale, const int &meshOrder){
    std::ostringstream str_factor;
    str_factor << std::setprecision(20) << factor;
    std::ostringstream str_scale;
    str_scale<< std::setprecision(20) << scale;

    std::string command = "gmsh " + meshName + " -2 -match ";
    command += " -nt " + std::to_string(nThreads);
    command += " -tol 1e-20 ";
    command += " -v 2 ";
    command += " -setnumber scale "+str_scale.str();
    command += " -setnumber factor "+str_factor.str();
    command += " -order " + std::to_string(meshOrder);
    if( meshOrder > 1 ) command += " -optimize_ho";
    command += " -o " + newName;
    std::cout<<"Generating mesh with: "<<std::endl<<command<<std::endl;

    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), 128, pipe.get()) != nullptr){
        result += buffer.data();
    }
    std::cout<<result<<std::endl;
}

void RunSimulation(SPZModalAnalysisData &simData,std::ostringstream &eigeninfo) {
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    std::cout<<"Creating GMesh..."<<std::endl;
    boost::posix_time::ptime t1_g =
        boost::posix_time::microsec_clock::local_time();
    TPZVec<int> matIdVec;
    ReadGMesh(gmesh, simData.pzOpts.meshFile, matIdVec, simData.pzOpts.prefix,
              simData.pzOpts.exportGMesh);
    boost::posix_time::ptime t2_g =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created!  "<<t2_g-t1_g<<std::endl;
    TPZVec<TPZCompMesh *> meshVec(1);
    std::cout<<"Creating CMesh..."<<std::endl;
    boost::posix_time::ptime t1_c =
        boost::posix_time::microsec_clock::local_time();

    CreateCMesh(meshVec, gmesh, simData.pzOpts.pOrder, matIdVec,
                simData.physicalOpts.urVec, simData.physicalOpts.erVec,
                simData.physicalOpts.lambda,simData.physicalOpts.isCutOff,
                simData.pzOpts.prefix,simData.pzOpts.exportCMesh,
                simData.pzOpts.scaleFactor,
                simData.physicalOpts.hasPML, simData.physicalOpts.alphaMax); // funcao para criar a malha computacional
    boost::posix_time::ptime t2_c =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_c-t1_c<<std::endl;
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
    solver.SetAbsoluteValue(simData.pzOpts.absVal);
    an.SetSolver(solver);

    if(simData.pzOpts.exportGMesh){
        TPZStack<std::string> scalnames, vecnames;
        vecnames.Push("Material");
        std::string plotfile = simData.pzOpts.prefix + "Material" + ".vtk";
        // estara na pasta debug
        const int dim = 2;
        an.DefineGraphMesh(dim, scalnames, vecnames,
                           plotfile);  // define malha grafica
        an.PostProcess(simData.pzOpts.vtkRes);
    }
    std::cout << "Assembling..." << std::endl;
    boost::posix_time::ptime t1 =
        boost::posix_time::microsec_clock::local_time();
    an.Assemble();
    boost::posix_time::ptime t2 =
        boost::posix_time::microsec_clock::local_time();
    std::cout << "Finished assembly." << std::endl;

    std::cout << "Solving..." << std::endl;
    boost::posix_time::ptime t3 =
        boost::posix_time::microsec_clock::local_time();

    an.Solve();
    boost::posix_time::ptime t4 =
        boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2 - t1 << " Time for solving "
              << t4 - t3 << std::endl;
    const TPZManVector<SPZAlwaysComplex<STATE>::type> eigenValues = an.GetEigenvalues();
    const TPZFMatrix<SPZAlwaysComplex<STATE>::type> eigenVectors = an.GetEigenvectors();

    typedef std::numeric_limits< double > dbl;
    std::cout.precision(dbl::max_digits10);
    for(int i = 0; i < eigenValues.size() ; i++){
        std::cout<<std::fixed<<eigenValues[i]<<std::endl;
    }
    if(simData.pzOpts.exportEigen){
        eigeninfo.precision(dbl::max_digits10);
        std::cout<<"Exporting eigen info..."<<std::endl;
        REAL hSize = 1e12;
        REAL tol = 0.000001;
        REAL elRadius = 0;
        TPZVec<REAL> qsi(2,0.25);
        TPZVec<REAL> x(3,0.);
        for (int j = 0; j < gmesh->NElements(); ++j) {
            TPZGeoEl &el = *(gmesh->ElementVec()[j]);
            for (int i = 0; i < el.NCornerNodes() ; ++i) {
                auto node = el.Node(i);
                if(node.Coord(0) < tol && node.Coord(1) < tol){
                    elRadius = el.ElementRadius();
                    hSize = elRadius < hSize ? elRadius : hSize;
                }
            }
//            TPZBndCond *mat = dynamic_cast<TPZBndCond *>(meshVec[0]->MaterialVec()[el.MaterialId()]);
//            if(!mat){
//                elRadius = el.ElementRadius();
//                hSize = elRadius < hSize ? elRadius : hSize;
//            }
        }
        eigeninfo << std::fixed << hSize << "," << simData.pzOpts.pOrder<<",";
        eigeninfo << std::fixed << simData.physicalOpts.lambda<<",";
        eigeninfo << eigenValues.size()<<",";
        for(int i = 0; i < eigenValues.size() ; i++){
            eigeninfo<<std::fixed<<std::real(eigenValues[i])<<",";
            eigeninfo<<std::fixed<<std::imag(eigenValues[i]);
            if(i == eigenValues.size() - 1 ) {
                eigeninfo << std::endl;
            }else{
                eigeninfo << ",";
            }
        }

        std::cout<<" Done!"<<std::endl;
    }

    if (simData.pzOpts.genVTK) {
        TPZMatModalAnalysis *matPointer =
                dynamic_cast<TPZMatModalAnalysis *>(meshVec[0]->MaterialVec()[1]);
        TPZVec<TPZCompMesh *> temporalMeshVec(2);
        temporalMeshVec[matPointer->H1Index()] = meshVec[1 + matPointer->H1Index()];
        temporalMeshVec[matPointer->HCurlIndex()] =
                meshVec[1 + matPointer->HCurlIndex()];

        std::cout << "Post Processing..." << std::endl;

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Ez");
        vecnames.Push("Et");
        std::string plotfile = simData.pzOpts.prefix + "fieldPlot" + ".vtk";
                                                        // estara na pasta debug
        const int dim = 2;
        an.DefineGraphMesh(dim, scalnames, vecnames,
                           plotfile);  // define malha grafica

        TPZFMatrix<SPZAlwaysComplex<STATE>::type> currentEigenvector(neq,1);
        TPZFMatrix<SPZAlwaysComplex<STATE>::type> scatteredEigen(neqOriginal,1);
        for (int iSol = 0; iSol < eigenValues.size(); iSol++) {
            for(int j = 0; j < eigenVectors.Rows(); j++){
                currentEigenvector(j,0) = simData.pzOpts.absVal ?
                                          std::abs(eigenVectors.GetVal(j,iSol)) :
                                          std::real(eigenVectors.GetVal(j,iSol));
            }
            strmtrx->EquationFilter().Scatter(currentEigenvector, scatteredEigen);
            an.LoadSolution(scatteredEigen);
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(temporalMeshVec,
                                                               cmesh);
            an.PostProcess(simData.pzOpts.vtkRes);
        }
    }
    gmesh->SetReference(nullptr);
    for (int k = 0; k < meshVec.size(); ++k) {
        meshVec[k]->SetReference(nullptr);
        delete meshVec[k];
        meshVec[k] = nullptr;
    }
    delete gmesh;
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
        TPZBndCond *mat = dynamic_cast<TPZBndCond *>(meshVec[0]->MaterialVec()[cel->Reference()->MaterialId()]);
        if (mat && mat->Type() == 0) {
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



void
ReadGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec, const std::string &prefix,
          const bool &print, const REAL &scale) {
    TPZGmshReader meshReader;
    meshReader.SetfDimensionlessL(scale);
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
    if(print){
        std::string meshFileName = prefix + "gmesh";
        const size_t strlen = meshFileName.length();
        meshFileName.append(".vtk");
        std::ofstream outVTK(meshFileName.c_str());
        meshFileName.replace(strlen, 4, ".txt");
        std::ofstream outTXT(meshFileName.c_str());

        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
        gmesh->Print(outTXT);
        outTXT.close();
        outVTK.close();        
    }

    return;
}

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh,
                 int pOrder, TPZVec<int> &matIdVec, TPZVec<STATE> &urVec,
                 TPZVec<STATE> &erVec, REAL f0, bool isCutOff, const std::string &prefix,
                 const bool &print, const REAL &scale, const bool &hasPML, const REAL &alphaMax) {
    enum {
      dirichlet = 0
    }; // tipo da condicao de contorno do problema
    const int dim = 2;   // dimensao do problema

    TPZVec<int> volMatId(matIdVec.size()-1);
    for(int i = 0; i< volMatId.size(); i++){
      volMatId[i] = matIdVec[i];
    }
    const int boundMatId = matIdVec[matIdVec.size()-1];
    const int nPML = hasPML ? 8 : 0;
    if(volMatId.size() - nPML !=urVec.size()) DebugStop();

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

    ///PML Materials
    TPZVec<int> pmlIds(0,0);
    int outerMaterial = 0;
    if(hasPML){
        pmlIds.Resize(8);
        for(int i = 0; i < 8 ; i++){
            pmlIds[i] = volMatId[volMatId.size()-8 + i];
        }
        volMatId.Resize(volMatId.size()-8);
        outerMaterial = volMatId.size()-1;

        for (int i = 0; i < cmeshH1->NElements(); ++i) {
            TPZCompEl *compel = cmeshH1->Element(i);
            TPZGeoEl *geo = cmeshH1->Element(i)->Reference();
            const int matId = geo->MaterialId();
            if (matId == pmlIds[0] ||
                matId == pmlIds[1] ||
                matId == pmlIds[2] ||
                matId == pmlIds[3] ||
                matId == pmlIds[4] ||
                matId == pmlIds[5] ||
                matId == pmlIds[6] ||
                matId == pmlIds[7]) {
                TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *>(compel);
                cel->PRefine(1);
            }
        }
        for (int i = 0; i < cmeshHCurl->NElements(); ++i) {
            TPZCompEl *compel = cmeshHCurl->Element(i);
            TPZGeoEl *geo = cmeshHCurl->Element(i)->Reference();
            const int matId = geo->MaterialId();
            if (matId == pmlIds[0] ||
                matId == pmlIds[1] ||
                matId == pmlIds[2] ||
                matId == pmlIds[3] ||
                matId == pmlIds[4] ||
                matId == pmlIds[5] ||
                matId == pmlIds[6] ||
                matId == pmlIds[7]) {
                TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *>(compel);
                cel->PRefine(1);
            }
        }

        cmeshH1->ExpandSolution();
        cmeshH1->ComputeNodElCon();
        cmeshH1->CleanUpUnconnectedNodes();

        cmeshHCurl->ExpandSolution();
        cmeshHCurl->ComputeNodElCon();
        cmeshHCurl->CleanUpUnconnectedNodes();

    }

    if (isCutOff) {
        for (int i = 0; i < volMatId.size(); ++i) {
            matMultiPhysics[i] =
                new TPZMatWaveguideCutOffAnalysis(volMatId[i], f0, urVec[i], erVec[i], 1./scale);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    } else {
        for (int i = 0; i < volMatId.size(); ++i) {
            matMultiPhysics[i] =
                new TPZMatModalAnalysis(volMatId[i], f0, urVec[i], erVec[i], 1./scale);
            cmeshMF->InsertMaterialObject(matMultiPhysics[i]);
        }
    }

    if(hasPML){
        REAL xR=1e20,xL=-1e20,yU=1e20,yD=-1e20,d=-1e20;
        TPZGeoMesh * gmesh = cmeshMF->Reference();
        const int urId = pmlIds[5];
        const int llId = pmlIds[7];
        for (int i = 0; i < gmesh->NElements(); ++i) {
            TPZGeoEl *geo = gmesh->Element(i);
            if(geo->MaterialId() == urId){
                for (int j = 0; j < geo->NCornerNodes(); ++j) {
                    TPZManVector<REAL,3> co(3);
                    geo->Node(j).GetCoordinates(co);
                    const REAL & xP = co[0];
                    const REAL & yP = co[1];
                    if( xP < xR ){
                        xR = xP;
                    }
                    if( yP < yU ){
                        yU = yP;
                    }
                    if( xP > d){
                        d = xP;
                    }
                }
            }
            else if(geo->MaterialId() == llId){
                for (int j = 0; j < geo->NCornerNodes(); ++j) {
                    TPZManVector<REAL,3> co(3);
                    geo->Node(j).GetCoordinates(co);
                    const REAL & xP = co[0];
                    const REAL & yP = co[1];
                    if( xP > xL ){
                        xL = xP;
                    }
                    if( yP > yD ){
                        yD = yP;
                    }
                }
            }
            else{
                continue;
            }
        }
        d = d - xR;
        REAL xP;
        REAL yP;
        bool attx = true;
        bool atty = false;
        for(int i = 0; i < nPML; i++){
            xP = (i / 2) % 2 ? xL : xR;
            yP = !((i / 2) % 2) != !(i%2)? yU : yD;
            if(i <4){
                attx = !(i % 2);
                atty = !(attx);
            }
            if(i == 4){
                attx = true;
                atty = true;
            }
            matMultiPhysics[volMatId.size() + i] =
                    new TPZMatWaveguidePml(pmlIds[i],
                                           *matMultiPhysics[outerMaterial],
                                           attx,xP,
                                           atty,yP,
                                           alphaMax,d);
            cmeshMF->InsertMaterialObject(matMultiPhysics[volMatId.size() + i]);
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
    for (int i = 0; i < nPML; ++i) {
        set.insert(pmlIds[i]);
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

    cmeshMF->ExpandSolution();
    cmeshMF->ComputeNodElCon();
    cmeshMF->CleanUpUnconnectedNodes();

    if(print){
        std::ofstream fileH1(prefix + "cmeshH1.txt");
        cmeshH1->Print(fileH1);
        std::ofstream fileHCurl(prefix + "cmeshHCurl.txt");
        cmeshHCurl->Print(fileHCurl);
        std::ofstream fileMF(prefix + "cmeshMFHCurl.txt");
        cmeshMF->Print(fileMF);
    }

    meshVecOut.resize(3);

    meshVecOut[0] = cmeshMF;
    meshVecOut[1 + matMultiPhysics[0]->H1Index()] = cmeshH1;
    meshVecOut[1 + matMultiPhysics[0]->HCurlIndex()] = cmeshHCurl;
    return;
}
