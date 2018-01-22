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
#include <TPZEigenSolver.h>
#include <TPZSpStructMatrix.h>
#include <tpzgeoelmapped.h>
#include <tpzarc3d.h>
#include "pzextractval.h"
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
#include "TPZMatWaveguideCutOffAnalysis.h"
#include "TPZMatModalAnalysis.h"
#include "TPZGmshReader.h"

void RunSimulation(const bool &isCutOff, const std::string &mshFileName, const int &pOrder,
                   const REAL &f0, const bool &genVTK,
                   const bool &l2error, const bool &exportEigen, const int &nThreads,
                   const bool &optimizeBandwidth, const bool &filterEquations);

void ReadGMesh(TPZGeoMesh *&gmesh, const std::string mshFileName, TPZVec<int> &matIdVec);

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh,
                 int pOrder, TPZVec<int> &matIdVec, TPZVec<STATE> &urVec,
                 TPZVec<STATE> &er, REAL f0, bool isCutOff);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> cmeshMF,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);


#define HARDCODED_WG

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif



    //Command-Line interactive mode
    #ifndef HARDCODED_WG
    std::cout<<"Input .msh file name: ";
    std::string mshFileName;
    std::cin >> mshFileName;
    if(!mshFileName.empty()) DebugStop();

    bool isCutOff = false;
    REAL fOp = 9e+9;
    std::cout<<"Cut-off analyser? Y/N (default: N) : ";
    {
        std::string input;
        std::getline( std::cin, input );
        if ( !input.empty() ) {
            if(input[0] == 'Y') isCutOff = true;
        }
    }

    if(!isCutOff){
        std::cout<<"Input operational frequency in Hz (default: 9GHz) : ";
        std::string input;
        std::getline( std::cin, input );
        if ( !input.empty() ) {
            std::istringstream stream( input );
            stream >> fOp;
        }
    }

    #else
    //hard-coded mode
    std::string mshFileName = "coarseMesh.msh";
    REAL fOp = 5e+9;
    bool isCutOff = false;//analysis of cutoff frequencies for eigenmodes
    #endif

    int pOrder = 1;           // polynomial order of basis functions
    bool genVTK = false;      // generate vtk for fields visualisation
    bool l2error = false;     // TODO: implement error analysis
    bool exportEigen = false; // export eigen values
    const int nThreads = 8;
    bool optimizeBandwidth = true; //whether to renumber equations (OFF for debugging purposes)
    bool filterEquations = true; //whether to impose dirichlet conditions removing boundary equations


    const int nSim = 1;
    for (int i = 0; i < nSim; i++) {
        std::cout << "iteration " << i + 1 << " of " << nSim << std::endl;

        RunSimulation(isCutOff, mshFileName, pOrder, fOp, genVTK, l2error, exportEigen, nThreads, optimizeBandwidth,
                      filterEquations);
    }

    return 0;
}

void RunSimulation(const bool &isCutOff, const std::string &mshFileName, const int &pOrder,
                   const REAL &f0, const bool &genVTK,
                   const bool &l2error, const bool &exportEigen, const int &nThreads,
                   const bool &optimizeBandwidth, const bool &filterEquations) {
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    std::cout<<"Creating GMesh...";
    #ifdef USING_BOOST
    boost::posix_time::ptime t1_g =
            boost::posix_time::microsec_clock::local_time();
    #endif
    TPZVec<int> matIdVec;
    ReadGMesh(gmesh, mshFileName, matIdVec);
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
    CreateCMesh(meshVec, gmesh, pOrder, matIdVec, f0,
                isCutOff); // funcao para criar a malha computacional
    #ifdef USING_BOOST
    boost::posix_time::ptime t2_c =
            boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_c-t1_c<<std::endl;
    #endif
    TPZCompMesh *cmesh = meshVec[0];
    TPZMatModalAnalysis *matPointer =
        dynamic_cast<TPZMatModalAnalysis *>(cmesh->MaterialVec()[1]);
    TPZVec<TPZCompMesh *> temporalMeshVec(2);
    temporalMeshVec[matPointer->H1Index()] = meshVec[1 + matPointer->H1Index()];
    temporalMeshVec[matPointer->HCurlIndex()] =
        meshVec[1 + matPointer->HCurlIndex()];

    TPZEigenAnalysis an(cmesh, optimizeBandwidth);
    
    TPZManVector<long, 1000> activeEquations;
    int neq = 0;
    int neqOriginal = 0;

  TPZAutoPointer<TPZStructMatrix> strmtrx;
  strmtrx = new TPZSpStructMatrix(cmesh);
//    strmtrx = new TPZFStructMatrix(cmesh);
  strmtrx->SetNumThreads(nThreads);
    if (filterEquations) {
      FilterBoundaryEquations(meshVec, activeEquations, neq, neqOriginal);
      strmtrx->EquationFilter().SetActiveEquations(activeEquations);
    }
    an.SetStructuralMatrix(strmtrx);

    const int nSolutions = neq >= 10 ? 10 : neq;
    TPZEigenSolver<STATE> solver;
    solver.SetAsGeneralised(true);
    solver.SetAbsoluteValue(false);
    solver.SetDesiredPartOfSpectrum(EDesiredEigen::MNE);//Most Negative Eigenvalues
    solver.SetHowManyEigenValues(nSolutions);
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

//    TPZMatrix<STATE> *stiffAPtr = NULL, *stiffBPtr = NULL;
//    stiffAPtr = new TPZFMatrix<STATE>(
//            *dynamic_cast<TPZFMatrix<STATE> *>(an.Solver().MatrixA().operator->()));
//    stiffBPtr = new TPZFMatrix<STATE>(
//        *dynamic_cast<TPZFMatrix<STATE> *>(an.Solver().MatrixB().operator->()));

    std::cout << "Solving..." << std::endl;
#ifdef USING_BOOST
    boost::posix_time::ptime t3 =
        boost::posix_time::microsec_clock::local_time();
#endif
//    TPZManVector<STATE,10> eValues = solver.GetEigenvalues();
//    TPZFMatrix<STATE> eVectors = solver.GetEigenvectors();
//    stiffAPtr->SolveGeneralisedEigenProblem( *stiffBPtr, eValues , eVectors);
    an.Solve();
#ifdef USING_BOOST
    boost::posix_time::ptime t4 =
        boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2 - t1 << " Time for solving "
              << t4 - t3 << std::endl;
#endif
    return;
    TPZManVector<SPZAlwaysComplex<STATE>::type> eValues = an.GetEigenvalues();
    TPZFMatrix<SPZAlwaysComplex<STATE>::type> eVectors = an.GetEigenvectors();
    std::set<std::pair<REAL, TPZFMatrix<STATE>>> eigenValuesRe;
    TPZFMatrix<SPZAlwaysComplex<STATE>::type> eVector;
    std::pair<REAL, TPZFMatrix<STATE>> duplet;

    for (int i = 0; i < eValues.size(); i++) {
        eVectors.GetSub(0, i, eVectors.Rows(), 1, eVector);
        duplet.first = eValues[i].real();
        duplet.second = eVector;
        eigenValuesRe.insert(duplet);
    }
    if (exportEigen) {
		int i = 0;
		std::string fileName;
		fileName = "ev";
		fileName.append(std::to_string(nDiv));
		fileName.append("_p");
		fileName.append(std::to_string(pOrder));
		fileName.append(".csv");
			std::ofstream fileEigenValues(fileName.c_str());
			for (std::set<std::pair<REAL, TPZFMatrix<STATE>>>::iterator iT =
					 eigenValuesRe.begin();
				 iT != eigenValuesRe.end(); iT++) {
				if (isCutOff) {
					if (std::abs(iT->first) < 1e-2)
						continue;
				}
				fileEigenValues << iT->first << std::endl;
				i++;
			}
    }
	int i = 0;

    for (std::set<std::pair<REAL, TPZFMatrix<STATE>>>::iterator iT =
             eigenValuesRe.begin();
         iT != eigenValuesRe.end(); iT++) {
        if (isCutOff) {
			if (std::abs(iT->first) < 1e-2) {
				continue;
			}
		}
		if (i >= nSolutions) {
			break;
		}
		std::cout << iT->first << std::endl;
		i++;
    }
//
//    if (isCutOff) {
//        std::cout << "FINISHED!" << std::endl;
//        return;
//    }
//    if (exportEigen) {
//        std::string fileName("evectors");
//        fileName.append(std::to_string(nDiv));
//        fileName.append("_p");
//        fileName.append(std::to_string(pOrder));
//        fileName.append(".csv");
//        std::ofstream fileA(fileName.c_str());
//        char number[256];
//        for (int i = 0; i < eigenValuesRe.begin()->second.Rows(); i++) {
//            for (int j = 0; j < eigenValuesRe.begin()->second.Cols(); j++) {
//                sprintf(number, "%32.32Lf",
//                        (long double)TPZExtractVal::val(std::real(
//                            eigenValuesRe.begin()->second.GetVal(i, j))));
//                fileA << number;
//
//                fileA << " + I * ";
//
//                sprintf(number, "%32.32Lf",
//                        (long double)TPZExtractVal::val(std::imag(
//                            eigenValuesRe.begin()->second.GetVal(i, j))));
//                fileA << number;
//                if (j != eigenValuesRe.begin()->second.Cols() - 1) {
//                    fileA << " , ";
//                }
//            }
//            fileA << std::endl;
//        }
//        fileA.close();
//    }
//    if (genVTK) {
//        std::cout << "Post Processing..." << std::endl;
//        TPZFMatrix<STATE> solMat(neqOriginal, 1, 0.);
//
//        TPZStack<std::string> scalnames, vecnames;
//        scalnames.Push("Ez"); // setando para imprimir u
//        vecnames.Push("Et");
//        std::string plotfile = "waveguideModes.vtk"; // arquivo de saida que
//                                                        // estara na pasta debug
//        int dim = 2;
//        an.DefineGraphMesh(dim, scalnames, vecnames,
//                           plotfile);  // define malha grafica
//        int postProcessResolution = 2; // define resolucao do pos processamento
//
//        std::set<std::pair<REAL, TPZFMatrix<STATE>>>::iterator iT =
//            eigenValuesRe.begin();
//        for (int iSol = 0; iSol < nSolutions; iSol++) {
//            if (iT == eigenValuesRe.end()) {
//                DebugStop();
//            }
//            for (int i = 0; i < neq; i++) {
//                solMat(activeEquations[i], 0) = (iT->second).GetVal(i, 0);
//            }
//            an.LoadSolution(solMat);
//            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(temporalMeshVec,
//                                                               cmesh);
//            an.PostProcess(postProcessResolution);
//            iT++;
//        }
//    }
//    if (l2error) {
//        DebugStop(); // yet to be implemented
//    }
    std::cout << "FINISHED!" << std::endl;
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

    TPZStack<int> matIds = meshReader.fMaterialDataVec.fMatID;
    matIdVec.Resize(0);
    while (matIds.size() != 0){
        matIdVec.Resize(matIdVec.size()+1);
        matIdVec[matIdVec.size()-1] = matIds.Pop();
        std::cout<<"Read material "<<matIdVec[matIdVec.size()-1]<<std::endl;
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

    const int dim = 2;   // dimensao do problema

    TPZVec<int> boundMatId, volMatId;
    for(int i = 0; i< matIdVec.size(); i++){
        if(matIdVec[i]>0){
            volMatId.Resize(volMatId.size()+1);
            volMatId[volMatId.size()-1] = matIdVec[i];
        }
        else{
            boundMatId.Resize(boundMatId.size()+1);
            boundMatId[boundMatId.size()-1] = matIdVec[i];
        }
    }

    urVec.Resize(volMatId.size());
    erVec.Resize(volMatId.size());

    for (int i = 0; i < volMatId.size(); ++i) {
        erVec[i] = 1.;
        std::cout<<"Input er for material "<<volMatId[i]<<" (default: 1): ";
        {
            std::string input;
            std::getline( std::cin, input );
            if ( !input.empty() ) {
                std::istringstream stream( input );
                stream >> erVec[i];
            }
        }

        urVec[i] = 1.;
        std::cout<<"Input ur for material "<<volMatId[i]<<" (default: 1): ";
        {
            std::string input;
            std::getline( std::cin, input );
            if ( !input.empty() ) {
                std::istringstream stream( input );
                stream >> urVec[i];
            }
        }
    }
    enum {
        dirichlet = 0,
        neumann,
        mixed
    }; // tipo da condicao de contorno do problema
    // Criando material

    /// criar malha computacional H1

    TPZCompMesh *cmeshH1 = new TPZCompMesh(gmesh);
    cmeshH1->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmeshH1->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha
    const int nState = 1;
    TPZVec<STATE> sol; // only for creating material. this material will not be
                       // used in reality
    TPZL2Projection *matH1 = new TPZL2Projection(volMatId[0], dim, nState, sol);
    cmeshH1->InsertMaterialObject(matH1);

    /// electrical conductor boundary conditions
    TPZFNMatrix<1, STATE> val1(1, 1, 0.), val2(1, 1, 0.);

    TPZMaterial *BCondH1Dir = matH1->CreateBC(
        matH1, boundMatId[0], dirichlet, val1, val2); // cria material que implementa a
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
    TPZMatHCurlProjection *matHCurl = new TPZMatHCurlProjection(volMatId[0]);
    cmeshHCurl->InsertMaterialObject(matHCurl);

    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *BCondHCurlDir =
        matHCurl->CreateBC(matHCurl, bc0, dirichlet, val1,
                           val2); // cria material que implementa a condicao de
                                  // contorno de dirichlet

    cmeshHCurl->InsertMaterialObject(BCondHCurlDir); // insere material na malha
    
    cmeshHCurl->SetAllCreateFunctionsHCurl(); // define espaco de aproximacao
    
    cmeshHCurl->AutoBuild();
    cmeshHCurl->CleanUpUnconnectedNodes();

    TPZMatModalAnalysis *matMultiPhysics = NULL;
    TPZVec<TPZCompMesh *> meshVec(2);
    if (isCutOff) {
        TPZMatWaveguideCutOffAnalysis *dummy =
            new TPZMatWaveguideCutOffAnalysis(matId, f0, ur, er);
        matMultiPhysics = dummy;
    } else {
        TPZMatModalAnalysis *dummy = NULL;
        dummy = new TPZMatModalAnalysis(matId, f0, ur, er);
        matMultiPhysics = dummy;
    }
    meshVec[matMultiPhysics->H1Index()] = cmeshH1;
    meshVec[matMultiPhysics->HCurlIndex()] = cmeshHCurl;

    TPZCompMesh *cmeshMF = new TPZCompMesh(gmesh);

    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *BCondMFDir =
        matMultiPhysics->CreateBC(matMultiPhysics, bc0, dirichlet, val1,
                                  val2); // cria material que implementa a
                                         // condicao de contorno de dirichlet

    cmeshMF->InsertMaterialObject(matMultiPhysics);
    cmeshMF->InsertMaterialObject(BCondMFDir); // insere material na malha
    std::set<int> set;
    set.insert(matId);
    set.insert(bc0);

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild(set);
    cmeshMF->CleanUpUnconnectedNodes();

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
    meshVecOut[1 + matMultiPhysics->H1Index()] = cmeshH1;
    meshVecOut[1 + matMultiPhysics->HCurlIndex()] = cmeshHCurl;
    return;
}
