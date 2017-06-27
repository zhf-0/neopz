/**
 * @file
 * @brief Afirst attempt at a multi physics hcurl/h1 formulation
 * @details Adequadte for problems with longitudinal axis symmetry
 * such as some setions of waveguides (closed waveguides).
 * it uses hcurl for tranverse componentes and h1 for longitudinal components
 *
 * @author Francisco Orlandini
 * @since 2015
 */

#include <iostream>
#include <fstream>
#include <string>
#include "pzextractval.h"
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "TPZTimer.h"
#include "pzanalysis.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "pzl2projection.h"
#include "pzfstrmatrix.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZMatMFHCurlH1.h"
#include "TPZMatWaveguideCutOffAnalysis.h"

enum meshTypeE { createRectangular = 1, createTriangular, createZigZag };

STATE ur(const TPZVec<REAL> &x) { return 1.; }
STATE er(const TPZVec<REAL> &x) { return 1.; }

void RunSimulation(bool isCutOff, const int meshType, int pOrder, int nDiv,
                   REAL hDomain, REAL wDomain, REAL f0, bool genVTK,
                   bool l2error, const int nThreads, bool optimizeBandwidth,
                   bool filterEquations);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> cmeshMF,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh, int pOrder,
                            STATE (&ur)(const TPZVec<REAL> &),
                            STATE (&er)(const TPZVec<REAL> &), REAL f0,
                            bool isCutOff);

void CreateGMesh(TPZGeoMesh *&gmesh, const int meshType, const REAL hDomain,
                 const REAL wDomain, const int xDiv, const int zDiv);

int main(int argc, char *argv[]) {

    HDivPiola = 1; // use piola mapping
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // PARAMETROS FISICOS DO PROBLEMA
    const REAL hDomain = 4 * 2.54 * 1e-3;
    const REAL wDomain = 9 * 2.54 * 1e-3;
    const REAL f0 = 25 * 1e+9;

    bool isCutOff = false;
    const enum meshTypeE meshType = createTriangular;
    int pOrder = 2; // polynomial order of basis functions
    bool genVTK = false;
    bool l2error = false;// TODO: implement error analysis
    const int nThreads = 0; // TODO: fix multithread issue
    bool optimizeBandwidth = true;
    bool filterEquations = true;

    int nDiv = 4;
    const int nSim = 2;
    for (int i = 0; i < nSim; i++) {
        std::cout << "iteration " << i + 1 << " of " << nSim << std::endl;
        RunSimulation(isCutOff, meshType, pOrder, nDiv, hDomain, wDomain, f0,
                      genVTK, l2error, nThreads, optimizeBandwidth,
                      filterEquations);
        nDiv *= 2;
    }

    return 0;
}

void RunSimulation(bool isCutOff, const int meshType, int pOrder, int nDiv,
				   REAL hDomain, REAL wDomain, REAL f0, bool genVTK,
				   bool l2error, const int nThreads, bool optimizeBandwidth,
				   bool filterEquations) {
    TPZTimer timer;

    timer.start();

    int xDiv = nDiv;
    int zDiv = nDiv;
    TPZGeoMesh *gmesh = new TPZGeoMesh();

    CreateGMesh(gmesh, meshType, hDomain, wDomain, xDiv, zDiv);

	TPZVec<TPZCompMesh *> meshVec(1);
    CreateCMesh(meshVec,gmesh, pOrder, ur, er, f0,
              isCutOff); // funcao para criar a malha computacional
    TPZCompMesh *cmeshMF = meshVec[0];
    TPZMatMFHCurlH1 *matPointer =
        dynamic_cast<TPZMatMFHCurlH1 *>(cmeshMF->MaterialVec()[1]);
    TPZVec<TPZCompMesh *> temporalMeshVec(2);
    temporalMeshVec[matPointer->H1Index()] = meshVec[1 + matPointer->H1Index()];
    temporalMeshVec[matPointer->HCurlIndex()] =
        meshVec[1 + matPointer->HCurlIndex()];

    TPZAnalysis an(cmeshMF, optimizeBandwidth);
    //configuring analysis object
	
    TPZManVector<long, 1000> activeEquations;
    int neq = 0;
    int neqOriginal = 0;

    TPZAutoPointer<TPZFStructMatrix> fmtrx;

    fmtrx = new TPZFStructMatrix(cmeshMF);
    fmtrx->SetNumThreads(nThreads);
	if(filterEquations){
		FilterBoundaryEquations(meshVec, activeEquations, neq, neqOriginal);
		fmtrx->EquationFilter().SetActiveEquations(activeEquations);
	}
    an.SetStructuralMatrix(fmtrx);

    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky); // caso simetrico
    an.SetSolver(step);
    int matId = 1;
    TPZMatMFHCurlH1 *matAlias =
        dynamic_cast<TPZMatMFHCurlH1 *>(cmeshMF->FindMaterial(matId));

    matAlias->SetMatrixA();
	std::cout << "entrando no assemble matrix A" << std::endl;
    an.Assemble();
	
    TPZMatrix<STATE> *stiffAPtr = NULL, *stiffBPtr = NULL;
    stiffAPtr = new TPZFMatrix<STATE>(
        *dynamic_cast<TPZFMatrix<STATE> *>(an.Solver().Matrix().operator->()));
	
    matAlias->SetMatrixB();
    std::cout << "entrando no assemble matrix B" << std::endl;
    an.Assemble();
	
    stiffBPtr = new TPZFMatrix<STATE>(
        *dynamic_cast<TPZFMatrix<STATE> *>(an.Solver().Matrix().operator->()));
    std::cout << "saindo do assemble" << std::endl;

    TPZVec<STATE> eValues;
    TPZFMatrix<STATE> eVectors;
    std::cout << "entrando no calculo dos autovalores" << std::endl;
    stiffAPtr->SolveGeneralisedEigenProblem(*stiffBPtr, eValues, eVectors);
    std::cout << "saindo do calculo dos autovalores" << std::endl;
    timer.stop();

    std::set<std::pair<REAL, TPZFMatrix<STATE>>> eigenValuesRe;
    TPZFMatrix<STATE> eVector(eVectors.Rows(), 1);
    std::pair<REAL, TPZFMatrix<STATE>> duplet;

    for (int i = 0; i < eValues.size(); i++) {
        eVectors.GetSub(0, i, eVectors.Rows(), 1, eVector);
        duplet.first = eValues[i].real();
        duplet.second = eVector;
        eigenValuesRe.insert(duplet);
    }
    int i = 0;

    std::string fileName;
    fileName = "../ev";
    fileName.append(std::to_string(nDiv));
	fileName.append("_p");
	fileName.append(std::to_string(pOrder));
    fileName.append(".csv");

	const int nSolutions = neq >= 10 ? 10 : neq;
    std::ofstream fileEigenValues(fileName.c_str());
    for (std::set<std::pair<REAL, TPZFMatrix<STATE>>>::iterator iT =
             eigenValuesRe.begin();
         iT != eigenValuesRe.end(); iT++) {
        if (isCutOff) {
            if (std::abs(iT->first) < 1e-2)
                continue;
            std::cout << iT->first << std::endl;
            fileEigenValues << iT->first << std::endl;
			if (i >= nSolutions){
                break;
			}
			i++;
        } else {
            fileEigenValues << iT->first << std::endl;
			if (i >= nSolutions){
                continue;
			}
			if (i > 50){
                break;
			}
            std::cout << iT->first << std::endl;
			i++;
        }
    }
    if (isCutOff) {
        std::cout << "FINISHED!" << std::endl;
        return;
    }

    {
        std::string fileName("../evectors");
        fileName.append(std::to_string(nDiv));
		fileName.append("_p");
		fileName.append(std::to_string(pOrder));
        fileName.append(".csv");
        std::ofstream fileA(fileName.c_str());
        char number[256];
        for (int i = 0; i < eigenValuesRe.begin()->second.Rows(); i++) {
            for (int j = 0; j < eigenValuesRe.begin()->second.Cols(); j++) {
                sprintf(number, "%32.32Lf",
                        (long double)TPZExtractVal::val(std::real(
                            eigenValuesRe.begin()->second.GetVal(i, j))));
                fileA << number;

                fileA << " + I * ";

                sprintf(number, "%32.32Lf",
                        (long double)TPZExtractVal::val(std::imag(
                            eigenValuesRe.begin()->second.GetVal(i, j))));
                fileA << number;
                if (j != eigenValuesRe.begin()->second.Cols() - 1) {
                    fileA << " , ";
                }
            }
            fileA << std::endl;
        }
        fileA.close();
    }
    if (genVTK) {
		std::cout << "Post Processing..." << std::endl;
        TPZFMatrix<STATE> solMat(neqOriginal, 1, 0.);

        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Ez"); // setando para imprimir u
        vecnames.Push("Et");
        std::string plotfile = "../waveguideModes.vtk"; // arquivo de saida que
                                                        // estara na pasta debug
        int dim = 2;
        an.DefineGraphMesh(dim, scalnames, vecnames,
                           plotfile);  // define malha grafica
        int postProcessResolution = 2; // define resolucao do pos processamento

        std::set<std::pair<REAL, TPZFMatrix<STATE>>>::iterator iT =
            eigenValuesRe.begin();
        for (int iSol = 0; iSol < nSolutions; iSol++) {
            if (iT == eigenValuesRe.end()) {
                DebugStop();
            }
            for (int i = 0; i < neq; i++) {
                solMat(activeEquations[i], 0) = (iT->second).GetVal(i, 0);
            }
            an.LoadSolution(solMat);
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(temporalMeshVec,
                                                               cmeshMF);
            an.PostProcess(postProcessResolution);
            iT++;
        }
    }
	if (l2error) {
		DebugStop();//yet to be implemented
	}
    std::cout << "FINISHED!" << std::endl;
}

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> meshVec,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal) {
    TPZCompMesh *cmeshMF = meshVec[0];
    TPZCompMesh *cmeshHCurl = meshVec[1];
    TPZCompMesh *cmeshH1 = meshVec[2];

    TPZManVector<long, 1000> allConnects;
    std::set<long> boundConnects;

    for (int iel = 0; iel < cmeshMF->NElements(); iel++) {
        TPZCompEl *cel = cmeshMF->ElementVec()[iel];
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
    for (int iCon = 0; iCon < cmeshMF->NConnects(); iCon++) {
        if (boundConnects.find(iCon) == boundConnects.end()) {
            TPZConnect &con = cmeshMF->ConnectVec()[iCon];
            if (con.HasDependency())
                continue;
            int seqnum = con.SequenceNumber();
            int pos = cmeshMF->Block().Position(seqnum);
            int blocksize = cmeshMF->Block().Size(seqnum);
            if (blocksize == 0)
                continue;

            int vs = activeEquations.size();
            activeEquations.Resize(vs + blocksize);
            for (int ieq = 0; ieq < blocksize; ieq++) {
                activeEquations[vs + ieq] = pos + ieq;
            }
        }
    }

    neqOriginal = cmeshMF->NEquations();

    int nHCurlEquations = 0, nH1Equations = 0;
    long nEq = 0;
    TPZMatMFHCurlH1 *mat =
        dynamic_cast<TPZMatMFHCurlH1 *>(cmeshMF->FindMaterial(1));
    for (int iCon = 0; iCon < cmeshMF->NConnects(); iCon++) {
        bool isH1;
        if (boundConnects.find(iCon) == boundConnects.end()) {
            if (cmeshMF->ConnectVec()[iCon].HasDependency())
                continue;
            int seqnum = cmeshMF->ConnectVec()[iCon].SequenceNumber();
            int blocksize = cmeshMF->Block().Size(seqnum);
            if (mat->H1Index() == 0 && iCon < cmeshH1->NConnects()) {
                isH1 = true;
            } else if (mat->H1Index() == 1 && iCon >= cmeshHCurl->NConnects()) {
                isH1 = true;
            } else {
                isH1 = false;
            }
            for (int ieq = 0; ieq < blocksize; ieq++) {
                nEq++;
                isH1 == true ? nH1Equations++ : nHCurlEquations++;
            }
        }
    }
    std::cout << "------\tactive eqs\t-------" << std::endl;
    std::cout << "# H1 equations: " << nH1Equations << std::endl;
    std::cout << "# HCurl equations: " << nHCurlEquations << std::endl;
    std::cout << "# equations: " << nEq << std::endl;
    std::cout << "------\t----------\t-------" << std::endl;
    neq = nEq;
    return;
}

void CreateGMesh(TPZGeoMesh *&gmesh, const int meshType, const REAL hDomain,
                 const REAL wDomain, const int xDiv, const int yDiv) {

    TPZManVector<int, 3> nx(3, 0);
    TPZManVector<REAL, 3> llCoord(3, 0.), ulCoord(3, 0.), urCoord(3, 0.),
        lrCoord(3, 0.);
    llCoord[0] = 0.;
	llCoord[1] = 0.;

	ulCoord[0] = 0.;
    ulCoord[1] = hDomain;

    urCoord[0] = wDomain;
    urCoord[1] = hDomain;

    lrCoord[0] = wDomain;
	lrCoord[1] = 0.;

    nx[0] = xDiv;
    nx[1] = yDiv;
    int numl = 1;
    TPZGenGrid *gengrid = NULL;
    switch (meshType) {
    case createRectangular: {
        DebugStop(); // No quadrilateral Hcurl elements yet!
        REAL rot = 0.0;
        gengrid = new TPZGenGrid(nx, llCoord, urCoord, numl, rot);
        gengrid->SetElementType(EQuadrilateral);
    } break;
    case createTriangular: {
        REAL rot = 0.0;
        gengrid = new TPZGenGrid(nx, llCoord, urCoord, numl, rot);
        gengrid->SetElementType(ETriangle);
    } break;
    case createZigZag: {
        DebugStop(); // No quadrilateral Hcurl elements yet!
        REAL rot = 0.0;
        gengrid = new TPZGenGrid(nx, llCoord, urCoord, numl, rot);
        gengrid->SetElementType(EQuadrilateral);
        gengrid->SetZigZagPattern();
    } break;
    default:
        DebugStop();
        break;
    }
    gmesh = new TPZGeoMesh();
    const int matId = 1; // define id para um material(formulacao fraca)
    const int bc0 = -1;  // define id para um material(cond contorno dirichlet)
    gengrid->Read(gmesh, matId);

    gengrid->SetBC(gmesh, ulCoord, llCoord, bc0);
    gengrid->SetBC(gmesh, urCoord, ulCoord, bc0);
    gengrid->SetBC(gmesh, lrCoord, urCoord, bc0);
    gengrid->SetBC(gmesh, llCoord, lrCoord, bc0);

    gmesh->BuildConnectivity();
#ifdef PZDEBUG
    std::ofstream outTxt, outVtk;
    switch (meshType) {
    case createRectangular: {
        outTxt.open("../gmeshRectangular.txt"); // define arquivo de saida para
                                                // impressao da malha no
        outVtk.open("../gmeshRectangular.vtk"); // define arquivo de saida para
                                                // impressao da malha no
                                                // paraview

    } break;
    case createTriangular:

    {
        outTxt.open("../gmeshTriangular.txt"); // define arquivo de saida para
                                               // impressao da malha no
        outVtk.open("../gmeshTriangular.vtk"); // define arquivo de saida para
                                               // impressao da malha no paraview

    } break;
    case createZigZag: {
        outTxt.open("../gmeshZigZag.txt"); // define arquivo de saida para
                                           // impressao da malha no
        outVtk.open("../gmeshZigZag.vtk"); // define arquivo de saida para
                                           // impressao da malha no paraview
    } break;
    default:
        DebugStop();
        break;
    }
    gmesh->Print(outTxt);
    outTxt.close();
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVtk,
                                 true); // imprime a malha no formato vtk
    outVtk.close();
#endif
    delete gengrid;
}

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh, int pOrder,
                            STATE (&ur)(const TPZVec<REAL> &),
                            STATE (&er)(const TPZVec<REAL> &), REAL f0,
                            bool isCutOff) {

    const int dim = 2;   // dimensao do problema
    const int matId = 1; // define id para um material(formulacao fraca)
    const int bc0 = -1;  // define id para um material(cond contorno dirichlet)
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
    TPZL2Projection *matH1 = new TPZL2Projection(matId, dim, nState, sol);
    cmeshH1->InsertMaterialObject(matH1);

    /// electrical conductor boundary conditions
    TPZFNMatrix<1, STATE> val1(1, 1, 0.), val2(1, 1, 0.);

    TPZMaterial *BCondH1Dir = matH1->CreateBC(
        matH1, bc0, dirichlet, val1, val2); // cria material que implementa a
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
    TPZMatHCurlProjection *matHCurl = new TPZMatHCurlProjection(matId);
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

    TPZMatMFHCurlH1 *matMultiPhysics = NULL;
    TPZVec<TPZCompMesh *> meshVec(2);
    if (isCutOff) {
        TPZMatWaveguideCutOffAnalysis *dummy =
            new TPZMatWaveguideCutOffAnalysis(matId, f0, ur, er);
        matMultiPhysics = dummy;
    } else {
        TPZMatMFHCurlH1 *dummy = new TPZMatMFHCurlH1(matId, f0, ur, er);
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
    std::ofstream fileH1("../cmeshH1.txt");
    cmeshH1->Print(fileH1);
    std::ofstream fileHCurl("../cmeshHCurl.txt");
    cmeshHCurl->Print(fileHCurl);
    std::ofstream fileMF("../cmeshMFHCurl.txt");
    cmeshMF->Print(fileMF);

    meshVecOut.resize(3);

    meshVecOut[0] = cmeshMF;
    meshVecOut[1 + matMultiPhysics->H1Index()] = cmeshH1;
    meshVecOut[1 + matMultiPhysics->HCurlIndex()] = cmeshHCurl;

    std::cout << "------\t----------\t-------" << std::endl;
    std::cout << "cmeshH1->NEquations()"
              << "\t" << cmeshH1->NEquations() << std::endl;
    std::cout << "cmeshHCurl->NEquations()"
              << "\t" << cmeshHCurl->NEquations() << std::endl;
    std::cout << "cmeshMF->NEquations()"
              << "\t" << cmeshMF->NEquations() << std::endl;
    std::cout << "------\t----------\t-------" << std::endl;
    return;
}
