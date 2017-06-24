/**
 * @file
 * @brief TODO:DESCRIBE IT
 * @details BLABLABLA
 *
 * @author Francisco Orlandini
 * @since 2017
 */

#include "TPZMatHelmholtz2D.h"
#include "TPZTimer.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzelchdiv.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzlog.h"
#include "pzmatred.h"
#include "pzsbndmat.h"
#include "pzsbstrmatrix.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzsubcmesh.h"
#include "pzvisualmatrix.h"
#include "tpzmatredstructmatrix.h"
#include <fstream>
#include <iostream>
#include <string>

enum meshTypeE { createRectangular = 1, createTriangular, createZigZag };

STATE constitutiveFunc(const TPZVec<REAL> &coord) { return 1.; }

void loadVec(const TPZVec<REAL> &coord, TPZVec<STATE> &val) {
    val.Resize(3, 0.);
    //  val[0] = 3 - coord[1] * coord[1];
    //  val[1] = 3 - coord[0] * coord[0];
    val[0] = (2 * M_PI * M_PI + constitutiveFunc(coord)) * M_PI *
             cos(M_PI * coord[0]) * sin(M_PI * coord[1]);
    val[1] = (2 * M_PI * M_PI + constitutiveFunc(coord)) * M_PI * (-1.) *
             sin(M_PI * coord[0]) * cos(M_PI * coord[1]);
}

void exactSol(const TPZVec<REAL> &coord, TPZVec<STATE> &result,
              TPZFMatrix<STATE> &curl) {
    result.Resize(3, 0.);
    result[0] = M_PI * cos(M_PI * coord[0]) * sin(M_PI * coord[1]);
    result[1] = M_PI * (-1.) * sin(M_PI * coord[0]) * cos(M_PI * coord[1]);

    curl.Resize(1, 1);
    curl(0, 0) = -2 * M_PI * M_PI * cos(M_PI * coord[0]) * cos(M_PI * coord[1]);
}

void FilterBoundaryEquations(TPZCompMesh *cmeshHCurl,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);

void CreateCMesh(TPZCompMesh *&cmesh, TPZGeoMesh *gmesh, int pOrder,
                 void (&func)(const TPZVec<REAL> &, TPZVec<STATE> &),
                 STATE (&constFunc)(const TPZVec<REAL> &));

void CreateGMesh(TPZGeoMesh *&gmesh, const int meshType, const REAL hDomain,
                 const REAL wDomain, const int xDiv, const int zDiv);

void RunSimulation(const int nDiv, const int pOrder,
                   const enum meshTypeE meshType, bool filterEquations,
                   bool usingFullMtrx, bool optimizeBandwidth,
                   const int nThreads, bool genVTK, bool l2error,
                   TPZVec<REAL> &errorVec);
int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    int pOrder = 1; // ordem polinomial de aproximacao
    int nDiv = 4;

    bool filterEquations = true;
    bool usingFullMtrx = true;
    bool optimizeBandwidth = true;
    const int nThreads = 0;//TODO: fix multithread issue
    bool genVTK = false;
    bool l2error = true;
    const enum meshTypeE meshType = createTriangular;

    TPZVec<REAL> errorVec(1, 0);
    for (int iDiv = 0; iDiv < 5; iDiv++) {
        std::cout << "beginning simulation with nEl = " << nDiv * nDiv * 2
                  << std::endl;
        RunSimulation(nDiv, pOrder, meshType, filterEquations, usingFullMtrx,
                      optimizeBandwidth, nThreads, genVTK, l2error, errorVec);
        nDiv *= 2;
    }

    return 0;
}

void RunSimulation(const int nDiv, const int pOrder,
                   const enum meshTypeE meshType, bool filterEquations,
                   bool usingFullMtrx, bool optimizeBandwidth,
                   const int nThreads, bool genVTK, bool l2error,
                   TPZVec<REAL> &errorVec) {
    // PARAMETROS FISICOS DO PROBLEMA
    const REAL hDomain = 2;
    const REAL wDomain = 2;

    TPZTimer timer;
    timer.start();

    TPZGeoMesh *gmesh = new TPZGeoMesh();
    CreateGMesh(gmesh, meshType, hDomain, wDomain, nDiv, nDiv);

    TPZCompMesh *cmeshHCurl = NULL;
    CreateCMesh(cmeshHCurl, gmesh, pOrder, loadVec,
                constitutiveFunc); // funcao para criar a malha computacional
	{
		std::string fileName("../cmeshHCurl");
		fileName.append(std::to_string(nDiv*nDiv*2));
		fileName.append(".txt");
		std::ofstream fileHCurl(fileName.c_str());
		cmeshHCurl->Print(fileHCurl);
	}
    TPZAnalysis an(cmeshHCurl, optimizeBandwidth);
    // configuracoes do objeto de analise
    TPZManVector<long, 1000> activeEquations;
    int neq = 0;
    int neqOriginal = 0;
    FilterBoundaryEquations(cmeshHCurl, activeEquations, neq, neqOriginal);

    TPZAutoPointer<TPZSBandStructMatrix> sbstr;

    TPZAutoPointer<TPZFStructMatrix> fmtrx;

    if (usingFullMtrx) {
        fmtrx = new TPZFStructMatrix(cmeshHCurl);
        fmtrx->SetNumThreads(nThreads);
        if (filterEquations) {
            fmtrx->EquationFilter().SetActiveEquations(activeEquations);
        }

        an.SetStructuralMatrix(fmtrx);
    } else {
        sbstr = new TPZSBandStructMatrix(cmeshHCurl);
        sbstr->SetNumThreads(nThreads);
        if (filterEquations) {
            sbstr->EquationFilter().SetActiveEquations(activeEquations);
        }
        an.SetStructuralMatrix(sbstr);
    }
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky); // caso simetrico
    an.SetSolver(step);

    an.Run();
    timer.stop();
    an.LoadSolution();

    if (genVTK) {
        TPZStack<std::string> scalnames, vecnames;
        vecnames.Push("E");
        std::string plotfile = "../sol";
        plotfile.append(std::to_string(cmeshHCurl->NElements()));
        plotfile.append(".vtk");

        an.DefineGraphMesh(2, scalnames, vecnames,
                           plotfile);  // define malha grafica
        int postProcessResolution = 5; // define resolucao do pos processamento
        an.PostProcess(postProcessResolution);
    }
    if (l2error) {
        an.SetExact(&exactSol);
        errorVec.Resize(2, 0.);
        an.PostProcessError(errorVec);
		std::string fileName = "../error";
		fileName.append(std::to_string(nDiv*nDiv*2));
		fileName.append(".csv");
		std::ofstream errorFile(fileName.c_str());
                errorFile << errorVec[0] << "," << errorVec[1] << ","
                          << errorVec[2] << std::endl;
		errorFile.close();
    }
}

void FilterBoundaryEquations(TPZCompMesh *cmeshHCurl,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal) {
    TPZManVector<long, 1000> allConnects;
    std::set<long> boundConnects;

    for (int iel = 0; iel < cmeshHCurl->NElements(); iel++) {
        TPZCompEl *cel = cmeshHCurl->ElementVec()[iel];
        if (cel == NULL) {
            continue;
        }
        if (cel->Reference() == NULL) {
            continue;
        }
        if (cel->Reference()->MaterialId() == -1) {
            std::set<long> boundConnectsEl;
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

    for (int iCon = 0; iCon < cmeshHCurl->NConnects(); iCon++) {
        if (boundConnects.find(iCon) == boundConnects.end()) {
            TPZConnect &con = cmeshHCurl->ConnectVec()[iCon];
            int seqnum = con.SequenceNumber();
            int pos = cmeshHCurl->Block().Position(seqnum);
            int blocksize = cmeshHCurl->Block().Size(seqnum);
            if (blocksize == 0)
                continue;

            int vs = activeEquations.size();
            activeEquations.Resize(vs + blocksize);
            for (int ieq = 0; ieq < blocksize; ieq++) {
                activeEquations[vs + ieq] = pos + ieq;
            }
        }
    }
    neqOriginal = cmeshHCurl->NEquations();

    //  std::cout << "activeEquations" << std::endl;
    long nEq = 0;
    for (int iCon = 0; iCon < cmeshHCurl->NConnects(); iCon++) {
        if (boundConnects.find(iCon) == boundConnects.end()) {
            int seqnum = cmeshHCurl->ConnectVec()[iCon].SequenceNumber();
            int blocksize = cmeshHCurl->Block().Size(seqnum);
            if (blocksize == 0)
                continue;
            nEq++;
        }
    }
    std::cout << "# equations(before): " << neqOriginal << std::endl;
    std::cout << "# equations(after): " << nEq << std::endl;
    neq = nEq;
    return;
}

void CreateGMesh(TPZGeoMesh *&gmesh, const int meshType, const REAL hDomain,
                 const REAL wDomain, const int xDiv, const int yDiv) {
    TPZManVector<int, 3> nx(3, 0);
    TPZManVector<REAL, 3> llCoord(3, 0.), ulCoord(3, 0.), urCoord(3, 0.),
        lrCoord(3, 0.);
    llCoord[0] = -wDomain / 2;
    llCoord[1] = -hDomain / 2;

    ulCoord[0] = -wDomain / 2;
    ulCoord[1] = hDomain / 2;

    urCoord[0] = wDomain / 2;
    urCoord[1] = hDomain / 2;

    lrCoord[0] = wDomain / 2;
    lrCoord[1] = -hDomain / 2;

    nx[0] = xDiv;
    nx[1] = yDiv;
    int numl = 1;
    TPZGenGrid *gengrid = NULL;
    switch (meshType) {
    case createRectangular: {
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

    // gmesh->ResetConnectivities();

    gmesh->BuildConnectivity();
#ifdef PZDEBUG
    std::ofstream outTxt, outVtk;
    switch (meshType) {
    case createRectangular: {
        outTxt.open("../gmeshRectangular.txt"); // define arquivo de saida para
                                                // impressao da malha no
        outVtk.open("../gmeshRectangular.vtk"); // define arquivo de saida para
        // impressao da malha no paraview

    } break;
    case createTriangular: {
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

void CreateCMesh(TPZCompMesh *&cmeshHCurl, TPZGeoMesh *gmesh, int pOrder,
                 void (&loadVec)(const TPZVec<REAL> &, TPZVec<STATE> &),
                 STATE (&constFunc)(const TPZVec<REAL> &)) {
    const int dim = 2;   // dimensao do problema
    const int matId = 1; // define id para um material(formulacao fraca)
    const int bc0 = -1;  // define id para um material(cond contorno dirichlet)
    enum {
        dirichlet = 0,
        neumann,
        mixed
    }; // tipo da condicao de contorno do problema
    // Criando material

    cmeshHCurl = new TPZCompMesh(gmesh);
    cmeshHCurl->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    cmeshHCurl->SetDimModel(dim);        // seta dimensao do modelo
    // Inserindo material na malha
    TPZMatHelmholtz2D *matHCurl = new TPZMatHelmholtz2D(matId, constFunc);
    matHCurl->SetForcingFunction(loadVec, 4);
    cmeshHCurl->InsertMaterialObject(matHCurl);

    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *BCondHCurlDir = matHCurl->CreateBC(
        matHCurl, bc0, dirichlet, val1, val2); // cria material que implementa a
                                               // condicao de contorno de
                                               // dirichlet

    cmeshHCurl->InsertMaterialObject(BCondHCurlDir); // insere material na malha

    cmeshHCurl->SetAllCreateFunctionsHCurl(); // define espaco de aproximacao
    cmeshHCurl->AutoBuild();
}
