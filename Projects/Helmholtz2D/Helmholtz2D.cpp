/**
 * @file
 * @brief TODO:DESCRIBE IT
 * @details BLABLABLA
 *
 * @author Francisco Orlandini
 * @since 2017
 */

#include <TPZSSpStructMatrix.h>
#include <TPZSpStructMatrix.h>
#include <pzskylstrmatrix.h>
#include <pzsbstrmatrix.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "TPZMatHelmholtz2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzelchdiv.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzstepsolver.h"

enum meshTypeE { createRectangular = 1, createTriangular, createZigZag };

void loadVec(const TPZVec<REAL> &coord, TPZVec<STATE> &val) {
    val.Resize(3, 0.);
    //  val[0] = 3 - coord[1] * coord[1];
    //  val[1] = 3 - coord[0] * coord[0];
    val[0] = (2 * M_PI * M_PI + 1.) * M_PI *
             cos(M_PI * coord[0]) * sin(M_PI * coord[1]);
    val[1] = (2 * M_PI * M_PI + 1.) * M_PI * (-1.) *
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

void
CreateCMesh(TPZCompMesh *&cmesh, TPZGeoMesh *gmesh, int pOrder, void (&func)(const TPZVec<REAL> &, TPZVec<STATE> &),
            const STATE &param, const std::string &prefix, const bool &print, const REAL &scale);

void CreateGMesh(TPZGeoMesh *&gmesh, const int meshType, const REAL hDomain, const REAL wDomain, const int xDiv,
                 const std::string &prefix, const bool &print, const REAL &scale);

void RunSimulation(const int &nDiv, const int &pOrder);


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    int pOrder = 1; //PARAMS
    const int nDivIni = 4; //PARAMS
    const int nPcycles = 3;
    const int nHcycles = 5;

    for (int iP = 0; iP < nPcycles; ++iP, ++pOrder) {
        int nDiv = nDivIni;
        for (int iH = 0; iH < nHcycles; ++iH) {
            std::cout << iH+1<<"/"<<nHcycles<<": Beginning simulation with nEl = "
                      << nDiv * nDiv * 2
                      <<" and p = "<<pOrder<< std::endl;
            RunSimulation(nDiv, pOrder);
            nDiv *= 2;
            std::cout << "************************************" << std::endl;
        }
    }
    return 0;
}

void RunSimulation(const int &nDiv, const int &pOrder) {
    // PARAMETROS FISICOS DO PROBLEMA
    const enum meshTypeE meshType = createTriangular;
    const REAL hDomain = 2;
    const REAL wDomain = 2;
    const int nThreads = 8; //PARAMS
    const bool l2error = true; //PARAMS
    const bool genVTK = false; //PARAMS
    const std::string prefix = "../";//PARAMS
    const bool print = true;//PARAMS
    const STATE param = 1.;//PARAMS
    const int postprocessRes = 0;//PARAMS
    const REAL scale = 1.;//PARAMS

    std::cout<<"Creating gmesh... ";
    boost::posix_time::ptime t1_g =
        boost::posix_time::microsec_clock::local_time();
    TPZGeoMesh *gmesh = NULL;
    CreateGMesh(gmesh, meshType, hDomain, wDomain, nDiv, prefix, print, scale);
    boost::posix_time::ptime t2_g =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_g-t1_g<<std::endl;

    std::cout<<"Creating cmesh... ";
    boost::posix_time::ptime t1_c =
        boost::posix_time::microsec_clock::local_time();
    TPZCompMesh *cmeshHCurl = NULL;
    CreateCMesh(cmeshHCurl, gmesh, pOrder, loadVec, param, prefix, print,
                scale); // funcao para criar a malha computacional
    boost::posix_time::ptime t2_c =
        boost::posix_time::microsec_clock::local_time();
    std::cout<<"Created! "<<t2_c-t1_c<<std::endl;

    TPZAnalysis an(cmeshHCurl);
    // configuracoes do objeto de analise
    TPZManVector<long, 1000> activeEquations;
    int neq = 0;
    int neqOriginal = 0;

//  TPZAutoPointer<TPZStructMatrix> strmtrx;
//  strmtrx = new TPZSpStructMatrix(cmeshHCurl);
//  strmtrx->SetNumThreads(nThreads);
//  FilterBoundaryEquations(cmeshHCurl, activeEquations, neq, neqOriginal);
//  strmtrx->EquationFilter().SetActiveEquations(activeEquations);
//  an.SetStructuralMatrix(strmtrx);

    TPZSymetricSpStructMatrix matrix(cmeshHCurl);
    //TPZSBandStructMatrix matrix(cmeshHCurl);
    matrix.SetNumThreads(nThreads);
    FilterBoundaryEquations(cmeshHCurl, activeEquations, neq, neqOriginal);
    matrix.EquationFilter().SetActiveEquations(activeEquations);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.SetStructuralMatrix(matrix);

    boost::posix_time::ptime t1_a =
        boost::posix_time::microsec_clock::local_time();
    an.Assemble();
    boost::posix_time::ptime t2_a =
        boost::posix_time::microsec_clock::local_time();
    boost::posix_time::ptime t1_s =
        boost::posix_time::microsec_clock::local_time();
    an.Solve();
    boost::posix_time::ptime t2_s =
        boost::posix_time::microsec_clock::local_time();


    std::cout<<"Assembly duration: "<<t2_a-t1_a
             <<"  Solver duration: "<<t2_s-t1_s<<std::endl;
    an.LoadSolution();
    if (l2error) {
        std::cout<<"Calculating errors..."<<std::endl;
        TPZVec<REAL> errorVec(1,0);
        an.SetExact(&exactSol);
        errorVec.Resize(2, 0.);
        an.SetThreadsForError(nThreads);
        an.PostProcessError(errorVec);
        std::string fileName = prefix + "errorHDiv";
        fileName.append(std::to_string(nDiv * nDiv * 2));
        fileName.append("_p");
        fileName.append(std::to_string(pOrder));
        fileName.append(".csv");
        std::ofstream errorFile(fileName.c_str());
        errorFile << errorVec[0] << "," << errorVec[1] << "," << errorVec[2]
                  << std::endl;
        errorFile.close();
        std::cout<<" Done!"<<std::endl;
    }
    if (genVTK) {
        std::cout<<"Post processing... ";
        TPZStack<std::string> scalnames, vecnames;
        vecnames.Push("E");
        std::string plotfile = prefix+"sol";
        plotfile.append(std::to_string(cmeshHCurl->NElements()));
        plotfile.append(".vtk");

        an.DefineGraphMesh(2, scalnames, vecnames,
                           plotfile);  // define malha grafica
        int postProcessResolution = postprocessRes; // define resolucao do pos processamento
        an.PostProcess(postProcessResolution);
        std::cout<<" Done!"<<std::endl;
    }
    gmesh->SetReference(nullptr);
    cmeshHCurl->SetReference(nullptr);
    delete cmeshHCurl;
    delete gmesh;
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
        TPZBndCond *mat = dynamic_cast<TPZBndCond *>(cmeshHCurl->MaterialVec()[cel->Reference()->MaterialId()]);
        if (mat && mat->Type() == 0) {
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

void CreateGMesh(TPZGeoMesh *&gmesh, const int meshType, const REAL hDomainOrig, const REAL wDomainOrig, const int nDiv,
                 const std::string &prefix, const bool &print, const REAL &scale) {
    const REAL wDomain = wDomainOrig/scale;
    const REAL hDomain = hDomainOrig/scale;

    TPZManVector<int, 3> nx(3, 0);
    TPZManVector<REAL, 3> llCoord(3, 0.), ulCoord(3, 0.), urCoord(3, 0.),lrCoord(3, 0.);
    llCoord[0] = -wDomain / 2;
    llCoord[1] = -hDomain / 2;

    ulCoord[0] = -wDomain / 2;
    ulCoord[1] = hDomain / 2;

    urCoord[0] = wDomain / 2;
    urCoord[1] = hDomain / 2;

    lrCoord[0] = wDomain / 2;
    lrCoord[1] = -hDomain / 2;

    nx[0] = nDiv;
    nx[1] = nDiv;
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

    gmesh->BuildConnectivity();
    if(print) {
        std::ofstream outTxt, outVtk;
        switch (meshType) {
            case createRectangular: {
                outTxt.open(prefix + "gmeshRec.txt"); // define arquivo de saida para
                // impressao da malha no
                outVtk.open("gmeshRec.vtk"); // define arquivo de saida para
                // impressao da malha no paraview

            }
                break;
            case createTriangular: {
                outTxt.open("gmeshTri.txt"); // define arquivo de saida para
                // impressao da malha no
                outVtk.open("gmeshTri.vtk"); // define arquivo de saida para
                // impressao da malha no paraview

            }
                break;
            case createZigZag: {
                outTxt.open("gmeshZig.txt"); // define arquivo de saida para
                // impressao da malha no
                outVtk.open("gmeshZig.vtk"); // define arquivo de saida para
                // impressao da malha no paraview
            }
                break;
            default:
                DebugStop();
                break;
        }
        gmesh->Print(outTxt);
        outTxt.close();
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVtk,
                                     true); // imprime a malha no formato vtk
        outVtk.close();
    }
    delete gengrid;
}

void CreateCMesh(TPZCompMesh *&cmeshHCurl, TPZGeoMesh *gmesh, int pOrder, void (&loadVec)(const TPZVec<REAL> &, TPZVec<STATE> &)
            , const STATE &param, const std::string &prefix, const bool &print, const REAL &scale) {
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
    TPZMatHelmholtz2D *matHCurl = NULL;
    matHCurl = new TPZMatHelmholtz2D(matId, param, scale);
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
    if(print){
        std::string fileName(prefix+ "cmesh");
        fileName.append(".txt");
        std::ofstream fileHCurl(fileName.c_str());
        cmeshHCurl->Print(fileHCurl);
    }
}
