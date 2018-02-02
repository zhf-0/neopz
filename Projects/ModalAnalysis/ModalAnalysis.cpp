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
#ifdef USING_SLEPC
#include <TPZSlepcEPSHandler.h>
#elif defined USING_LAPACK
#include <TPZLapackWrapper.h>
#endif
#include "TPZRefPatternDataBase.h"
#include "tpzgeoelrefpattern.h"
#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "TPZMatWaveguideCutOffAnalysis.h"
#include "TPZMatModalAnalysis.h"
#include "TPZMatMFHDivRotH1.h"
#include "TPZGmshReader.h"

enum meshTypeE { createRectangular = 1, createTriangular, createZigZag };

void RunSimulation(bool isRectangularWG, bool isCutOff, const meshTypeE meshType, int pOrder,
                   int nDiv, const TPZVec<REAL> geoParams, REAL f0, bool genVTK,
                   bool l2error, bool exportEigen, const int nThreads,
                   bool optimizeBandwidth, bool filterEquations);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> cmeshMF,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);

void CreateCMesh(TPZVec<TPZCompMesh *> &meshVecOut, TPZGeoMesh *gmesh,
                 int pOrder, const STATE &ur,
                 const STATE &er, REAL f0, bool isCutOff);

void CreateGMeshRectangularWaveguide(TPZGeoMesh *&gmesh,
                                     const meshTypeE meshType,
                                     const REAL wDomain, const REAL hDomain, 
                                     const int nDiv);

void CreateGMeshCircularWaveguide(TPZGeoMesh *&gmesh, const meshTypeE meshType,
                                  const REAL rDomain, const int nDiv);

const bool usingGMSH = false;
const bool usingHDivRot = false;

int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
	
	  const bool isRectangularWG = true;//true = rectangular , false = circular
    const bool isCutOff = false;//analysis of cutoff frequencies for eigenmodes
    const meshTypeE meshType = createTriangular;
    int pOrder = 1;           // polynomial order of basis functions
    bool genVTK = false;      // generate vtk for fields visualisation
    bool l2error = false;     // TODO: implement error analysis
    bool exportEigen = false; // export eigen values
    const int nThreads = 8;
    bool optimizeBandwidth = true; //whether to renumber equations (OFF for debugging purposes)
    bool filterEquations = true; //whether to impose dirichlet conditions removing boundary equations
	
	TPZManVector<REAL, 2> geoParams(1,-1);
    REAL fOp = -1;//operational frequency
    int nDiv = -1;//number of mesh divisions
	if (isRectangularWG) { //WR-90 waveguide
		geoParams.Resize(2, 0.);
		geoParams[0] = 9 * 2.54 * 1e-3;//width
		geoParams[1] = 4 * 2.54 * 1e-3;//height
        fOp = 25e+9;
      nDiv = 25;
	}
	else{
		geoParams.Resize(1, 0.);
        geoParams[0] = 1.;//radius
        fOp = 9e+9;
      nDiv = 2;
//        geoParams[0] = 1.00;//radius
//        fOp = 25e+9;
	}

    const int nSim = 1;
    for (int i = 0; i < nSim; i++) {
        std::cout << "iteration " << i + 1 << " of " << nSim << std::endl;
        RunSimulation(isRectangularWG, isCutOff, meshType, pOrder, nDiv, geoParams, fOp,
                      genVTK, l2error, exportEigen, nThreads, optimizeBandwidth,
                      filterEquations);
        nDiv *= 2;
    }

    return 0;
}

void RunSimulation(bool isRectangularWG, bool isCutOff, const meshTypeE meshType, int pOrder,
                   int nDiv, const TPZVec<REAL> geoParams, REAL f0, bool genVTK,
                   bool l2error, bool exportEigen, const int nThreads,
                   bool optimizeBandwidth, bool filterEquations) {
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    std::cout<<"Creating GMesh...";
    #ifdef USING_BOOST
    boost::posix_time::ptime t1_g =
            boost::posix_time::microsec_clock::local_time();
    #endif
	if (isRectangularWG) {
		CreateGMeshRectangularWaveguide(gmesh, meshType, geoParams[0], geoParams[1], nDiv);
	}
	else{
		CreateGMeshCircularWaveguide(gmesh, meshType, geoParams[0], nDiv);
	}
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
    const STATE ur = 1.0;
    const STATE er = 1.0;
    CreateCMesh(meshVec, gmesh, pOrder, ur, er, f0,
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
    #ifdef USING_SLEPC
    TPZSlepcEPSHandler<STATE> solver;
    #elif defined USING_LAPACK
    TPZLapackWrapper<STATE> solver;
    #endif
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

void CreateGMeshRectangularWaveguide(TPZGeoMesh *&gmesh,
                                     const meshTypeE meshType,
                                     const REAL wDomain, const REAL hDomain,
                                     const int nDiv) {

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

    nx[0] = nDiv + 1;
    nx[1] = nDiv + 1;
//    nx[0] = 2;//REFINEMENT TEST
//    nx[1] = 1;//REFINEMENT TEST
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
    
//    //REFINEMENT TEST                                                       //
//    TPZManVector<TPZGeoEl *, 3> sons;                                       //
//    TPZManVector<REAL,3> qsi(3,0.), x(3,0.);                                //
//    qsi[0]= 0.5;                                                            //
//    qsi[1]= 0.5;                                                            //
//    TPZManVector<REAL,3> refPointsX(2,0.);                                  //
//    TPZManVector<REAL,2> refPointsY(2,0.);                                  //
//    refPointsX[0] = wDomain/2;                                              //
//    refPointsX[1] = wDomain/2;                                              //
//    refPointsY[0] = 0;                                                      //
//    refPointsY[1] = hDomain/2;                                              //
//    for (int iref = 0; iref < 2; iref++) {                                  //
//        int nel = gmesh->NElements();                                       //
//        for (int iel = 0; iel < nel; iel++) {                               //
//            TPZGeoEl *gel = gmesh->ElementVec()[iel];                       //
//            gel->X(qsi, x);//gets center of element                         //
//                                                                            //
//            if(x[0]>refPointsX[iref] && x[1] > refPointsY[iref]){           //
//                if (gel->HasSubElement()) {                                 //
//                    continue;                                               //
//                }                                                           //
//                gel->Divide(sons);                                          //
//            }                                                               //
//        }                                                                   //
//    }                                                                       //
    
#ifdef PZDEBUG
    std::ofstream outTxt, outVtk;
	std::string meshFileName("gmesh");
    switch (meshType) {
    case createRectangular: {
		meshFileName.append("Rectangular");
    } break;
    case createTriangular: {
        meshFileName.append("Triangular");
    } break;
    case createZigZag: {
        meshFileName.append("ZigZag");
	} break;
    default:
        DebugStop();
        break;
    }
	meshFileName.append(std::to_string(nDiv));
	const size_t strlen = meshFileName.length();
	meshFileName.append(".vtk");
	outVtk.open(meshFileName.c_str());
	meshFileName.replace(strlen, 4, ".txt");
	outTxt.open(meshFileName.c_str());
	
    gmesh->Print(outTxt);
    outTxt.close();
#ifdef PZDEBUG
	TPZCheckGeom * Geometrytest = new TPZCheckGeom(gmesh);
	int isBadMeshQ = Geometrytest->PerformCheck();
	
	if (isBadMeshQ) {
		DebugStop();
	}
#endif
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVtk,
                                 true); //prints vtk with geomesh
    outVtk.close();
#endif
    delete gengrid;
}

void CreateGMeshCircularWaveguide(TPZGeoMesh *&gmesh, const meshTypeE meshType,
                                  const REAL rDomain, const int nDiv) {
    if(usingGMSH == false){
        if (meshType != meshTypeE::createTriangular) {
            DebugStop(); // HCurl quadrilateral elements not implemented!
        }
        
        const int nNodes = 4 + 4 + 1;
        const int matId = 1; // define id para um material(formulacao fraca)
        const int bc0 = -1;  // define id para um material(cond contorno dirichlet)
        
        gmesh = new TPZGeoMesh;
        gmesh->NodeVec().Resize(nNodes);
        
        // create 4 nodes which will be triangle vertices
        // r arg 0, r arg pi/2, r arg pi, r arg 3pi/2
        int nodeId = 0;
        for (int iNode = 0; iNode < 4; iNode++) {
            TPZManVector<REAL, 3> node(3, 0.);
            const int c0 = (1+(iNode/2)*(-2))*((iNode+1)%2);//expected: 1 0 -1 0
            const int c1 = (1+((iNode-1)/2)*(-2))*(iNode%2);//expected: 0 1 0 -1
            node[0] = c0*rDomain;
            node[1] = c1*rDomain;
            gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
            nodeId++;
        }
        // create 4 nodes which will be triangle midsides points @ boundary
        // r arg pi/4, r arg 3pi/4, r arg 5pi/4, r arg 7pi/4
        for (int iNode = 0; iNode < 4; iNode++) {
            TPZManVector<REAL, 3> node(3, M_SQRT1_2 * rDomain);
            const int c0 = (1+(iNode/2)*(-2))*((iNode+1)%2);
            const int c1 = (1+((iNode-1)/2)*(-2))*(iNode%2);
            node[0] *= c0 - c1;//expected: +1 -1 -1 +1
            node[1] *= c0 + c1;//expected: +1 +1 -1 -1
            node[2] = 0.;
            gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
            nodeId++;
        }
        // create center node
        {
            TPZManVector<REAL, 3> node(3, 0.);
            gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
        }

        TPZManVector<long, 3> nodesIdVec(3, 0);
        // creates volumetric elements
        for (int iTri = 0; iTri < 4; iTri++) {
            nodesIdVec[0] = (iTri) % 4;
            nodesIdVec[1] = (iTri + 1) % 4;
            nodesIdVec[2] = 8;
          TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> > * triangulo =
                  new TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > > (nodesIdVec,matId,*gmesh);
        }
        // creates boundary elements
        for (int iArc = 0; iArc < 4; iArc++) {
            nodesIdVec[0] = (iArc) % 4;
            nodesIdVec[1] = (iArc + 1) % 4;
            nodesIdVec[2] = iArc + 4;//midpoint
          TPZGeoElRefPattern< pzgeom::TPZArc3D > *arc =
                  new TPZGeoElRefPattern< pzgeom::TPZArc3D > (nodesIdVec,bc0,*gmesh);
        }
        
        gmesh->BuildConnectivity();

        TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
        qsi[0]=0.49;
        qsi[1]=0.49;
        gmesh->ElementVec()[0]->X(qsi,xqsi);
        std::cout << "qsi = { ";
        std::cout << qsi[0] << " , " << qsi[1];
        std::cout << " }  ;  x(qsi) = { ";
        std::cout << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";

      TPZVec<TPZGeoEl *> sons;

      const int nref = nDiv;
      for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
          TPZGeoEl *gel = gmesh->ElementVec()[iel];
          if(!gel->HasSubElement()) gel->Divide(sons);
        }
      }
    }
    else{
        TPZGmshReader meshReader;
        gmesh = meshReader.GeometricGmshMesh("circlequad.msh");
    }
#ifdef PZDEBUG
//	TPZCheckGeom * Geometrytest = new TPZCheckGeom(gmesh);
//	int isBadMeshQ = Geometrytest->PerformCheck();
//
//	if (isBadMeshQ) {
//		DebugStop();
//	}
#endif

    std::string meshFileName("gmeshCirc");
    meshFileName.append(std::to_string(nDiv));
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
                 int pOrder, const STATE &ur,
                 const STATE &er, REAL f0, bool isCutOff) {

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
    if(!usingHDivRot) cmeshH1->SetDefaultOrder(pOrder); // seta ordem polimonial de aproximacao
    else cmeshH1->SetDefaultOrder(pOrder+1);
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
    TPZVecL2 *matHCurl = new TPZVecL2(matId);
    cmeshHCurl->InsertMaterialObject(matHCurl);

    val1(0, 0) = 0.;
    val2(0, 0) = 0.;
    TPZMaterial *BCondHCurlDir =
        matHCurl->CreateBC(matHCurl, bc0, dirichlet, val1,
                           val2); // cria material que implementa a condicao de
                                  // contorno de dirichlet

    cmeshHCurl->InsertMaterialObject(BCondHCurlDir); // insere material na malha
    
    if(!usingHDivRot)cmeshHCurl->SetAllCreateFunctionsHCurl(); // define espaco de aproximacao
    else cmeshHCurl->SetAllCreateFunctionsHDiv();
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
        if(!usingHDivRot) dummy = new TPZMatModalAnalysis(matId, f0, ur, er);
        else dummy = new TPZMatMFHDivRotH1(matId, f0, ur, er);
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

    //    std::cout << "------\t----------\t-------" << std::endl;
    //    std::cout << "cmeshH1->NEquations()"
    //              << "\t" << cmeshH1->NEquations() << std::endl;
    //    std::cout << "cmeshHCurl->NEquations()"
    //              << "\t" << cmeshHCurl->NEquations() << std::endl;
    //    std::cout << "cmeshMF->NEquations()"
    //              << "\t" << cmeshMF->NEquations() << std::endl;
    //    std::cout << "------\t----------\t-------" << std::endl;
    return;
}
