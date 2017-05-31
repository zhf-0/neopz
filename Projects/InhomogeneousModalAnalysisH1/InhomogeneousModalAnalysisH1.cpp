/**
 * @file
 * @brief Modal analysis for inhomogeneous waveguides with coupled Ez and Hz fields
 * @details Contains spurious modes with must be dealt with
 * @author Francisco Orlandini
 * @since 2017
 */

#include "TPZMatInhomogeneousModalAnalysisH1.h"
#include "TPZTimer.h"//TPZTimer
#include "pzanalysis.h"//TPZAnalysis
#include "pzfstrmatrix.h"//TPZFStructMatrix
#include "pzsbstrmatrix.h"//TPZSBandStructMatrix
#include "pzsbndmat.h"//TPZSBMatrix
#include "pzstepsolver.h"//TPZStepSolver
#include "pzextractval.h"//TPZExtractVal
#include "pzbuildmultiphysicsmesh.h"//TPZBuildMultiphysicsMesh
#include "pzgengrid.h"//TPZGenGrid
#include "TPZVTKGeoMesh.h"//TPZVTKGeoMesh
#include "pzl2projection.h"//TPZL2Projection
#include "pzbndcond.h"//TPZBndCond

enum meshTypeE{ createRectangular=1, createTriangular, createZigZag};

STATE ur( const TPZVec<REAL> & x){
    return 1.;
}
STATE er( const TPZVec<REAL> & x){
    if( x[1] > 0 )
        return 1.;
    else
        return 4.;
}
void RunSimulation( const int meshType , int pOrder, int nDiv, REAL hDomain, REAL wDomain, REAL kzOverk02);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> cmeshMF , TPZVec<long> &activeEquations , int &neq, int &neqOriginal);

TPZVec<TPZCompMesh *>CreateCMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL kzOverk02);

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int nDiv);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    //PARAMETROS FISICOS DO PROBLEMA
    
    REAL wDomain = 1;
    REAL hDomain = wDomain/2;
    //wDomain = 0.5; //plane of simmetry as an electric wall
    
    const int meshType = createTriangular;
    int pOrder = 1; //ordem polinomial de aproximacao
    REAL kzOverk02 = 0.0;
    int nDiv = 15;
    int nSim = 10;
    for (int i = 0 ; i < nSim; i++) {
        std::cout<<"iteration "<<i+1<<" of "<<nSim<<std::endl;
        RunSimulation( meshType, pOrder, nDiv, hDomain, wDomain, kzOverk02);
        kzOverk02+=0.1;
    }
    
    
    
    return 0;
}

void RunSimulation( const int meshType , int pOrder, int nDiv, REAL hDomain, REAL wDomain, REAL kzOverk02){
    TPZTimer timer;
    
    timer.start();
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    CreateGMesh(gmesh, meshType, hDomain, wDomain,  nDiv);
    
    TPZVec<TPZCompMesh *>meshVec = CreateCMesh(gmesh, pOrder, ur , er , kzOverk02); //funcao para criar a malha computacional
    TPZCompMesh *cmeshMF = meshVec[0];
    TPZMatInhomogeneousModalAnalysisH1 *matPointer = dynamic_cast<TPZMatInhomogeneousModalAnalysisH1 * >( cmeshMF->MaterialVec()[1]);
    TPZVec<TPZCompMesh *>temporalMeshVec(2);
    temporalMeshVec[ matPointer->EzIndex() ] = meshVec[1 + matPointer->EzIndex()];
    temporalMeshVec[ matPointer->HzIndex() ] = meshVec[1 + matPointer->HzIndex()];
    
    bool optimizeBandwidth = true;
    TPZAnalysis an(cmeshMF,optimizeBandwidth);
    //configuracoes do objeto de analise
    TPZManVector<long,1000>activeEquations;
    int neq = 0;
    int neqOriginal = 0;
    FilterBoundaryEquations(meshVec, activeEquations , neq , neqOriginal);
    
    int nSolutions = neq >= 10 ? 10 : neq;
    
//    TPZAutoPointer<TPZSBandStructMatrix> strMtrx;
//    strMtrx = new TPZSBandStructMatrix(cmeshMF);
    TPZAutoPointer<TPZFStructMatrix> strMtrx;
    strMtrx = new TPZFStructMatrix(cmeshMF);
    
    strMtrx->EquationFilter().SetActiveEquations(activeEquations);
    an.SetStructuralMatrix(strMtrx);
  
    const int matId = 1;
    TPZMatInhomogeneousModalAnalysisH1 *matAlias = dynamic_cast<TPZMatInhomogeneousModalAnalysisH1 *> (cmeshMF->FindMaterial( matId ) );
    
    std::cout<<"entrando no assemble matrix A"<<std::endl;
    matAlias->SetMatrixA();
    
    an.Assemble();
    
//    TPZSBMatrix<STATE> *stiffAPtr = NULL , *stiffBPtr = NULL;
//    
//    stiffAPtr = new TPZSBMatrix<STATE> ( *dynamic_cast< TPZSBMatrix<STATE> *>( an.Solver().Matrix().operator->() ) );
//    matAlias->SetMatrixB();
//    std::cout<<"entrando no assemble matrix B"<<std::endl;
//    an.Assemble();
//    stiffBPtr = new TPZSBMatrix<STATE> ( *dynamic_cast< TPZSBMatrix<STATE> *>( an.Solver().Matrix().operator->() ) );
    TPZFMatrix<STATE> *stiffAPtr = NULL , *stiffBPtr = NULL;
    
    stiffAPtr = new TPZFMatrix<STATE> ( *dynamic_cast< TPZFMatrix<STATE> *>( an.Solver().Matrix().operator->() ) );
    matAlias->SetMatrixB();
    std::cout<<"entrando no assemble matrix B"<<std::endl;
    an.Assemble();
    stiffBPtr = new TPZFMatrix<STATE> ( *dynamic_cast< TPZFMatrix<STATE> *>( an.Solver().Matrix().operator->() ) );
    
    std::cout<<"saindo do assemble"<<std::endl;
    
#ifdef PZDEBUG
    std::ofstream fileA("../stiffA.csv");
    std::ofstream fileB("../stiffB.csv");
    char number[256];
    for (int i = 0; i<stiffAPtr->Rows(); i++) {
        for(int j = 0 ; j<stiffAPtr->Rows();j++){
            sprintf(number, "%16.16Lf",(long double) TPZExtractVal::val(std::real(stiffAPtr->GetVal(i,j))) );
            fileA<<number;
            
            fileA<<" + I * ";
            
            sprintf(number, "%16.16Lf",(long double) TPZExtractVal::val(std::imag(stiffAPtr->GetVal(i,j))) );
            fileA<<number;
            if( j != stiffAPtr->Rows() - 1){
                fileA<<" , ";
            }
            
            sprintf(number, "%16.16Lf",(long double) TPZExtractVal::val(std::real(stiffBPtr->GetVal(i,j))) );
            fileB<<number;
            
            fileB<<" + I * ";
            
            sprintf(number, "%16.16Lf",(long double) TPZExtractVal::val(std::imag(stiffBPtr->GetVal(i,j))) );
            fileB<<number;
            if( j != stiffAPtr->Rows() - 1){
                fileB<<" , ";
            }
        }
        fileA<<std::endl;
        fileB<<std::endl;
    }
    fileA.close();
    fileB.close();
#endif
    
    
    TPZVec<STATE> eValues;
    TPZFMatrix< STATE > eVectors;
    std::cout<<"entrando no calculo dos autovalores"<<std::endl;
    
    stiffAPtr->SolveGeneralisedEigenProblem( *stiffBPtr, eValues , eVectors);
    
    std::cout<<"saindo do calculo dos autovalores"<<std::endl;
    timer.stop();
    std::cout<<"salvando resultados"<<std::endl;
    std::set<std::pair<REAL,TPZFMatrix<STATE> > > eigenSolutions;//OTHERWISE IT WONT BE POSSIBLE TO USE ITERATORS
    TPZFMatrix<STATE> eVector( eVectors.Rows() , 1);
    std::pair<REAL,TPZFMatrix<STATE> > eigenPair;
    
    for(int i = 0 ; i< eValues.size();i++){
        eVectors.GetSub(0, i, eVectors.Rows(), 1 , eVector);
        eigenPair.first = eValues[i].real();
        eigenPair.second = eVector;
        eigenSolutions.insert(eigenPair);
    }
    
    std::string pathName;
    pathName = "..";
    
    std::string fileName = pathName;
    fileName.append("/ev");
    fileName.append(std::to_string(nDiv));
    fileName.append(".csv");
    
    std::ofstream fileEigenValues(fileName.c_str());
    fileName = pathName;
    fileName.append("/k0a");
    fileName.append(std::to_string(kzOverk02));
    fileName.append(".csv");
    std::ofstream filek0a(fileName.c_str());
    int i = 0;
    for (std::set<std::pair<REAL,TPZFMatrix<STATE> > > ::iterator iT = eigenSolutions.begin(); iT != eigenSolutions.end(); iT++) {
        fileEigenValues<<iT->first<<std::endl;
        if(iT->first > 0.01 && i <= nSolutions){
            filek0a<<sqrt(iT->first)*wDomain<<std::endl;
            std::cout<<sqrt(iT->first)*wDomain<<std::endl;
            i++;
        }
    }
//    std::cout << "Post Processing..." << std::endl;
//    
//    
//    TPZFMatrix<STATE> solMat(neqOriginal,1,0.);
//    
//    TPZStack<std::string> scalnames, vecnames;
//    scalnames.Push("Ez");
//    scalnames.Push("Hz");
//    std::string plotfile= "../waveguideModes.vtk";//arquivo de saida que estara na pasta debug
//    int dim = 2;
//    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
//    int postProcessResolution = 0;//define resolucao do pos processamento
//    
//    
//    std::set<std::pair<REAL,TPZFMatrix<STATE> > > ::iterator iT = eigenSolutions.begin();
//    for (int iSol = 0; iSol < nSolutions; iSol ++) {
//        if( iT == eigenSolutions.end() ){
//            DebugStop();
//        }
//        for (int i = 0 ; i < neq; i++) {
//            solMat( activeEquations[i] ,0) = (iT->second).GetVal(i, 0);
//        }
//        an.LoadSolution( solMat );
//        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(temporalMeshVec, cmeshMF);
//        an.PostProcess(postProcessResolution);
//        iT++;
//    }
    std::cout << "FINISHED!" << std::endl;
}

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> meshVec , TPZVec<long> &activeEquations , int &neq, int &neqOriginal)
{
    TPZCompMesh *cmeshMF = meshVec[0];
    
    const int matId = 1;
    TPZMatInhomogeneousModalAnalysisH1 *matAlias = dynamic_cast<TPZMatInhomogeneousModalAnalysisH1 *> (cmeshMF->FindMaterial( matId ) );
    TPZCompMesh *cmeshHZ = meshVec[1 + matAlias->HzIndex()];
    TPZCompMesh *cmeshEZ = meshVec[1 + matAlias->EzIndex()];
    
    TPZManVector<long,1000> allConnects;
    std::set<long> boundConnects;
    
    
    for (int iel = 0; iel < cmeshMF->NElements(); iel++) {
        TPZCompEl *cel = cmeshMF->ElementVec()[iel];
        if ( cel == NULL) {
            continue;
        }
        if ( cel->Reference() == NULL) {
            continue;
        }
        if(cel->Reference()->MaterialId() == -1 ){//boundaryEl
            std::set<long> boundConnectsEl;
            std::set<long> depBoundConnectsEl;
            std::set<long> indepBoundConnectsEl;
            cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
            cel->BuildConnectList(boundConnectsEl);
            for (std::set<long>::iterator iT = boundConnectsEl.begin(); iT != boundConnectsEl.end(); iT++) {
                const long val = *iT;
                //must select only EZ connects
                if( matAlias->HzIndex() == 0 && val >= cmeshHZ->NConnects() ){
                    if( boundConnects.find(val) == boundConnects.end() ){
                        boundConnects.insert(val);
                    }
                }
                else if( matAlias->HzIndex() == 1 && val < cmeshEZ->NConnects() ){
                    if( boundConnects.find(val) == boundConnects.end() ){
                        boundConnects.insert(val);
                    }
                }
            }
        }
    }
    //add active equations to vec
    for (int iCon = 0; iCon < cmeshMF->NConnects(); iCon++) {
        if( boundConnects.find( iCon ) == boundConnects.end() ){
            TPZConnect &con = cmeshMF->ConnectVec()[iCon];
            int seqnum = con.SequenceNumber();
            int pos = cmeshMF->Block().Position(seqnum);
            int blocksize = cmeshMF->Block().Size(seqnum);
            if(blocksize == 0) continue;
            
            int vs = activeEquations.size();
            activeEquations.Resize(vs+blocksize);
            for(int ieq = 0; ieq<blocksize; ieq++)
            {
                activeEquations[vs+ieq] = pos+ieq;
            }
        }
    }
    
    std::cout<<"------ ---Before Removal of Boundary Connects--- -------"<<std::endl;
    std::cout<<"cmeshEZ->NEquations()"<<"\t"<<cmeshEZ->NEquations()<<std::endl;
    std::cout<<"cmeshHHZ->NEquations()"<<"\t"<<cmeshHZ->NEquations()<<std::endl;
    std::cout<<"cmeshMF->NEquations()"<<"\t"<<cmeshMF->NEquations()<<std::endl;
    neqOriginal = cmeshMF->NEquations();
    
    int nHZEquations = 0 , nEZEquations = 0;
    std::cout<<"------ ---After Removal of Boundary Connects--- -------"<<std::endl;
    long nEq = 0;
    for (int iCon = 0; iCon < cmeshMF->NConnects(); iCon++) {
        bool isEZ;
        if( boundConnects.find( iCon ) == boundConnects.end() ){
            if( cmeshMF->ConnectVec()[iCon].HasDependency() ) continue;
            int seqnum = cmeshMF->ConnectVec()[iCon].SequenceNumber();
            int blocksize = cmeshMF->Block().Size(seqnum);
            if( matAlias->HzIndex() == 0 && iCon < cmeshHZ->NConnects() ){
                isEZ = false;
            }
            else if( matAlias->HzIndex() == 1 && iCon >= cmeshEZ->NConnects() ){
                isEZ = false;
            }
            else{
                isEZ = true;
            }
            for(int ieq = 0; ieq<blocksize; ieq++)
            {
                nEq++;
                isEZ == true ? nEZEquations ++ : nHZEquations++;
            }
        }
    }
    std::cout<<"cmeshEZ->NEquations()"<<"\t"<< nEZEquations<<std::endl;
    std::cout<<"cmeshHZ->NEquations()"<<"\t"<<nHZEquations<<std::endl;
    std::cout<<"cmeshMF->NEquations()"<<"\t"<<nEq<<std::endl;
    //std::cout<<activeEquations<<std::endl;
    neq = nEq;
    return;
}

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int nDiv)
{
    
    TPZManVector<int,3> nx(3,0);
    TPZManVector<REAL,3> llCoord(3,0.) , ulCoord(3,0.) , urCoord(3,0.) , lrCoord(3,0.);
    llCoord[0] = -wDomain/2;
    llCoord[1] = -hDomain/2;
    
    ulCoord[0] = -wDomain/2;
    ulCoord[1] = hDomain/2;
    
    urCoord[0] = wDomain/2;
    urCoord[1] = hDomain/2;
    
    lrCoord[0] = wDomain/2;
    lrCoord[1] = -hDomain/2;
    
    nx[0]=nDiv;
    nx[1]=nDiv;
    int numl = 1;
    TPZGenGrid *gengrid = NULL;
    switch (meshType) {
        case createRectangular:
        {
            REAL rot = 0.0;
            gengrid = new TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
            gengrid->SetElementType(EQuadrilateral);
        }
            break;
        case createTriangular:
        {
            REAL rot = 0.0;
            gengrid = new TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
            gengrid->SetElementType(ETriangle);
        }
            break;
        case createZigZag:
        {
            REAL rot = 0.0;
            gengrid = new TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
            gengrid->SetElementType(EQuadrilateral);
            gengrid->SetZigZagPattern();
        }
            break;
        default:
            DebugStop();
            break;
    }
    gmesh = new TPZGeoMesh();
    const int matId = 1; //define id para um material(formulacao fraca)
    const int bc0 = -1; //define id para um material(cond contorno dirichlet)
    gengrid->Read(gmesh , matId);
    
    
    gengrid->SetBC(gmesh, ulCoord, llCoord, bc0);
    gengrid->SetBC(gmesh, urCoord, ulCoord, bc0);
    gengrid->SetBC(gmesh, lrCoord, urCoord, bc0);
    gengrid->SetBC(gmesh, llCoord, lrCoord, bc0);
    
    //gmesh->ResetConnectivities();
    
    gmesh->BuildConnectivity();
#ifdef PZDEBUG
    std::ofstream outTxt , outVtk;
    switch (meshType) {
        case createRectangular:
        {
            outTxt.open("../gmeshRectangular.txt"); //define arquivo de saida para impressao da malha no
            outVtk.open("../gmeshRectangular.vtk"); //define arquivo de saida para impressao da malha no paraview
            
        }
            break;
        case createTriangular:
        {
            outTxt.open("../gmeshTriangular.txt"); //define arquivo de saida para impressao da malha no
            outVtk.open("../gmeshTriangular.vtk"); //define arquivo de saida para impressao da malha no paraview
            
        }
            break;
        case createZigZag:
        {
            outTxt.open("../gmeshZigZag.txt"); //define arquivo de saida para impressao da malha no
            outVtk.open("../gmeshZigZag.vtk"); //define arquivo de saida para impressao da malha no paraview
        }
            break;
        default:
            DebugStop();
            break;
    }
    gmesh->Print(outTxt);
    outTxt.close();
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVtk, true); //imprime a malha no formato vtk
    outVtk.close();
#endif
    delete gengrid;
}

TPZVec<TPZCompMesh *>CreateCMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL kzOverk02)
{
    
    const int dim = 2; //dimensao do problema
    const int matId = 1; //define id para um material(formulacao fraca)
    const int bc0 = -1; //define id para um material(cond contorno dirichlet)
    enum{ dirichlet = 0, neumann, mixed}; //tipo da condicao de contorno do problema
    // Criando material
    
    
    ///criar malha computacional HZ
    
    TPZCompMesh * cmeshHz = new TPZCompMesh(gmesh);
    cmeshHz->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmeshHz->SetDimModel(dim);//seta dimensao do modelo
    // Inserindo material na malha
    const int nState = 1;
    TPZVec<STATE> sol;//only for creating material. this material will not be used in reality
    TPZL2Projection *matHz = new TPZL2Projection( matId, dim, nState, sol);
    cmeshHz->InsertMaterialObject(matHz);
    
    
    TPZFNMatrix<1,STATE> val1(1,2,0.), val2(1,2,0.);
    TPZMaterial * BCondHZNeumann = matHz->CreateBC(matHz, bc0, neumann, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
    
    cmeshHz->InsertMaterialObject(BCondHZNeumann);//insere material na malha
    cmeshHz->SetAllCreateFunctionsContinuous();//define espaço de aproximação
    cmeshHz->AutoBuild();
    cmeshHz->CleanUpUnconnectedNodes();
    
    ///criar malha computacional EZ
    
    TPZCompMesh * cmeshEz = new TPZCompMesh(gmesh);
    cmeshEz->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmeshEz->SetDimModel(dim);//seta dimensao do modelo
    // Inserindo material na malha
    TPZL2Projection *matEz = new TPZL2Projection( matId, dim, nState, sol);
    cmeshEz->InsertMaterialObject(matEz);
    
    TPZMaterial * BCondEzDir = matEz->CreateBC(matEz, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
    
    cmeshEz->InsertMaterialObject(BCondEzDir);//insere material na malha
    cmeshEz->SetAllCreateFunctionsContinuous();//define espaço de aproximação
    cmeshEz->AutoBuild();
    cmeshEz->CleanUpUnconnectedNodes();
    
    TPZMatInhomogeneousModalAnalysisH1 *matMultiPhysics = new TPZMatInhomogeneousModalAnalysisH1(matId , kzOverk02 , ur , er);
    TPZVec<TPZCompMesh *> meshVec(2);
    
    meshVec[ matMultiPhysics->HzIndex() ] = cmeshHz;
    meshVec[ matMultiPhysics->EzIndex() ] = cmeshEz;
    
    
    TPZCompMesh *cmeshMF = new TPZCompMesh( gmesh );
    
    TPZMaterial * BCondMFDir = matMultiPhysics->CreateBC(matMultiPhysics, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
    
    
    cmeshMF->InsertMaterialObject(matMultiPhysics);
    cmeshMF->InsertMaterialObject(BCondMFDir);//insere material na malha
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
    std::ofstream fileHz("../cmeshHz.txt");
    cmeshHz->Print(fileHz);
    std::ofstream fileEz("../cmeshEz.txt");
    cmeshEz->Print(fileEz);
    std::ofstream fileMF("../cmeshMFHCurl.txt");
    cmeshMF->Print(fileMF);
    
    TPZVec< TPZCompMesh *>meshVecOut(3);
    
    meshVecOut[0] = cmeshMF;
    meshVecOut[1 + matMultiPhysics->HzIndex()] = cmeshHz;
    meshVecOut[1 + matMultiPhysics->EzIndex()] = cmeshEz;
    return meshVecOut;
}
