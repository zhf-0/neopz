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
#include "TPZMatHCurlProjection.h"
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzlog.h"
#include "TPZTimer.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "pzsbstrmatrix.h"
#include "pzmatred.h"
#include "pzsubcmesh.h"
#include "tpzmatredstructmatrix.h"
#include "pzsbndmat.h"
#include "pzfstrmatrix.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzvisualmatrix.h"


enum meshTypeE{ createRectangular=1, createTriangular, createZigZag};

void desiredSol(const TPZVec<REAL> &coord, TPZVec<STATE> &val){
    val.Resize(3, 0.);
    REAL a = 9 * 2.54 * 1e-3;
    val[1] = 800 * sin(M_PI * (coord[0] + a/2) / a);
    
}

void FilterBoundaryEquations(TPZCompMesh * cmeshHCurl , TPZVec<long> &activeEquations , int &neq, int &neqOriginal);

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, void (& func)( const TPZVec<REAL> & , TPZVec<STATE> &) );

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int xDiv, const int zDiv);

int main(int argc, char *argv[])
{
    InitializePZLOG();
    
    
    HDivPiola = 1;//use piola mapping
    TPZTimer timer;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    //PARAMETROS FISICOS DO PROBLEMA
    const int dim = 2;
    REAL hDomain = 4 * 2.54 * 1e-3;
    REAL wDomain = 9 * 2.54 * 1e-3;

    int pOrder = 1; //ordem polinomial de aproximacao
    int xDiv = 12;
    int zDiv = 4;

    bool filterEquations = false;
    bool usingFullMtrx = true;
    bool optimizeBandwidth = false;
    
    const int meshType = createTriangular;
    timer.start();
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    CreateGMesh(gmesh, meshType, hDomain, wDomain,  xDiv, zDiv);
    
    TPZCompMesh *cmeshHCurl = CMesh(gmesh, pOrder, desiredSol); //funcao para criar a malha computacional
    
    TPZAnalysis an(cmeshHCurl,optimizeBandwidth);
    //configuracoes do objeto de analise
    TPZManVector<long,1000>activeEquations;
    int neq = 0;
    int neqOriginal = 0;
    FilterBoundaryEquations(cmeshHCurl, activeEquations , neq , neqOriginal);
    
    TPZAutoPointer<TPZSBandStructMatrix> sbstr;
    
    TPZAutoPointer<TPZFStructMatrix> fmtrx;
    
    if(usingFullMtrx){
        fmtrx = new TPZFStructMatrix(cmeshHCurl);
        fmtrx->SetNumThreads(0);
        if (filterEquations) {
            fmtrx->EquationFilter().SetActiveEquations(activeEquations);
        }
        
        an.SetStructuralMatrix(fmtrx);
    }
    else{
        sbstr = new TPZSBandStructMatrix(cmeshHCurl);
        sbstr->SetNumThreads(0);
        if (filterEquations) {
            sbstr->EquationFilter().SetActiveEquations(activeEquations);
        }
        an.SetStructuralMatrix(sbstr);
    }
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky); //caso simetrico
    an.SetSolver(step);
    
    an.Run();
    
    
    an.LoadSolution();
    TPZStack<std::string> scalnames, vecnames;
    vecnames.Push("Et");
    std::string plotfile= "../waveguideModes.vtk";//arquivo de saida que estara na pasta debug

    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
    int postProcessResolution = 2;//define resolucao do pos processamento
    an.PostProcess(postProcessResolution);
    
    
    return 0;
}

void FilterBoundaryEquations(TPZCompMesh *cmeshHCurl , TPZVec<long> &activeEquations , int &neq , int &neqOriginal)
{
    
    TPZManVector<long,1000> allConnects;
    std::set<long> boundConnects;
    
    
    for (int iel = 0; iel < cmeshHCurl->NElements(); iel++) {
        TPZCompEl *cel = cmeshHCurl->ElementVec()[iel];
        if ( cel == NULL) {
            continue;
        }
        if ( cel->Reference() == NULL) {
            
            continue;
        }
        if(cel->Reference()->MaterialId() == -1 ){
            std::set<long> boundConnectsEl;
            std::set<long> depBoundConnectsEl;
            std::set<long> indepBoundConnectsEl;
            cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
            cel->BuildConnectList(boundConnectsEl);
            
            for (std::set<long>::iterator iT = depBoundConnectsEl.begin(); iT != depBoundConnectsEl.end(); iT++) {
                const long val = *iT;
                //                std::cout<<"val"<<std::endl;
                //                std::cout<<val<<std::endl;
                //                std::cout<<"cmeshMF->ConnectVec()[val].FirstDepend()->fDepConnectIndex"<<std::endl;
                //                std::cout<<cmeshMF->ConnectVec()[val].FirstDepend()->fDepConnectIndex<<std::endl;
                
                if( boundConnects.find(cmeshHCurl->ConnectVec()[val].FirstDepend()->fDepConnectIndex) == boundConnects.end() ){
                    boundConnects.insert(cmeshHCurl->ConnectVec()[val].FirstDepend()->fDepConnectIndex);
                }
            }
        }
    }
    
    
    for (int iCon = 0; iCon < cmeshHCurl->NConnects(); iCon++) {
        if( boundConnects.find( iCon ) == boundConnects.end() ){
            TPZConnect &con = cmeshHCurl->ConnectVec()[iCon];
            if( con.HasDependency() ) continue;
            int seqnum = con.SequenceNumber();
            int pos = cmeshHCurl->Block().Position(seqnum);
            int blocksize = cmeshHCurl->Block().Size(seqnum);
            if(blocksize == 0) continue;
            
            int vs = activeEquations.size();
            activeEquations.Resize(vs+blocksize);
            for(int ieq = 0; ieq<blocksize; ieq++)
            {
                activeEquations[vs+ieq] = pos+ieq;
            }
        }
    }
    std::cout<<"------\t------\t-------"<<std::endl;
    std::cout<<"cmeshHCurl->NEquations()"<<"\t"<<cmeshHCurl->NEquations()<<std::endl;
    neqOriginal = cmeshHCurl->NEquations();
    //    std::cout<<"bound connects"<<std::endl;
    //    for (std::set<long>::iterator iT = boundConnects.begin(); iT != boundConnects.end(); iT++) {
    //        const long iCon = *iT;
    //        int seqnum = cmeshMF->ConnectVec()[iCon].SequenceNumber();
    //        int pos = cmeshMF->Block().Position(seqnum);
    //        int blocksize = cmeshMF->Block().Size(seqnum);
    //        std::cout<<"connect #"<<"\t";
    //        std::cout<<iCon<<"\t";
    //        std::cout<<"seq#"<<"\t";
    //        std::cout<<seqnum<<"\t";
    //        std::cout<<"pos"<<"\t";
    //        std::cout<<pos<<"\t";
    //        std::cout<<"block size"<<"\t";
    //        std::cout<<blocksize<<std::endl;
    //    }
    
    std::cout<<"activeEquations"<<std::endl;
    long nEq = 0;
    for (int iCon = 0; iCon < cmeshHCurl->NConnects(); iCon++) {
        if( boundConnects.find( iCon ) == boundConnects.end() ){
            if( cmeshHCurl->ConnectVec()[iCon].HasDependency() ) continue;
            int seqnum = cmeshHCurl->ConnectVec()[iCon].SequenceNumber();
            int blocksize = cmeshHCurl->Block().Size(seqnum);
            if(blocksize == 0) continue;
            nEq++;
        }
    }
    std::cout<<"# equations: "<<nEq<<std::endl;
    neq = nEq;
    return;
}

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int xDiv, const int yDiv)
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
    
    nx[0]=xDiv;
    nx[1]=yDiv;
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

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, void (& func)( const TPZVec<REAL> & , TPZVec<STATE> &) )
{
    
    const int dim = 2; //dimensao do problema
    const int matId = 1; //define id para um material(formulacao fraca)
    const int bc0 = -1; //define id para um material(cond contorno dirichlet)
    enum{ dirichlet = 0, neumann, mixed}; //tipo da condicao de contorno do problema
    // Criando material
    
    
    TPZCompMesh * cmeshHCurl = new TPZCompMesh(gmesh);
    cmeshHCurl->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmeshHCurl->SetDimModel(dim);//seta dimensao do modelo
    // Inserindo material na malha
    TPZMatHCurlProjection *matHCurl = new TPZMatHCurlProjection(matId);
    matHCurl->SetForcingFunction(func, 4);
    cmeshHCurl->InsertMaterialObject(matHCurl);
    
    TPZFMatrix<STATE> val1(1,1,0.) , val2(1,1,0.);
    val1( 0, 0 ) = 0.;
    val2( 0, 0 ) = 0.;
    TPZMaterial * BCondHCurlDir = matHCurl->CreateBC(matHCurl, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
    
    cmeshHCurl->InsertMaterialObject(BCondHCurlDir);//insere material na malha
    
    cmeshHCurl->SetAllCreateFunctionsHDiv();//define espaco de aproximacao
    cmeshHCurl->AutoBuild();
    
    
    TPZAdmChunkVector< TPZCompEl* > elVec = cmeshHCurl->ElementVec();
    
    for (int i = 0; i < cmeshHCurl->NElements(); i++) {
        TPZCompElHDiv < pzshape::TPZShapeQuad > *el = dynamic_cast<TPZCompElHDiv <pzshape::TPZShapeQuad > *>( elVec[i] );
        if ( el == NULL) {
            continue;
        }
        el->SetSideOrient(4,  1);
        el->SetSideOrient(5,  1);
        el->SetSideOrient(6, -1);
        el->SetSideOrient(7, -1);
    }
    
    for (int i = 0; i < cmeshHCurl->NElements(); i++) {
        TPZCompElHDiv < pzshape::TPZShapeTriang > *el = dynamic_cast<TPZCompElHDiv <pzshape::TPZShapeTriang > *>( elVec[i] );
        if ( el == NULL) {
            continue;
        }
        if ( i % 2 == 1) {
            el->SetSideOrient(3,  1);
            el->SetSideOrient(4, -1);
            el->SetSideOrient(5, -1);
        }
        else{
            el->SetSideOrient(3,  1);
            el->SetSideOrient(4,  1);
            el->SetSideOrient(5, -1);
        }
    }
    
    
    if (pOrder == 1) {
        //cmesh->CleanUpUnconnectedNodes();
        TPZCreateApproximationSpace::MakeRaviartThomas(*cmeshHCurl);
        cmeshHCurl->CleanUpUnconnectedNodes();
    }
    else{//for now only lowest order elements are avaliable
        DebugStop();
    }
    return cmeshHCurl;
}