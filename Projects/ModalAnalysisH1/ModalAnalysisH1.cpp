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
#include "TPZMatModalAnalysisH1.h"
#include "pzextractval.h"
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzlog.h"
#include "TPZTimer.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzsbstrmatrix.h"
#include "pzl2projection.h"
#include "pzsbndmat.h"
#include "pzfstrmatrix.h"
#include "pzvisualmatrix.h"


enum meshTypeE{ createRectangular=1, createTriangular, createZigZag};

STATE ur( const TPZVec<REAL> & x){
    //return ( 2.-imaginary*0.1 );
    return 1.;
}
STATE er( const TPZVec<REAL> & x){
    return 1;
}
void RunSimulation( bool filterEquations, const int meshType ,bool usingFullMtrx, bool optimizeBandwidth, int pOrder, int nDiv, REAL hDomain, REAL wDomain, REAL f0, modeType teortm);

void FilterBoundaryEquations(TPZCompMesh * cMesh , TPZVec<long> &activeEquations , int &neq, int &neqOriginal);

TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL f0);

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int xDiv, const int zDiv);

int main(int argc, char *argv[])
{
    
    HDivPiola = 1;//use piola mapping
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    //PARAMETROS FISICOS DO PROBLEMA
    REAL hDomain = 4 * 2.54 * 1e-3;
    REAL wDomain = 9 * 2.54 * 1e-3;
    const modeType teortm = modesTE;
    REAL f0 = 15 * 1e+9;
    int pOrder = 1; //ordem polinomial de aproximacao
    

    bool filterEquations = true;
    bool usingFullMtrx = true;
    bool optimizeBandwidth = false;
    
    const int meshType = createTriangular;
    
    int nDiv = 10;
    int nSim = 1;
    for (int i = 0 ; i < nSim; i++) {
        std::cout<<"iteration "<<i+1<<" of "<<nSim<<std::endl;
        RunSimulation( filterEquations, meshType, usingFullMtrx, optimizeBandwidth, pOrder, nDiv, hDomain, wDomain, f0 , teortm);
        nDiv += 5;
    }
    
    
    
    return 0;
}

void RunSimulation( bool filterEquations, const int meshType ,bool usingFullMtrx, bool optimizeBandwidth, int pOrder, int nDiv, REAL hDomain, REAL wDomain, REAL f0, modeType teortm){
    TPZTimer timer;
    
    timer.start();
    
    int xDiv = nDiv;
    int zDiv = nDiv;
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    CreateGMesh(gmesh, meshType, hDomain, wDomain,  xDiv, zDiv);
    
    TPZCompMesh *cMesh = CreateCMesh(gmesh, pOrder, ur , er , f0); //funcao para criar a malha computacional
    
    
    
    TPZAnalysis an(cMesh,optimizeBandwidth);
    //configuracoes do objeto de analise
    TPZManVector<long,1000>activeEquations;
    int neq = 0;
    int neqOriginal = 0;
    FilterBoundaryEquations( cMesh, activeEquations , neq , neqOriginal);
    
    int nSolutions = neq >= 10 ? 10 : neq;
    
    TPZAutoPointer<TPZSBandStructMatrix> sbstr;
    
    TPZAutoPointer<TPZFStructMatrix> fmtrx;
    
    if(usingFullMtrx){
        fmtrx = new TPZFStructMatrix(cMesh);
        fmtrx->SetNumThreads(0);
        if (filterEquations) {
            fmtrx->EquationFilter().SetActiveEquations(activeEquations);
        }
        an.SetStructuralMatrix(fmtrx);
    }
    else{
        sbstr = new TPZSBandStructMatrix(cMesh);
        sbstr->SetNumThreads(0);
        if (filterEquations) {
            sbstr->EquationFilter().SetActiveEquations(activeEquations);
        }
        an.SetStructuralMatrix(sbstr);
    }
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky); //caso simetrico
    an.SetSolver(step);
    int matId = 1;
    TPZMatModalAnalysisH1 *matAlias = dynamic_cast<TPZMatModalAnalysisH1 *> (cMesh->FindMaterial( matId ) );
    
    
    matAlias->SetTEOrTM( teortm);
    std::cout<<"entrando no assemble matrix A"<<std::endl;
    matAlias->SetMatrixA();
    
    an.Assemble();
    TPZMatrix<STATE> *  stiffAPtr = NULL , *stiffBPtr = NULL;
    if (usingFullMtrx) {
        stiffAPtr = new TPZFMatrix<STATE> ( *dynamic_cast< TPZFMatrix<STATE> *>( an.Solver().Matrix().operator->() ) );
        matAlias->SetMatrixB();
        std::cout<<"entrando no assemble matrix B"<<std::endl;
        an.Assemble();
        stiffBPtr = new TPZFMatrix<STATE> ( *dynamic_cast< TPZFMatrix<STATE> *>( an.Solver().Matrix().operator->() ) );
        std::cout<<"saindo do assemble"<<std::endl;
        
        
    }
    else{
        stiffAPtr = new TPZSBMatrix<STATE> ( *dynamic_cast< TPZSBMatrix<STATE> *>( an.Solver().Matrix().operator->() ) );
        matAlias->SetMatrixB();
        std::cout<<"entrando no assemble matrix B"<<std::endl;
        an.Assemble();
        stiffBPtr = new TPZSBMatrix<STATE> ( *dynamic_cast< TPZSBMatrix<STATE> *>( an.Solver().Matrix().operator->() ) );
        std::cout<<"saindo do assemble"<<std::endl;
    }
    
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
    if (usingFullMtrx) {
        TPZFMatrix<STATE> *stiffA = dynamic_cast<TPZFMatrix<STATE> *>(stiffAPtr);
        TPZFMatrix<STATE> *stiffB = dynamic_cast<TPZFMatrix<STATE> *>(stiffBPtr);
        stiffA->SolveGeneralisedEigenProblem( *stiffB, eValues , eVectors);
    }
    else{
        TPZSBMatrix<STATE> *stiffA = dynamic_cast<TPZSBMatrix<STATE> *>(stiffAPtr);
        TPZSBMatrix<STATE> *stiffB = dynamic_cast<TPZSBMatrix<STATE> *>(stiffBPtr);
        stiffA->SolveGeneralisedEigenProblem( *stiffB, eValues , eVectors);
    }
    
    std::cout<<"saindo do calculo dos autovalores"<<std::endl;
    timer.stop();
    
    std::set<std::pair<REAL,TPZFMatrix<STATE> > > eigenValuesRe;
    TPZFMatrix<STATE> eVector( eVectors.Rows() , 1);
    std::pair<REAL,TPZFMatrix<STATE> > duplet;
    
    for(int i = 0 ; i< eValues.size();i++){
        eVectors.GetSub(0, i, eVectors.Rows(), 1 , eVector);
        duplet.first = eValues[i].real();
        duplet.second = eVector;
        eigenValuesRe.insert(duplet);
    }
    int i = 0;
    std::string fileName("../ev");
    fileName.append(std::to_string(nDiv));
    fileName.append(".txt");
    std::ofstream fileEigenValues(fileName.c_str());
    for (std::set<std::pair<REAL,TPZFMatrix<STATE> > > ::iterator iT = eigenValuesRe.begin(); iT != eigenValuesRe.end(); iT++) {
        std::cout<< iT->first <<std::endl;
        i++;
        if( i >= nSolutions)
            break;
        
    }

    std::cout << "Post Processing..." << std::endl;
    
    
    {
        std::ofstream fileA("../EV.csv");
        char number[256];
        std::cout<<eigenValuesRe.begin()->first<<std::endl;
        for (int i = 0; i<eigenValuesRe.begin()->second.Rows(); i++) {
            for(int j = 0 ; j<eigenValuesRe.begin()->second.Cols();j++){
                sprintf(number, "%32.32Lf",(long double) TPZExtractVal::val(std::real(eigenValuesRe.begin()->second.GetVal(i,j))) );
                fileA<<number;
                
                fileA<<" + I * ";
                
                sprintf(number, "%32.32Lf",(long double) TPZExtractVal::val(std::imag(eigenValuesRe.begin()->second.GetVal(i,j))) );
                fileA<<number;
                if( j != eigenValuesRe.begin()->second.Cols() - 1){
                    fileA<<" , ";
                }
            }
            fileA<<std::endl;
        }
        fileA.close();
    }
    TPZFMatrix<STATE> solMat(neqOriginal,1,0.);
    
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Ez");//setando para imprimir u
    std::string plotfile= "../waveguideModes.vtk";//arquivo de saida que estara na pasta debug
    int dim = 2;
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
    int postProcessResolution = 2;//define resolucao do pos processamento
    
    
    std::set<std::pair<REAL,TPZFMatrix<STATE> > > ::iterator iT = eigenValuesRe.begin();
    for (int iSol = 0; iSol < nSolutions; iSol ++) {
        if( iT == eigenValuesRe.end() ){
            DebugStop();
        }
        for (int i = 0 ; i < neq; i++) {
            solMat( activeEquations[i] ,0) = (iT->second).GetVal(i, 0);
        }
        an.LoadSolution( solMat );
        an.PostProcess(postProcessResolution);
        iT++;
    }
    std::cout << "FINISHED!" << std::endl;
}

void FilterBoundaryEquations(TPZCompMesh * cMesh , TPZVec<long> &activeEquations , int &neq , int &neqOriginal)
{
    
    TPZManVector<long,1000> allConnects;
    std::set<long> boundConnects;

    
    for (int iel = 0; iel < cMesh->NElements(); iel++) {
        TPZCompEl *cel = cMesh->ElementVec()[iel];
        if ( cel == NULL) {
            continue;
        }
        if ( cel->Reference() == NULL) {
            
            continue;
        }
        if(cel->Reference()->MaterialId() == -1 ){
            TPZMatModalAnalysisH1 *mat = dynamic_cast<TPZMatModalAnalysisH1*>(cMesh->FindMaterial(1));
            std::set<long> boundConnectsEl;
            std::set<long> depBoundConnectsEl;
            std::set<long> indepBoundConnectsEl;
            cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
            cel->BuildConnectList(boundConnectsEl);
            for (std::set<long>::iterator iT = boundConnectsEl.begin(); iT != boundConnectsEl.end(); iT++) {
                const long val = *iT;
                if( boundConnects.find(val) == boundConnects.end() ){
                    boundConnects.insert(val);
                }
            }
        }
    }
    for (int iCon = 0; iCon < cMesh->NConnects(); iCon++) {
        if( boundConnects.find( iCon ) == boundConnects.end() ){
            TPZConnect &con = cMesh->ConnectVec()[iCon];
            if( con.HasDependency() ) continue;
            int seqnum = con.SequenceNumber();
            int pos = cMesh->Block().Position(seqnum);
            int blocksize = cMesh->Block().Size(seqnum);
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
    std::cout<<"cMesh->NEquations()"<<"\t"<<cMesh->NEquations()<<std::endl;
    neqOriginal = cMesh->NEquations();
    
    int nHCurlEquations = 0 , nH1Equations = 0;
    std::cout<<"activeEquations"<<std::endl;
    long nEq = 0;
    TPZMatModalAnalysisH1 *mat = dynamic_cast<TPZMatModalAnalysisH1*>(cMesh->FindMaterial(1));
    for (int iCon = 0; iCon < cMesh->NConnects(); iCon++) {
        if( boundConnects.find( iCon ) == boundConnects.end() ){
            int seqnum = cMesh->ConnectVec()[iCon].SequenceNumber();
            int blocksize = cMesh->Block().Size(seqnum);
            for(int ieq = 0; ieq<blocksize; ieq++)
            {
                nEq++;
            }
        }
    }
    std::cout<<"# equations: "<< nEq<<std::endl;
    //std::cout<<activeEquations<<std::endl;
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

TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL f0)
{
    
    const int dim = 2; //dimensao do problema
    const int matId = 1; //define id para um material(formulacao fraca)
    const int bc0 = -1; //define id para um material(cond contorno dirichlet)
    enum{ dirichlet = 0, neumann, mixed}; //tipo da condicao de contorno do problema
    // Criando material
    
    
    ///criar malha computacional H1
    
    TPZCompMesh * cmeshH1 = new TPZCompMesh(gmesh);
    cmeshH1->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmeshH1->SetDimModel(dim);//seta dimensao do modelo
    // Inserindo material na malha
    
    TPZVec<STATE> sol;//only for creating material. this material will not be used in reality
    TPZMatModalAnalysisH1 *matH1 = new TPZMatModalAnalysisH1( matId, f0, ur, er);
    cmeshH1->InsertMaterialObject(matH1);
    
    ///electrical conductor boundary conditions
    TPZFNMatrix<1,STATE> val1(1,1,0.), val2(1,1,0.);
    
    TPZMaterial * BCondH1Dir = matH1->CreateBC(matH1, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
    
    cmeshH1->InsertMaterialObject(BCondH1Dir);//insere material na malha
    //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
    cmeshH1->SetAllCreateFunctionsContinuous();
    cmeshH1->AutoBuild();
    cmeshH1->CleanUpUnconnectedNodes();
    
    return cmeshH1;
}