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
#include "TPZMatModalAnalysisHDiv.h"
#include "pzextractval.h"
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzlog.h"
#include "TPZTimer.h"
#include "TPZMatHCurlProjection.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzsbstrmatrix.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzl2projection.h"
#include "pzsbndmat.h"
#include "pzfstrmatrix.h"
#include "pzvisualmatrix.h"
#include <sys/types.h>
#include <sys/stat.h>


enum meshTypeE{ createRectangular=1, createTriangular, createZigZag};

STATE ur( const TPZVec<REAL> & x){
    //return ( 2.-imaginary*0.1 );
    return 1.;
}
STATE er( const TPZVec<REAL> & x){
    return 1;
}
void RunSimulation(const int meshType ,bool usingFullMtrx, bool optimizeBandwidth, int pOrder, int nDiv, REAL hDomain, REAL wDomain, REAL f0, modeType teortm, int nSolutions, bool generatingResults);

void FilterBoundaryEquations(TPZVec<TPZCompMesh *> meshVec , TPZVec<long> &activeEquations , int &neq, int &neqOriginal , const modeType teortm);

TPZVec<TPZCompMesh *>CreateCMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL f0, modeType teortm);
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
    modeType teortm = modesTE;
    REAL f0 = 25 * 1e+9;
    int nSolutions = 4;
    const int meshType = createTriangular;
    
    
    int pOrder = 1; //ordem polinomial de aproximacao
    bool usingFullMtrx = true;
    bool optimizeBandwidth = true;
    bool generatingResults = true;
    
    int nDiv = 5;
    int nSim = 5;
    
    std::string fileName;
    if(generatingResults){
        if(teortm == modesTM){
            fileName = "../resultsQuali/TM/ev";
        }
        else{
            fileName = "../resultsQuali/TE/ev";
        }
        for(int i = 1 ; i < 100 ; i++){
            fileName.append(std::to_string(i));
            fileName.append(".csv");
            if( std::ifstream( fileName.c_str() ) ){
                std::remove( fileName.c_str() );
            }
        }
    }
    
    std::cout<<"MODOS TE"<<std::endl;
    for (int i = 0 ; i < nSim; i++) {
        std::cout<<"iteration "<<i+1<<" of "<<nSim<<std::endl;
        RunSimulation( meshType, usingFullMtrx, optimizeBandwidth, pOrder, nDiv, hDomain, wDomain, f0 , teortm , nSolutions, generatingResults);
        nDiv += 5;
    }
    std::cout<<"MODOS TM"<<std::endl;
    nDiv = 5;
    teortm = modesTM;
    if(generatingResults){
        if(teortm == modesTM){
            fileName = "../resultsQuali/TM/ev";
        }
        else{
            fileName = "../resultsQuali/TE/ev";
        }
        for(int i = 1 ; i < 100 ; i++){
            fileName.append(std::to_string(i));
            fileName.append(".csv");
            if( std::ifstream( fileName.c_str() ) ){
                std::remove( fileName.c_str() );
            }
        }
    }
    
    
    for (int i = 0 ; i < nSim; i++) {
        std::cout<<"iteration "<<i+1<<" of "<<nSim<<std::endl;
        RunSimulation( meshType, usingFullMtrx, optimizeBandwidth, pOrder, nDiv, hDomain, wDomain, f0 , teortm , nSolutions, generatingResults);
        nDiv += 5;
    }
    
    
    
    
    return 0;
}

void RunSimulation( const int meshType ,bool usingFullMtrx, bool optimizeBandwidth, int pOrder, int nDiv, REAL hDomain, REAL wDomain, REAL f0, modeType teortm , int nSolutions, bool generatingResults){
    TPZTimer timer;
    
    /*delete previous results*/
    {
        if(std::ifstream("../stiffA.csv"))std::remove("../stiffA.csv");
        if(std::ifstream("../stiffB.csv"))std::remove("../stiffB.csv");
        if(std::ifstream("../EV.csv"))std::remove("../EV.csv");
        if(std::ifstream("../gmeshTriangular.txt"))std::remove("../gmeshTriangular.txt");
        if(std::ifstream("../gmeshTriangular.vtk"))std::remove("../gmeshTriangular.vtk");
        if(std::ifstream("../gmeshRectangular.txt"))std::remove("../gmeshRectangular.txt");
        if(std::ifstream("../gmeshRectangular.vtk"))std::remove("../gmeshRectangular.vtk");
        
        std::string vtkName("../waveguideModes.scal_vec.");
        for (int i = 0; i < 100 ; i++) {
            std::string testVtk = vtkName;
            testVtk.append( std::to_string(i) );
            testVtk.append( ".vtk" );
            if (std::ifstream( testVtk.c_str() ) ) {
                std::remove( testVtk.c_str() );
            }
            std::string fileName("../ev");
            fileName.append(std::to_string(nDiv));
            fileName.append(".csv");
            if( std::ifstream( fileName.c_str() ) ){
                std::remove( fileName.c_str() );
            }
        }
        
        
    }
    
    timer.start();
    
    int xDiv = nDiv;
    int zDiv = nDiv;
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    CreateGMesh(gmesh, meshType, hDomain, wDomain,  xDiv, zDiv);
    
    TPZVec<TPZCompMesh *> cMeshVec = CreateCMesh(gmesh, pOrder, ur , er , f0, teortm); //funcao para criar a malha computacional
    TPZCompMesh *cMesh = cMeshVec[0];
    
    
    TPZAnalysis an(cMesh,false);
    //configuracoes do objeto de analise
    TPZManVector<long,1000>activeEquations;
    int neq = 0;
    int neqOriginal = 0;
    
    
    
    
    TPZAutoPointer<TPZSBandStructMatrix> sbstr;
    
    TPZAutoPointer<TPZFStructMatrix> fmtrx;
    
    if(usingFullMtrx){
        fmtrx = new TPZFStructMatrix(cMesh);
        fmtrx->SetNumThreads(0);
        FilterBoundaryEquations( cMeshVec, activeEquations , neq , neqOriginal, teortm);
        fmtrx->EquationFilter().SetActiveEquations(activeEquations);

        an.SetStructuralMatrix(fmtrx);
    }
    else{
        sbstr = new TPZSBandStructMatrix(cMesh);
        sbstr->SetNumThreads(0);
        FilterBoundaryEquations( cMeshVec, activeEquations , neq , neqOriginal, teortm);
        sbstr->EquationFilter().SetActiveEquations(activeEquations);
        
        an.SetStructuralMatrix(sbstr);
    }
    nSolutions = neq >= nSolutions ? nSolutions : neq;
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky); //caso simetrico
    an.SetSolver(step);
    if (optimizeBandwidth) {
        an.OptimizeBandwidth();
    }
    int matId = 1;
    TPZMatModalAnalysisHDiv *matAlias = dynamic_cast<TPZMatModalAnalysisHDiv *> (cMesh->FindMaterial( matId ) );
    
    
    
    matAlias->SetMatrixA();
    std::cout<<"entrando no assemble matrix A"<<std::endl;
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
    
    std::string pathName;
    if(generatingResults){
        struct stat sb;
        std::string command;
        pathName = "../resultsQuali";
        if (!(stat(pathName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
        {
            command = "mkdir";
            command.append(" ../resultsQuali");
            std::system(command.c_str());
        }
        if (teortm == modesTM) {
            pathName = " ../resultsQuali/TM";
        }
        else{
            pathName = " ../resultsQuali/TE";
        }
        
        if (!(stat(pathName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))){
            command = "mkdir ";
            command.append(pathName.c_str());
            std::system(command.c_str());
        }
    }
    else{
        pathName = "..";
    }
    
    int i = 0;
    
    std::string fileName = pathName;
    fileName.append("/ev");
    fileName.append(std::to_string(nDiv));
    fileName.append(".csv");
    std::ofstream fileEigenValues(fileName.c_str());
    for (std::set<std::pair<REAL,TPZFMatrix<STATE> > > ::iterator iT = eigenValuesRe.begin(); iT != eigenValuesRe.end(); iT++) {
        fileEigenValues<< iT->first <<std::endl;
        if(iT->first < 1e-3) continue;
        i++;
        if( i >= nSolutions)
            continue;
        std::cout<< iT->first <<std::endl;
    }
    if (generatingResults) {
        std::cout << "FINISHED!" << std::endl;
        return;
    }
    std::cout << "Post Processing..." << std::endl;
    
    
    {
        std::ofstream fileA("../EV.csv");
        char number[256];
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
    if (teortm == modeType::modesTM) {
        scalnames.Push("Ez");//setando para imprimir u
    }
    else{
        scalnames.Push("Hz");//setando para imprimir u
    }
    
    vecnames.Push("Et");
    vecnames.Push("Ht");
    std::string plotfile= "../waveguideModes.vtk";//arquivo de saida que estara na pasta debug
    int dim = 2;
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
    int postProcessResolution = 2;//define resolucao do pos processamento
    
    TPZVec<TPZCompMesh *>temporalMeshVec(2);
    temporalMeshVec[ matAlias->L2Index() ] = cMeshVec[1 + matAlias->L2Index()];
    temporalMeshVec[ matAlias->HDivIndex() ] = cMeshVec[1 + matAlias->HDivIndex()];
    
    std::set<std::pair<REAL,TPZFMatrix<STATE> > > ::iterator iT = eigenValuesRe.begin();
    for (int iSol = 0; iSol < nSolutions; iSol ++) {
        if(iT->first < 1e-3) {
            iSol--;
            iT++;
            continue;
        }
        if( iT == eigenValuesRe.end() ){
            DebugStop();
        }
        for (int i = 0 ; i < neq; i++) {
            solMat( activeEquations[i] ,0) = (iT->second).GetVal(i, 0);
        }
        an.LoadSolution( solMat );
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(temporalMeshVec, cMesh);
        const REAL ktSquared = iT->first;
        matAlias->SetKtSquared(ktSquared);
        an.PostProcess(postProcessResolution);
        iT++;
    }
    std::cout << "FINISHED!" << std::endl;
}

void FilterBoundaryEquations(TPZVec <TPZCompMesh *> meshVec , TPZVec<long> &activeEquations , int &neq , int &neqOriginal , const modeType teortm)
{
    
    TPZCompMesh *cmeshMF = meshVec[0];
    TPZMatModalAnalysisHDiv *mat = dynamic_cast<TPZMatModalAnalysisHDiv*>(cmeshMF->FindMaterial(1));
    TPZCompMesh *cmeshHDiv = meshVec[1 + mat->HDivIndex()];
    TPZCompMesh *cmeshL2 = meshVec[1 + mat->L2Index()];
    
    TPZManVector<long,1000> allConnects;
    std::set<long> boundaryConnects;
    if(teortm == modeType::modesTE){
        for (int iel = 0; iel < cmeshMF->NElements(); iel++) {
            TPZCompEl *cel = cmeshMF->ElementVec()[iel];
            if ( cel == NULL) {
                continue;
            }
            if ( cel->Reference() == NULL) {
                
                continue;
            }
            if(cel->Reference()->MaterialId() == -1 ){
                std::set<long> boundConnectsEl;
                cel->BuildConnectList(boundConnectsEl);
                
                for (std::set<long>::iterator iT = boundConnectsEl.begin(); iT != boundConnectsEl.end(); iT++) {
                    const long val = *iT;
                    boundaryConnects.insert(val);
                }
            }
        }
    }
    
    for (int iCon = 0; iCon < cmeshMF->NConnects(); iCon++) {
        if( boundaryConnects.find( iCon ) == boundaryConnects.end() ){
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
    
    std::cout<<"------\t------\t-------"<<std::endl;
    std::cout<<"cmeshL2->NEquations()"<<"\t"<<cmeshL2->NEquations()<<std::endl;
    std::cout<<"cmeshHDiv->NEquations()"<<"\t"<<cmeshHDiv->NEquations()<<std::endl;
    std::cout<<"cmeshMF->NEquations()"<<"\t"<<cmeshMF->NEquations()<<std::endl;
    neqOriginal = cmeshMF->NEquations();
    
    int nHDivEquations = 0 , nL2Equations = 0;
    std::cout<<"activeEquations"<<std::endl;
    long nEq = 0;
    for (int iCon = 0; iCon < cmeshMF->NConnects(); iCon++) {
        bool isL2;
        if( boundaryConnects.find( iCon ) == boundaryConnects.end() ){
            
            int seqnum = cmeshMF->ConnectVec()[iCon].SequenceNumber();
            int blocksize = cmeshMF->Block().Size(seqnum);
            if( mat->L2Index() == 0 && iCon < cmeshL2->NConnects() ){
                isL2 = true;
            }
            else if( mat->L2Index() == 1 && iCon >= cmeshHDiv->NConnects() ){
                isL2 = true;
            }
            else{
                isL2 = false;
            }
            for(int ieq = 0; ieq<blocksize; ieq++)
            {
                nEq++;
                isL2 == true ? nL2Equations ++ : nHDivEquations++;
            }
        }
    }
    std::cout<<"# L2 equations: "<< nL2Equations<<std::endl;
    std::cout<<"# HDiv equations: "<<nHDivEquations<<std::endl;
    std::cout<<"# equations: "<<nEq<<std::endl;
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

TPZVec<TPZCompMesh *>CreateCMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL f0,modeType teortm)
{
    
    const int dim = 2; //dimensao do problema
    const int matId = 1; //define id para um material(formulacao fraca)
    const int bc0 = -1; //define id para um material(cond contorno dirichlet)
    enum{ dirichlet = 0, neumann, mixed}; //tipo da condicao de contorno do problema
    // Criando material
    
    
    ///criar malha computacional L2
    
    TPZCompMesh * cmeshL2 = new TPZCompMesh(gmesh);
    cmeshL2->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmeshL2->SetDimModel(dim);//seta dimensao do modelo
    // Inserindo material na malha
    
    TPZVec<STATE> sol;//only for creating material. this material will not be used in reality
    int nState = 1;
    TPZL2Projection *matL2 = new TPZL2Projection( matId, dim, nState, sol);
    cmeshL2->InsertMaterialObject(matL2);
    
    ///electrical conductor boundary conditions
    TPZFNMatrix<1,STATE> val1(1,1,0.), val2(1,1,0.);
    
    int bcType = teortm == modesTM ? dirichlet : neumann;
    if( teortm == modesTM){
        TPZMaterial * BCondL2Dir = matL2->CreateBC(matL2, bc0, bcType, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
        cmeshL2->InsertMaterialObject(BCondL2Dir);//insere material na malha
    }
    //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
    cmeshL2->SetAllCreateFunctionsDiscontinuous();
    cmeshL2->AutoBuild();
    cmeshL2->CleanUpUnconnectedNodes();
    
    
    ///criar malha computacional HDiv
    
    TPZCompMesh * cmeshHDiv = new TPZCompMesh(gmesh);
    cmeshHDiv->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmeshHDiv->SetDimModel(dim);//seta dimensao do modelo
    // Inserindo material na malha
    TPZVecL2 *matHDiv = new TPZVecL2(matId);
    cmeshHDiv->InsertMaterialObject(matHDiv);
    
    val1( 0, 0 ) = 0.;
    val2( 0, 0 ) = 0.;
    if(teortm == modesTE){
        TPZMaterial * BCondHDivDir = matHDiv->CreateBC(matHDiv, bc0, bcType, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
    
        cmeshHDiv->InsertMaterialObject(BCondHDivDir);//insere material na malha
    }
    
    cmeshHDiv->SetAllCreateFunctionsHDiv();//define espaco de aproximacao
    cmeshHDiv->AutoBuild();
    cmeshHDiv->CleanUpUnconnectedNodes();
    
    
    TPZMatModalAnalysisHDiv *matMultiPhysics = new TPZMatModalAnalysisHDiv(matId , f0 , ur , er);
    TPZVec<TPZCompMesh *> meshVec(2);
    
    meshVec[ matMultiPhysics->L2Index() ] = cmeshL2;
    meshVec[ matMultiPhysics->HDivIndex() ] = cmeshHDiv;
    
    //  TPZMatHCurl2D *matMultiPhysics = new TPZMatHCurl2D(matId , lambda , ur , er , e0 , theta, scale);
    
    
    TPZCompMesh *cmeshMF = new TPZCompMesh( gmesh );
    
    val1( 0, 0 ) = 0.;
    val2( 0, 0 ) = 0.;
    TPZMaterial * BCondMFDir = matMultiPhysics->CreateBC(matMultiPhysics, bc0, bcType, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
    
    
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
    
    TPZVec< TPZCompMesh *>meshVecOut(3);
    
    meshVecOut[0] = cmeshMF;
    meshVecOut[1 + matMultiPhysics->L2Index()] = cmeshL2;
    meshVecOut[1 + matMultiPhysics->HDivIndex()] = cmeshHDiv;
    return meshVecOut;    
}