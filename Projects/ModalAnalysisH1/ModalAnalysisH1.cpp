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
#include <sstream>
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
#include <sys/types.h>
#include <sys/stat.h>

enum meshTypeE{ createRectangular=1, createTriangular, createZigZag};

STATE ur( const TPZVec<REAL> & x){
    //return ( 2.-imaginary*0.1 );
    return 1.;
}
STATE er( const TPZVec<REAL> & x){
    return 1.;
}
void RunSimulation(const int meshType ,bool usingFullMtrx, bool optimizeBandwidth, int pOrder, int nDiv, REAL hDomain, REAL wDomain, REAL f0, modeType teortm, bool generatingResults);

void CheckEigenValues(TPZMatrix<STATE> *stiffA , TPZMatrix<STATE> *stiffB , TPZVec<std::pair<STATE,TPZFMatrix<STATE> > > &eigenValuesComplex , bool usingFullMtrx);

void FilterBoundaryEquations(TPZCompMesh * cMesh , TPZVec<long> &activeEquations , int &neq, int &neqOriginal);

TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL f0, modeType teortm);

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain,  const int xDiv, const int zDiv);

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    //PARAMETROS FISICOS DO PROBLEMA
    REAL hDomain = 4 * 2.54 * 1e-3;
    REAL wDomain = 9 * 2.54 * 1e-3;
    const modeType teortm = modesTE;
    REAL f0 = 25 * 1e+9;
    int pOrder = 1; //ordem polinomial de aproximacao

    bool usingFullMtrx = false;
    bool optimizeBandwidth = false;
    bool generatingResults = true;
    const int meshType = createTriangular;
    
    int nDiv = 5;
    int nSim = 5;
    
    
    
    if(generatingResults){
        std::string fileName;
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
        RunSimulation(meshType, usingFullMtrx, optimizeBandwidth, pOrder, nDiv, hDomain, wDomain, f0 , teortm, generatingResults);
        nDiv += 5;
    }
    
    
    
    return 0;
}

void RunSimulation(const int meshType ,bool usingFullMtrx, bool optimizeBandwidth, int pOrder, int nDiv, REAL hDomain, REAL wDomain, REAL f0, modeType teortm, bool generatingResults){
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
        for (int i = 0; ; i++) {
            std::string testVtk = vtkName;
            testVtk.append( std::to_string(i) );
            testVtk.append( ".vtk" );
            if (std::ifstream( testVtk.c_str() ) ) {
                std::remove( testVtk.c_str() );
            }else{
                break;
            }
        }
        std::string fileName("../ev");
        fileName.append(std::to_string(nDiv));
        fileName.append(".txt");
        if( std::ifstream( fileName.c_str() ) ){
            std::remove( fileName.c_str() );
        }
        
    }

    
    timer.start();
    
    int xDiv = nDiv;
    int zDiv = nDiv;
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    CreateGMesh(gmesh, meshType, hDomain, wDomain,  xDiv, zDiv);
    
    TPZCompMesh *cMesh = CreateCMesh(gmesh, pOrder, ur , er , f0, teortm); //funcao para criar a malha computacional
    
    
    
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
        if (teortm == modeType::modesTM) {
            FilterBoundaryEquations( cMesh, activeEquations , neq , neqOriginal);
            fmtrx->EquationFilter().SetActiveEquations(activeEquations);
            if( nDiv < 8) {
                std::cout<<"activeEquations"<<std::endl;
                for (int i = 0; i< activeEquations.size(); i++) {
                    std::cout<<activeEquations[i]<<std::endl;
                }
            }
        }
        else{
            neq = neqOriginal = cMesh->NEquations();
            std::cout<<"nEq: "<<neq<<std::endl;
        }
        an.SetStructuralMatrix(fmtrx);
    }
    else{
        sbstr = new TPZSBandStructMatrix(cMesh);
        sbstr->SetNumThreads(0);
        if (teortm == modeType::modesTM) {
            FilterBoundaryEquations( cMesh, activeEquations , neq , neqOriginal);
            sbstr->EquationFilter().SetActiveEquations(activeEquations);
            if( nDiv < 8) {
                std::cout<<"activeEquations"<<std::endl;
                for (int i = 0; i< activeEquations.size(); i++) {
                    std::cout<<activeEquations[i]<<std::endl;
                }
            }
        }
        else{
            neq = neqOriginal = cMesh->NEquations();
            std::cout<<"nEq: "<<neq<<std::endl;
        }
        an.SetStructuralMatrix(sbstr);
    }
    int nSolutions = neq >= 10 ? 10 : neq;
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky); //caso simetrico
    an.SetSolver(step);
    if (optimizeBandwidth) {
        an.OptimizeBandwidth();
    }
    int matId = 1;
    TPZMatModalAnalysisH1 *matAlias = dynamic_cast<TPZMatModalAnalysisH1 *> (cMesh->FindMaterial( matId ) );
    
    
    
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
    
    
    TPZVec<std::pair<STATE,TPZFMatrix<STATE> > > eigenValuesComplex( stiffAPtr->Rows() );
    std::pair<STATE,TPZFMatrix<STATE> > dupletComplex;
    
    for(int i = 0 ; i< eValues.size();i++){
        eVectors.GetSub(0, i, eVectors.Rows(), 1 , eVector);
        duplet.first = eValues[i].real();
        duplet.second = eVector;
        eigenValuesRe.insert(duplet);
        
        dupletComplex.first = eValues[i];
        dupletComplex.second = eVector;
        eigenValuesComplex[i] = dupletComplex;
        
    }
    int i = 0;
//    std::string command = "mkdir";
//    command.append(" ../test");
//    std::system(command.c_str());

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
//    std::cout << "Checking eigenvalues" << std::endl;
//    CheckEigenValues(stiffAPtr , stiffBPtr , eigenValuesComplex, usingFullMtrx);
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
    int postProcessResolution = 1;//define resolucao do pos processamento
    
    
    std::set<std::pair<REAL,TPZFMatrix<STATE> > > ::iterator iT = eigenValuesRe.begin();
    for (int iSol = 0; iSol < nSolutions; iSol ++) {
        if( nDiv > 10)
            if(iT->first < 1e-3) {
                iSol--;
                iT++;
                continue;
            }
        if( iT == eigenValuesRe.end() ){
            DebugStop();
        }
        for (int i = 0 ; i < neq; i++) {
            if (teortm == modeType::modesTM) {
                solMat( activeEquations[i] ,0) = (iT->second).GetVal(i, 0);
            }
            else{
                solMat( i ,0) = (iT->second).GetVal(i, 0);
            }
            
        }
        an.LoadSolution( solMat );
        const REAL ktSquared = iT->first;
        matAlias->SetKtSquared(ktSquared);
        an.PostProcess(postProcessResolution);
        iT++;
    }
    std::cout << "FINISHED!" << std::endl;
}
void CheckEigenValues(TPZMatrix<STATE> *stiffA , TPZMatrix<STATE> *stiffB , TPZVec<std::pair<STATE,TPZFMatrix<STATE> > > &eigenValuesComplex , bool usingFullMtrx)
{

    
    const int dim = stiffA->Rows();
    
    int nWrong = 0;
    TPZManVector<STATE> res(dim,0.);
    
    for(int i = 0 ; i < 10 ; i++){
        for(int j = 0; j < 10; j++){
            std::cout<<std::real(stiffA->GetVal(i, j))<<"\t";
        }
        std::cout<<std::endl;
    }
    
    
    for(int i = 0 ; i < 10 ; i++){
        for(int j = 0; j < 10; j++){
            std::cout<<std::real(stiffB->GetVal(i, j))<<"\t";
        }
        std::cout<<std::endl;
    }
    
    for (int iEigen = 0; iEigen < dim; iEigen++) {
        STATE lambda = eigenValuesComplex[iEigen].first;
        TPZFMatrix<STATE> v = eigenValuesComplex[iEigen].second;
        if(std::norm(lambda) < 1e-3) continue;
        STATE maxVal = 0.;
        for(int i = 0 ; i < dim ; i++){
            STATE ax = 0.;
            STATE bx = 0.;
            for(int j = 0; j < dim; j++){
                ax += stiffA->GetVal(i, j)*v(j,0);
                bx += stiffB->GetVal(i, j)*v(j,0);
            }
            res[i] = ax - lambda * bx;
            if ( std::norm(res[i]) > std::norm(maxVal) ) {
                maxVal = res[i];
            }
        }
        if( std::norm(maxVal) > 1e-1){
            nWrong++;
        }
    }
    std::cout<<"# of wrong (lambda, v) pairs: "<< nWrong<<std::endl;
}

void FilterBoundaryEquations(TPZCompMesh * cMesh , TPZVec<long> &activeEquations , int &neq , int &neqOriginal)
{
    
    TPZManVector<long,1000> allConnects;
    std::set<long> boundaryConnects;

    
    for (int iel = 0; iel < cMesh->NElements(); iel++) {
        TPZCompEl *cel = cMesh->ElementVec()[iel];
        if ( cel == NULL) {
            continue;
        }
        if ( cel->Reference() == NULL) {
            DebugStop();
        }
        if(cel->Reference()->MaterialId() == -1 ){
            std::set<long> connectsEl;
            cel->BuildConnectList(connectsEl);
            for (std::set<long>::iterator iT = connectsEl.begin(); iT != connectsEl.end(); iT++) {
                const long val = *iT;
                if( boundaryConnects.find(val) == boundaryConnects.end() ){
                    boundaryConnects.insert(val);
                }
            }
        }
    }
    for (int iCon = 0; iCon < cMesh->NConnects(); iCon++) {
        if( boundaryConnects.find( iCon ) == boundaryConnects.end() ){
            TPZConnect &con = cMesh->ConnectVec()[iCon];
            
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
    
    std::cout<<"activeEquations"<<std::endl;
    long nEq = 0;
    for (int iCon = 0; iCon < cMesh->NConnects(); iCon++) {
        if( boundaryConnects.find( iCon ) == boundaryConnects.end() ){
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
    llCoord[0] = 0;
    llCoord[1] = 0;
    
    ulCoord[0] = 0;
    ulCoord[1] = hDomain;
    
    urCoord[0] = wDomain;
    urCoord[1] = hDomain;
    
    lrCoord[0] = wDomain;
    lrCoord[1] = 0;
    
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

TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( const TPZVec<REAL> &),STATE (& er)( const TPZVec<REAL> &), REAL f0,modeType teortm)
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
    
    int bcType = teortm == modesTM ? dirichlet : neumann;
    
    TPZMaterial * BCondH1Dir = matH1->CreateBC(matH1, bc0, bcType, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
    
    cmeshH1->InsertMaterialObject(BCondH1Dir);//insere material na malha
    //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
    cmeshH1->SetAllCreateFunctionsContinuous();
    cmeshH1->AutoBuild();
    cmeshH1->CleanUpUnconnectedNodes();
    
    return cmeshH1;
}