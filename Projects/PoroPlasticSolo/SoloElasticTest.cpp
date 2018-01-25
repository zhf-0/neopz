/*
 *  SoloElasticTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "SoloElasticTest.h"
#include "pzcheckgeom.h"
#include "pzstack.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "AnalyticalFunctions.h"
#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"
#include "pzelasmat.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZGmshReader.h"

#define gmshmesh

using namespace std;

const REAL Pi=M_PI;

SoloElasticTest::SoloElasticTest()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    
    //Materiais das condições de contorno
    fmatBCbott=2;
    fmatBCtop=4;
    fmatBCleft=5;
    fmatBCright=3;
    
    //Material do elemento de interface
    fmatInterface=4;
    
    //Materiais das condições de contorno (elementos de interface)
    fmatIntBCbott=-11;
    fmatIntBCtop=-12;
    fmatIntBCleft=-13;
    fmatIntBCright=-14;
    
    //Materia de um ponto
    fmatPoint=-5;
    
    //Condições de contorno do problema
    fdirichlet=0;
    fneumann=1;

    
    
    fquadmat1=1; //Parte inferior do quadrado
    fquadmat2=2; //Parte superior do quadrado
    fquadmat3=3; //Material de interface
    
    fviscosity=1.;
    fpermeability=1.;
    ftheta=-1.;
    
}

SoloElasticTest::~SoloElasticTest()
{
    
}

void SoloElasticTest::Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE theta, STATE sigma)
{
    
    
    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy); //Função para criar a malha geometrica
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    
    TPZCompMesh *cmesh_v = this->CMesh_v(gmesh, Space, pOrder); //Função para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = this->CMesh_p(gmesh, Space, pOrder); //Função para criar a malha computacional da pressão
    TPZCompMesh *cmesh_m = this->CMesh_m(gmesh, Space, pOrder, visco, theta, sigma); //Função para criar a malha computacional multifísica
    
#ifdef PZDEBUG
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
        
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
#endif
    
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    //AddMultiphysicsInterfaces(*cmesh_m,fmatInterface,fmatID);
    //AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCbott,fmatBCbott);
    //AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCtop,fmatBCtop);
    //AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCleft,fmatBCleft);
    //AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCright,fmatBCright);
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
//    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    

    TPZAnalysis PoroelasticAnalysis(cmesh_m);
    
    
    //TPZParSkylineStructMatrix matskl(cmesh_m, numthreads);
    TPZSkylineNSymStructMatrix matskl(cmesh_m); //OK para Hdiv
    //TPZFStructMatrix matskl(cmesh_m); //caso nao simetrico *** //OK para discont.
    matskl.SetNumThreads(numthreads);
    PoroelasticAnalysis.SetStructuralMatrix(matskl);
    //TPZStepSolver<STATE> step;
    //step.SetDirect(ELU);
    //an.SetSolver(step);
    
    TPZStepSolver<STATE> step;		// Create Solver object
    step.SetDirect(ELDLt);			//	Symmetric case
    PoroelasticAnalysis.SetSolver(step); //	Set solver
    TPZStepSolver<STATE> & temp = dynamic_cast<TPZStepSolver<STATE> &> (PoroelasticAnalysis.Solver());
    std::list <long> & zeropivot = temp.Singular();
    if (zeropivot.size())
    {
        int eq = * zeropivot.begin();
        PoroelasticAnalysis.Rhs().Zero();
        PoroelasticAnalysis.Rhs()(eq,0) = -10000.0;
        PoroelasticAnalysis.Solve();
        TPZFMatrix<STATE> TempSolution = PoroelasticAnalysis.Solution();
        std::string output;
        output = "DumpFolder/SingularNodes";
        std::stringstream outputfiletemp;
        outputfiletemp << output << ".vtk";
        std::string plotfile = outputfiletemp.str();
        //PostProcessPoroeasticity(docHandle,ComputationalMeshVector,ComputationalMeshPoroelasticityReCurrent,PoroelasticAnalysis,plotfile,2);
        
        DebugStop();
    }
    
    
    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
    
    PoroelasticAnalysis.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global

    
    std::cout << "Solving Matrix " << std::endl;
    
    PoroelasticAnalysis.Solve();
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        PoroelasticAnalysis.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        PoroelasticAnalysis.Rhs().Print("R = ",filerhs,EMathematicaInput);
    }
#endif
    
#ifdef PZDEBUG
    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
        
        std::ofstream fileAlpha("alpha.txt");
        PoroelasticAnalysis.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    }
#endif
    
    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,6> Errors;
    ofstream ErroOut("Error_results.txt", std::ofstream::app);
    PoroelasticAnalysis.SetExact(Sol_exact);
    PoroelasticAnalysis.PostProcessError(Errors);
    
    ErroOut <<"Sigma = "<< sigma/(pOrder*pOrder*(nx-1)) << "  //  Ordem = "<< pOrder << "  //  Tamanho da malha = "<< nx-1 <<" x "<< ny-1 << std::endl;
    ErroOut <<" " << std::endl;
    //ErroOut <<"Norma H1/HDiv - V = "<< Errors[0] << std::endl;
    ErroOut <<"Norma L2 - V = "<< Errors[1] << std::endl;
    ErroOut <<"Semi-norma H1/Hdiv - V = "<< Errors[2] << std::endl;
    //ErroOut <<"Norma L2 - P = "<< Errors[4] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut.flush();
    
    //Pós-processamento (paraview):
    std::cout << "Post Processing " << std::endl;
    std::string plotfile("PoroElastSolo.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("SigmaX");
    scalnames.Push("SigmaY");
    scalnames.Push("DisplacementX");
    scalnames.Push("DisplacementY");
    scalnames.Push("FluidPressure");
    
    vecnames.Push("Displacement");
    vecnames.Push("FluxVector");
    
    
    int postProcessResolution = 3; //  keep low as possible
    
    int dim = gmesh->Dimension();
    PoroelasticAnalysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    PoroelasticAnalysis.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
}

TPZGeoMesh *SoloElasticTest::CreateGMesh(int nx, int ny, double hx, double hy)
{
    
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    
#ifdef gmshmesh
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    //Aqui é implementado um método para malhas criadas no GMSH
    

    std::string dirname = PZSOURCEDIR;
    std::string grid;
    
    grid = dirname + "/Projects/PoroPlasticSolo/gmsh_meshes/msh/Geometry1DElast.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    gmesh = Geometry.GeometricGmshMesh(grid);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    
    int n_div = 0;
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
#else
    
    int i,j;
    long id, index;
    
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    //Vetor auxiliar para armazenar coordenadas:
    
    TPZVec <REAL> coord (3,0.);
    
    
    //Inicialização dos nós:
    
    for(i = 0; i < ny; i++){
        for(j = 0; j < nx; j++){
            id = i*nx + j;
            coord[0] = (j)*hx/(nx - 1);
            coord[1] = (i)*hy/(ny - 1);
            //using the same coordinate x for z
            coord[2] = 0.;
            //cout << coord << endl;
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
        }
    }
    
    //Ponto 1
    TPZVec<long> pointtopology(1);
    pointtopology[0] = 0;
    
    gmesh->CreateGeoElement(EPoint,pointtopology,fmatPoint,id);
    
    
    //Vetor auxiliar para armazenar as conecções entre elementos:
    
    TPZVec <long> connect(4,0);
    
    
    //Conectividade dos elementos:
    
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            index = (i)*(nx - 1)+ (j);
            connect[0] = (i)*ny + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,fmatID,id);
        }
    }
    
    
    //Gerando informação da vizinhança:
    
    gmesh->BuildConnectivity();
    
    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    long el, numelements = gmesh->NElements();
    
    TPZManVector <long> TopolPlate(4);
    
    for (el=0; el<numelements; el++)
    {
        long totalnodes = gmesh->ElementVec()[el]->NNodes();
        TPZGeoEl *plate = gmesh->ElementVec()[el];
        for (int i=0; i<4; i++){
            TopolPlate[i] = plate->NodeIndex(i);
        }
        
        //Colocando as condicoes de contorno:
        TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
        TPZManVector <REAL,3> nodecoord(3);
        
        //Na face x = 1
        TPZVec<long> ncoordzbottVec(0); long sizeOfbottVec = 0;
        TPZVec<long> ncoordztopVec(0); long sizeOftopVec = 0;
        TPZVec<long> ncoordzleftVec(0); long sizeOfleftVec = 0;
        TPZVec<long> ncoordzrightVec(0); long sizeOfrightVec = 0;
        
        for (long i = 0; i < totalnodes; i++)
        {
            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (nodecoord[2] == 0. & nodecoord[1] == 0.)
            {
                sizeOfbottVec++;
                ncoordzbottVec.Resize(sizeOfbottVec);
                ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[1] == hy)
            {
                sizeOftopVec++;
                ncoordztopVec.Resize(sizeOftopVec);
                ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == 0.)
            {
                sizeOfleftVec++;
                ncoordzleftVec.Resize(sizeOfleftVec);
                ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == hx)
            {
                sizeOfrightVec++;
                ncoordzrightVec.Resize(sizeOfrightVec);
                ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
            }
        }
        
        if (sizeOfbottVec == 2) {
            int sidesbott = plate->WhichSide(ncoordzbottVec);
            TPZGeoElSide platesidebott(plate, sidesbott);
            TPZGeoElBC(platesidebott,fmatBCbott);
            TPZGeoElBC(platesidebott,fmatIntBCbott);
        }
        
        if (sizeOftopVec == 2) {
            int sidestop = plate->WhichSide(ncoordztopVec);
            TPZGeoElSide platesidetop(plate, sidestop);
            TPZGeoElBC(platesidetop,fmatBCtop);
            TPZGeoElBC(platesidetop,fmatIntBCtop);
        }
        
        if (sizeOfleftVec == 2) {
            int sidesleft = plate->WhichSide(ncoordzleftVec);
            TPZGeoElSide platesideleft(plate, sidesleft);
            TPZGeoElBC(platesideleft,fmatBCleft);
            TPZGeoElBC(platesideleft,fmatIntBCleft);
        }
        
        if (sizeOfrightVec == 2) {
            int sidesright = plate->WhichSide(ncoordzrightVec);
            TPZGeoElSide platesideright(plate, sidesright);
            TPZGeoElBC(platesideright,fmatBCright);
            TPZGeoElBC(platesideright,fmatIntBCright);
        }
        
        
        ncoordzbottVec.Resize(0);
        sizeOfbottVec = 0;
        ncoordztopVec.Resize(0);
        sizeOftopVec = 0;
        ncoordzleftVec.Resize(0);
        sizeOfleftVec = 0;
        ncoordzrightVec.Resize(0);
        sizeOfrightVec = 0;
        
    }
    
    // Criando e inserindo elemento de interfação:
    //    TPZVec<long> nodind3(2);
    //
    //    nodind3[0]=1;
    //    nodind3[1]=4;
    //
    //    gmesh->CreateGeoElement(EOned, nodind3, matInterface, index); //Criando elemento de interface (GeoElement)
    
    
    //Criando interface (Geralizado):
    
    TPZVec<long> nodint(2);
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            if(j>0&&j<(nx-1)){
                nodint[0]=j+nx*i;
                nodint[1]=j+nx*(i+1);
                gmesh->CreateGeoElement(EOned, nodint, fmatInterface, index); //Criando elemento de interface (GeoElement)
                
            }
            if(i>0&&j<(ny-1)){
                nodint[0]=j+ny*i;
                nodint[1]=j+ny*i+1;
                gmesh->CreateGeoElement(EOned, nodint, fmatInterface, index); //Criando elemento de interface (GeoElement)
                
            }
            
        }
    }
    
    
    //new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodind3,matInterface,*gmesh); //Criando elemento de interface (RefPattern)
    id++;
    
    gmesh->AddInterfaceMaterial(fquadmat1, fquadmat2, fquadmat3);
    gmesh->AddInterfaceMaterial(fquadmat2, fquadmat1, fquadmat3);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //Impressão da malha geométrica:
    
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
#endif
    

}

TPZCompEl *SoloElasticTest::CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}



//TPZGeoMesh *SoloElasticTest::GMeshDeformed(int dim, bool ftriang, int ndiv)
//{
//
//    DebugStop();
//
//}


void SoloElasticTest::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
    
    dsol.Resize(3,2);
    sol.Resize(3);
    
    REAL xv = x[0];
    REAL yv = x[1];
    
    STATE v_x =  cos(2.*Pi*yv)*sin(2.*Pi*xv);
    STATE v_y =  -(cos(2.*Pi*xv)*sin(2.*Pi*yv));
//    STATE pressure= xv*xv+yv*yv;
    STATE pressure=  -(1./Pi)*cos(Pi*xv)*exp(yv);
    
    sol[0]=v_x;
    sol[1]=v_y;
    sol[2]=pressure;
    
    // vx direction
    dsol(0,0)= 2.*Pi*cos(2.*Pi*xv)*cos(2.*Pi*yv);
    dsol(0,1)= 2.*Pi*sin(2.*Pi*xv)*sin(2.*Pi*yv);
    
    // vy direction
    dsol(1,0)= -2.*Pi*sin(2.*Pi*xv)*sin(2.*Pi*yv);
    dsol(1,1)= -2.*Pi*cos(2.*Pi*xv)*cos(2.*Pi*yv);
    
    // Gradiente pressão
    
//    dsol(2,0)= 2.*xv;
//    dsol(2,1)= 2.*yv;
    
    dsol(2,0)= exp(yv)*sin(Pi*xv);
    dsol(2,1)= -exp(yv)*cos(Pi*xv)*(1./Pi);
    
}

void SoloElasticTest::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    f.resize(2);
    
    REAL xv = x[0];
    REAL yv = x[1];
    //    STATE zv = x[2];
    
    //STATE f_x = 2.0*xv + 8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
    //STATE f_y = 2.0*yv - 8.0*Pi*Pi*cos(2.0*Pi*xv)*sin(2.0*Pi*yv);
    
    STATE f_x = exp(yv)*sin(Pi*xv) + 8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
    STATE f_y = -exp(yv)*cos(Pi*xv)*(1./Pi) - 8.0*Pi*Pi*cos(2.0*Pi*xv)*sin(2.0*Pi*yv);
    
    f[0] = f_x; // x direction
    f[1] = f_y; // y direction
    
    
}




TPZCompMesh *SoloElasticTest::CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    
    // Plane strain assumption
    int planestress = 0;
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    
    // Data for elastic problem
    REAL Eyoung			=	0.0;
    REAL PoissonRatio	=	0.0;
    REAL Lambda			=	0.0;
    REAL G				=	0.0;
    REAL RockDensity	=	0.0;
    REAL BodyForceX		=	0.0;
    REAL BodyForceY		=	0.0;
    
    TPZElasticityMaterial * MaterialElastic;

    
    //Definição do espaço de aprximação:
    
    
    //TPZMat2dLin *material = new TPZMat2dLin(fmatID); //Criando material que implementa a formulação fraca do problema modelo
    
    //cmesh->InsertMaterialObject(material); //Insere material na malha
    
    BodyForceX		= 0.;
    BodyForceY		= 0.;
    Lambda			= 4000000000.;
    G				= 6000000000.;
    RockDensity		= 2300.;
    
    Eyoung = (G*(3*Lambda+2*G))/(Lambda+G);
    PoissonRatio = (Lambda)/(2*(Lambda+G));
    MaterialElastic = new TPZElasticityMaterial(fmatID, Eyoung, PoissonRatio, BodyForceX, BodyForceY, planestress);
    TPZMaterial * Material(MaterialElastic);
    cmesh->InsertMaterialObject(Material);
    
    
    if (Space==1) {
        cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
        cmesh->ApproxSpace().CreateDisconnectedElements(true); //HDIV-Full:
        
        //Dimensões do material (para HDiv):
   //     TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
   //     material->SetMaterial(xkin, xcin, xfin);
        
    }else if(Space==2){
        
        cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
        
        //Dimensões do material (para H1 e descontinuo):
     //   TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
     //   MaterialElastic->SetMaterial(xkin, xcin, xfin);
        
        
    }else if(Space==3){
        
        cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
        //Criando elementos com graus de liberdade differentes para cada elemento (descontínuo):
        cmesh->ApproxSpace().CreateDisconnectedElements(true); //Criando elementos desconectados (descontínuo)
        
        //Dimensões do material (para H1 e descontinuo):
    //    TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //    material->SetMaterial(xkin, xcin, xfin);
        
    }
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(3,1,0.), val2(3,1,0.);
    
    TPZMaterial * BCond0 = MaterialElastic->CreateBC(Material, fmatBCbott, fneumann, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = MaterialElastic->CreateBC(Material, fmatBCtop, fneumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = MaterialElastic->CreateBC(Material, fmatBCleft, fneumann, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = MaterialElastic->CreateBC(Material, fmatBCright, fneumann, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    
    return cmesh;
   
    
}


TPZCompMesh *SoloElasticTest::CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    
    if (Space==2||Space==3) {
        pOrder--;
    }
    
    // Pure Diffusion problem
    TPZVec <REAL> convdir(3,0.);
    REAL diff	=	0.0;
    REAL conv	=	0.0;
    REAL flux	=	0.0;
    
    //Criando malha computacional:
    
    // Aproximation Space of order -> pOrder
    TPZCompEl::SetgOrder(pOrder);
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(fdim);
    cmesh->SetAllCreateFunctionsContinuous();
    
    
    //Criando material:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    //TPZMat2dLin *material = new TPZMat2dLin(fmatID);//criando material que implementa a formulacao fraca do problema modelo
    
    // Data for elastic problem
    REAL Eyoung			=	0.0;
    REAL PoissonRatio	=	0.0;
    REAL Permeability	=	0.0;
    REAL Viscosity		=	0.0;
    REAL FluidDensity	=	0.0;
    REAL BodyForceX		=	0.0;
    REAL BodyForceY		=	0.0;
    TPZMatPoisson3d * MaterialDiffusion;
    
    // Using Gravity field
    BodyForceX		= 0.;
    BodyForceY		= 0.;
    Permeability	= 0.0000000000003;
    Viscosity		= 0.001;
    FluidDensity	= 1000.;
    
    diff = (Permeability)/(Viscosity);
    
    MaterialDiffusion = new TPZMatPoisson3d(fmatID,fdim);
    MaterialDiffusion->SetParameters(diff, conv, convdir);
    MaterialDiffusion->SetInternalFlux(flux);
    MaterialDiffusion->NStateVariables();
    TPZMaterial * Material(MaterialDiffusion);
    cmesh->InsertMaterialObject(Material);
    
    //Dimensões do material (para H1 e descontínuo):
 //   TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
 //   material->SetMaterial(xkin, xcin, xfin);
    
    
    //Condições de contorno
    
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    

        TPZMaterial * BCond0 = MaterialDiffusion->CreateBC(Material, fmatBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
        cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    //
        TPZMaterial * BCond1 = MaterialDiffusion->CreateBC(Material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
        cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    //
        TPZMaterial * BCond2 = MaterialDiffusion->CreateBC(Material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
        cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //
        TPZMaterial * BCond3 = MaterialDiffusion->CreateBC(Material, fmatBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //    Ponto de pressao:
    //
  //  TPZFMatrix<STATE> val3(3,1,0.), val4(3,1,0.);
    ////
  //  TPZMaterial * BCPoint = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressao
  //  cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    //    //    Ponto de pressao2:
    //    //
    //    TPZFMatrix<STATE> val5(1,1,0.), val6(1,1,0.);
    //    ////
    //    TPZMaterial * BCPoint2 = material->CreateBC(material, matPoint2, pointtype, val5, val6); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPoint2); //Insere material na malha
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(fmatID);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *SoloElasticTest::CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE theta, STATE sigma)
{
    
    int Dimensionless = 1;
    
    // Time Data
    
    REAL Delta = 0.0001;
    REAL Theta = 0.5;
    
    // Data for poroelastic problem
    REAL LambdaU		=	0.0;
    REAL Lambda			=	0.0;
    REAL G				=	0.0;
    REAL alpha			=	0.0;
    REAL MixtureDensity	=	0.0;
    REAL RockDensity	=	0.0;
    REAL RockPorosity	=	0.0;
    REAL FluidDensity	=	0.0;
    REAL FluidViscosity	=	0.0;
    REAL Permeability	=	0.0;
    REAL BodyForceX		=	0.0;
    REAL BodyForceY		=	0.0;
    REAL diff			=	0.0;
    REAL S				=	0.0;
    REAL Se				=	0.0;
    REAL SFluid			=	0.0;
    
    // Plane strain assumption
    int planestress = 0;
    
    //Criando malha computacional:
    int bc_inte_order = 3;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    
    // Criando material:
    
    TPZPoroElastic2d * MaterialPoroElastic;
    RockDensity		= 2300.0;
    FluidDensity	= 1000.0;
    RockPorosity	= 0.3;
    Permeability	= 3.0e-13;
    FluidViscosity	= 0.001;
    Lambda			= 4.0e9;
    G				= 6.0e9;
    LambdaU			= 4.0e9;
    alpha			= 1.0;
    
    diff = (Permeability)/(FluidViscosity);
    MixtureDensity = (1 - RockPorosity) * RockDensity + RockPorosity * FluidDensity;
    SFluid = LambdaU;
    // Using Gravity field
    BodyForceX		= 0.;
    BodyForceY		= 0.;

    if (Lambda != LambdaU)
    {
        if (alpha != 0.0) {
            S = ((pow(alpha,2))/((LambdaU-Lambda)))*((LambdaU+2.0*G)/(Lambda+2.0*G));
            Se = (pow(alpha,2))/((LambdaU-Lambda));
        }
        else
        {
            S = (1)/(3*Lambda+2.0*G);
            Se = RockPorosity*(SFluid + S);
        }
        
    }
    else
    {
        
        if (alpha != 0.0) {
            S = (pow(alpha,2)/(Lambda+2.0*G));
            Se = 0;
        }
        else
        {
            S = (1.0)/(3*Lambda+2.0*G);
            Se = 0;
        }
        
    }
    
    if (Dimensionless)
    {
        if (alpha != 0.0) {
            Se				= Se/S;
        }
        else
        {
            if (Se != 0.0) {
                Se				= 1.0;
            }
        }
        Lambda			= Lambda*S;
        LambdaU			= LambdaU*S;
        G				= G*S;
        Permeability	= 1.0;
        FluidViscosity	= 1.0;
    }
    
   // TPZDummyFunction<STATE> * TimeDepFExact;
   // TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolutionSemiInfiniteColumn1D);
   // TimeDepFExact->SetPolynomialOrder(bc_inte_order);
    
    
    MaterialPoroElastic = new TPZPoroElastic2d (fmatID, fdim);
    MaterialPoroElastic->SetParameters(Lambda,G,LambdaU,BodyForceX,BodyForceY);
    MaterialPoroElastic->SetParameters(Permeability,FluidViscosity);
    MaterialPoroElastic->SetfPlaneProblem(planestress);
    MaterialPoroElastic->SetBiotParameters(alpha,Se);
  //  MaterialPoroElastic->SetTimeStep(Delta,Theta);
  //  MaterialPoroElastic->SetTimeDependentFunctionExact(TimeDepFExact);
    TPZMaterial * Material(MaterialPoroElastic);
    cmesh->InsertMaterialObject(Material);

    //TPZDummyFunction<STATE> * TimeDepFExact_fake;
    //TimeDepFExact_fake = new TPZDummyFunction<STATE> (F_source);
    //TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact);
    //material->SetForcingFunction(TimeDepFExact_fake);
    //TimeDepFExact_fake->SetPolynomialOrder(bc_inte_order);
   // material->SetForcingFunctionExact(solp);

    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    TPZMaterial * BCond0 = MaterialPoroElastic->CreateBC(Material, fmatBCbott, 101, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    val2(1,0)=-1.;
    
    TPZMaterial * BCond1 = MaterialPoroElastic->CreateBC(Material, fmatBCtop, 110, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    val2(1,0)=0.;
    
    TPZMaterial * BCond2 = MaterialPoroElastic->CreateBC(Material, fmatBCleft, 11, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = MaterialPoroElastic->CreateBC(Material, fmatBCright, 11, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    //Ponto
    
  //  TPZFMatrix<STATE> val3(3,1,0.), val4(3,1,0.);
    
 //   TPZMaterial * BCPoint = MaterialPoroElastic->CreateBC(Material, fmatPoint, fpointtype, val3, val4);//Cria material que implementa um ponto para a pressão
 //   cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    
    
    

    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }

    
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
    
}


void SoloElasticTest::AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    long nel = gmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != matfrom) {
            continue;
        }
        
        int nsides= gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        if (celstack.size() != 2) {
            DebugStop();
        }
        gel->SetMaterialId(mattarget);
        long index;
        new TPZMultiphysicsInterfaceElement(cmesh,gel,index,celstack[1],celstack[0]);
    }
    
}





