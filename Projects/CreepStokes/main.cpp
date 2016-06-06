

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZStokesMaterial.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiPhysicsInterfaceEl.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"

//------------------STOKES Creep Of Concrete------------------------




/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidimensional formada por nel elementos de tamanho elsize
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy);

/**
 * @brief Funcao para criar a malha computacional da velocidade a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da pressão a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional multi-fisica ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);


//Função para criar interface entre elmentos:

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
const int quadmat1=  1; // Parte inferior do quadrado
const int quadmat2 =  2; // Parte superior do quadrado
const int matInterface = 4;
const int quadmat3=3;// Material de interface

void AddMultiphysicsInterfaces(TPZCompMesh &cmesh);

using namespace std;

// definition of f
void f_source(const TPZVec<REAL> & x, TPZVec<STATE>& f);

// definition of v analytic
void v_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f);

//Função principal do programa:

int main(int argc, char *argv[])
{

    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    double hx=1.,hy=0.5; //Dimensões em x e y do domínio
    int nelx=2, nely=1; //Número de elementos em x e y
    int nx=nelx+1 ,ny=nely+1; //Número de nos em x  y
    int pOrder = 1; //Ordem polinomial de aproximação
    //int dim = 2; //Dimensão do problema
    //double elsizex=hx/nelx, elsizey=hy/nely; //Tamanho dos elementos
    //int nel = elsizex*elsizey; //Número de elementos a serem utilizados
    
    
    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy); //Função para criar a malha geometrica
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    
    TPZCompMesh *cmesh_v = CMesh_v(gmesh, pOrder); //Função para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = CMesh_p(gmesh, pOrder); //Função para criar a malha computacional da pressão
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); //Função para criar a malha computacional multifísica

#ifdef PZDEBUG
    std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
    std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
    cmesh_v->Print(filecv);
    cmesh_p->Print(filecp);
#endif
    
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    AddMultiphysicsInterfaces(*cmesh_m);

#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif
    
    //Resolvendo o Sistema:
    
    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineNSymStructMatrix matskl(cmesh_m); //caso simetrico
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equações

#ifdef PZDEBUG
    //Imprimindo vetor solução:

    TPZFMatrix<REAL> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
    
    solucao.Print("Sol",cout,EMathematicaInput);//Imprime na formatação do Mathematica
#endif
    

    //Imprimir Matriz de rigidez Global:
    
    std::ofstream filestiff("stiffness.txt");
    an.Solver().Matrix()->Print("K = ",filestiff,EMathematicaInput);

    std::ofstream filerhs("rhs.txt");
    an.Rhs().Print("R = ",filerhs,EMathematicaInput);
    
    std::ofstream fileAlpha("alpha.txt");
    an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    
    an.Solve();
    
    //Pós-processamento (paraview):

    std::string plotfile("Stokes.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Velocity");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    
    
    int postProcessResolution = 4; //  keep low as possible
    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}

// definition of f
void f_source(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    f.resize(2);
    
    STATE xv = x[0];
    STATE yv = x[1];
//    STATE zv = x[2];
    
    STATE f_x = 4.*((xv*xv*xv)*(6. - 12.*yv) + (xv*xv*xv*xv)*(-3. + 6.*yv) +
                    yv*(1. - 3.*yv + 2.*(yv*yv)) - 6.*xv*yv*(1. - 3.*yv + 2.*(yv*yv)) +
                    3.*(xv*xv)*(-1. + 4.*yv - 6.*(yv*yv) + 4.*(yv*yv*yv)));
    STATE f_y = -4.*(-3.*((-1. + yv)*(-1. + yv))*(yv*yv) - 3.*(xv*xv)*(1. - 6.*yv + 6.*(yv*yv)) +
                     2.*(xv*xv*xv)*(1. - 6.*yv + 6.*(yv*yv)) +
                     xv*(1. - 6.*yv + 12.*(yv*yv) - 12.*(yv*yv*yv) + 6.*(yv*yv*yv*yv)));
    
    f[0] = f_x; // x direction
    f[1] = f_y; // y direction
}

// definition of v analytic
void v_exact(const TPZVec<REAL> & x, TPZVec<STATE>& f){
    f.resize(2);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE v_x =  -2.*((-1. + xv)*(-1. + xv))*(xv*xv)*(-1. + yv)*yv*(-1. + 2.*yv);
    STATE v_y =  +2.*(-1. + xv)*xv*(-1. + 2.*xv)*((-1. + yv)*(-1. + yv))*(yv*yv);
    
    f[0] = v_x; // x direction
    f[1] = v_y; // y direction
    
}

TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy)
{
 
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
            gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
        }
    }
    
    
    //Gerando informação da vizinhança:
    
    gmesh->BuildConnectivity();
    
    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    long el, numelements = gmesh->NElements();
    int  dirbottID = -1, dirtopID = -2, dirleftID = -3,dirrightID = -4;
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
            TPZGeoElBC(platesidebott,dirbottID);
        }
        
        if (sizeOftopVec == 2) {
            int sidestop = plate->WhichSide(ncoordztopVec);
            TPZGeoElSide platesidetop(plate, sidestop);
            TPZGeoElBC(platesidetop,dirtopID);
        }
        if (sizeOfleftVec == 2) {
            int sidesleft = plate->WhichSide(ncoordzleftVec);
            TPZGeoElSide platesideleft(plate, sidesleft);
            TPZGeoElBC(platesideleft,dirleftID);
        }
        
        if (sizeOfrightVec == 2) {
            int sidesright = plate->WhichSide(ncoordzrightVec);
            TPZGeoElSide platesideright(plate, sidesright);
            TPZGeoElBC(platesideright,dirrightID);
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
    TPZVec<long> nodind3(2);
    
    nodind3[0]=1;
    nodind3[1]=4;
    
    //gmesh->CreateGeoElement(EOned, nodind3, matInterface, index); //Criando elemento de interface (GeoElement)
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodind3,matInterface,*gmesh); //Criando elemento de interface (RefPattern)
    id++;
    
    gmesh->AddInterfaceMaterial(quadmat1, quadmat2, quadmat3);
    gmesh->AddInterfaceMaterial(quadmat2, quadmat1, quadmat3);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //Impressão da malha geométrica:
    
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;

    
    
}

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}


TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder)
{
    const int dim = 2; //Dimensão do problema
    const int matId = 1, bc0 = -1, bc1 = -2, bc2=-3, bc3=-4; //Mesmos ids da malha geometrica
    const int dirichlet = 0, neumann = 1, mixed = 2; //Definindo condições de contorno do problema ->default dirichlet na esquerda e na direita
    REAL visco=1.; //Coeficiente de viscosidade
    
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim);//Insere dimensão do modelo
  
    
    //Definição do espaço de aprximação:
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
    
    //cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
    
    
    //Criando elementos com graus de liberdade differentes para cada elemento (descontínuo):
    
    //cmesh->ApproxSpace().CreateDisconnectedElements(true); //Criando elementos desconectados (descontínuo)
    
    
    //Criando material:
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *material = new TPZMat2dLin(matId); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    material->SetMaterial(xkin, xcin, xfin);
    
    //Dimensões do material (para HDiv):
    //TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    //material->SetMaterial(xkin, xcin, xfin);

    
    //Condições de contorno:

    TPZFMatrix<REAL> val1(1,1,0.), val2(1,1,0.);
    
    TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2); //Cria material que implementa a condição de contorno da esquerda
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno da direita
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno inferior
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
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

TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{
    const int dim = 2; //Dimensão do problema
    const int matId = 1, bc0 = -1, bc1 = -2, bc2=-3, bc3=-4; //Mesmos ids da malha geométrica
    const int dirichlet = 0, neumann = 1, mixed = 2; //Definindo condição de contorno do problema ->default Dirichlet na esquerda e na direita
    REAL visco=1.; //Coeficiente de viscosidade
    
    // @omar::
    
    //pOrder--; // Space restriction apapapa
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    // @omar::
   // cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    //Criando material:
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *material = new TPZMat2dLin(matId);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    
    
//    ///Inserir condicao de contorno esquerda
//    TPZFMatrix<REAL> val1(1,1,0.), val2(1,1,0.);
//    TPZMaterial * BCond0 = material->CreateBC(material, bc0, neumann, val1, val2);//cria material que implementa a condicao de contorno da esquerda
//    
//    cmesh->InsertMaterialObject(BCond0);//insere material na malha
//    
//    // Condicao de contorno da direita
//    TPZMaterial * BCond1 = material->CreateBC(material, bc1, neumann, val1, val2);//cria material que implementa a condicao de contorno da direita
//    
//    cmesh->InsertMaterialObject(BCond1);//insere material na malha
//    
//    val2(0,0) = 1.0;//potencial na placa inferior
//    // Condicao de contorno da placa inferior
//    TPZMaterial * BCond2 = material->CreateBC(material, bc2, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da placa inferior
//    
//    cmesh->InsertMaterialObject(BCond2);//insere material na malha
//    
//    val2(0,0) = 1.5;//potencial na placa superior
//    // Condicao de contorno da placa superior
//    TPZMaterial * BCond3 = material->CreateBC(material, bc3, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da placa superior
//    
//    cmesh->InsertMaterialObject(BCond3);//insere material na malha
    

    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha

    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }

    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AutoBuild();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder)
{
    const int dim = 2; //Dimensão do problema
    const int matId = 1, bc0 = -1, bc1 = -2, bc2=-3, bc3=-4; //Mesmos ids da malha geometrica
    const int dirichlet = 0, neumann = 1, mixed = 2; //Definindo condições de contorno do problema ->default Dirichlet na esquerda e na direita
    REAL visco=1.,theta=-1.;
    
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    
    // Criando material:
    
    TPZStokesMaterial *material = new TPZStokesMaterial(matId,dim,visco,theta);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (f_source);
    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (v_exact);
    material->SetForcingFunction(fp);
    material->SetForcingFunctionExact(vp);
    cmesh->InsertMaterialObject(material);
    
    
    //Condições de contorno:
    
    TPZFMatrix<REAL> val1(1,1,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
    TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2); //Cria material que implementa a condição de contorno da esquerda
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno da direita
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno inferior
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    
#ifdef PZDEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

void AddMultiphysicsInterfaces(TPZCompMesh &cmesh)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    long nel = gmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != matInterface) {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        if (celstack.size() != 2) {
            DebugStop();
        }
        gel->SetMaterialId(1);
        long index;
        new TPZMultiphysicsInterfaceElement(cmesh,gel,index,celstack[1],celstack[0]);
    }
}
