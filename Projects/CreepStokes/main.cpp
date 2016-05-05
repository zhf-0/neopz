

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

///Adiciona namespace std
using namespace std;

///Funcao principal do programa
int main(int argc, char *argv[])
{

//    int dim = 2; //dimensao do problema
    double hx=2.,hy=2.; //dimensoes em x e y do dominio
    
    int nelx=2, nely=1; //nuemero de elementos em x e y
    int nx=nelx+1 ,ny=nely+1; //numero de nos em x  y
    //double elsizex=hx/nelx, elsizey=hy/nely; //tamanho dos elementos
    //int nel = elsizex*elsizey; //numero de elementos a serem utilizados
    int pOrder = 1; //ordem polinomial de aproximacao
    
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy); //funcao para criar a malha geometrica
    
    std::ofstream file("MalhaGeo.txt");
    std::ofstream filevtk("MalhaGeo.vtk");
    
    gmesh->Print(file);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filevtk,true);

    TPZCompMesh *cmesh_v = CMesh_v(gmesh, pOrder); //funcao para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = CMesh_p(gmesh, pOrder); //funcao para criar a malha computacional da pressão
    
    std::ofstream filecv("MalhaC_v.txt");
    std::ofstream filecp("MalhaC_p.txt");
    cmesh_v->Print(filecv);
    cmesh_p->Print(filecp);
    
    
    
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); //funcao para criar a malha computacional da pressão
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    
    std::ofstream filecm("MalhaC_m.txt");
    cmesh_m->Print(filecm);
    
    // Resolvendo o Sistema
    bool optimizeBandwidth = false; //impede a renumeracao das equacoes do problema(para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //cria objeto de analise que gerenciaria a analise do problema
    an.Run();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
    TPZFMatrix<REAL> solucao=cmesh_m->Solution();//Pegando o vetor de solucao, alphaj
    solucao.Print("Sol",cout,EMathematicaInput);//imprime na formatacao do Mathematica
    
    //fazendo pos processamento para paraview

    
    std::string plotfile("placaaf.vtk");
//    TPZVec <std::string> scalnames(3), vecnames(1);
//    vecnames[0] = "Displacement";
//    scalnames[0] = "Mn1";
//    scalnames[1] = "Mn2";
//    scalnames[2] = "Mn1n2";
    
    
    int postProcessResolution = 3;//define resolucao do pos processamento
    an.PostProcess(postProcessResolution);//realiza pos processamento
    
    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}



TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy)
{
 
    int i,j;
    long id, index;
    
    //Creates the geometric mesh... The nodes and elements
    //will be inserted into mesh object during initilize process
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    //Auxiliar vector to store a coordinate
    TPZVec <REAL> coord (3,0.);
    
    //Nodes initialization
    for(i = 0; i < nx; i++){
        for(j = 0; j < ny; j++){
            id = i*ny + j;
            coord[0] = (i)*hx/(nx - 1);
            coord[1] = (j)*hy/(ny - 1);
            //using the same coordinate x for z
            coord[2] = 0.;
            //cout << coord << endl;
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
        }
    }
    
    //Auxiliar vector to store a element connectivities
    TPZVec <long> connect(4,0);
    
    //Element connectivities
    for(i = 0; i < (nx - 1); i++){
        for(j = 0; j < (ny - 1); j++){
            index = (i)*(ny - 1)+ (j);
            connect[0] = (i)*ny + (j);
            connect[1] = connect[0]+(ny);
            connect[2] = connect[1]+1;
            connect[3] = connect[0]+1;
            gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
        }
    }
    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
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
        
        // Colocando as condicoes de contorno
        TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
        TPZManVector <REAL,3> nodecoord(3);
        // na face x = 1
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
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;

    
    
}

TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder)
{
    const int dim = 2; //dimensao do problema
    const int matId = 1, bc0 = -1, bc1 = -2, bc2=-3, bc3=-4; //MESMOS ids da malha geometrica
    const int dirichlet = 0, neumann = 1, mixed = 2; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita
    REAL visco=1.;
    
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    cmesh->SetAllCreateFunctionsDiscontinuous(); // Setting up h1 approximation space
    
    // Criando material
    TPZStokesMaterial *material = new TPZStokesMaterial(matId,dim,visco);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
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
    
    //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{
    const int dim = 2; //dimensao do problema
    const int matId = 1, bc0 = -1, bc1 = -2, bc2=-3, bc3=-4; //MESMOS ids da malha geometrica
    const int dirichlet = 0, neumann = 1, mixed = 2; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita
    REAL visco=1.;
    
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    cmesh->SetAllCreateFunctionsDiscontinuous(); // Setting up h1 approximation space
    
    // Criando material
    TPZStokesMaterial *material = new TPZStokesMaterial(matId,dim,visco);//criando material que implementa a formulacao fraca do problema modelo
    
        // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
//    
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
    
    //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
    cmesh->AutoBuild();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder)
{
    const int dim = 2; //dimensao do problema
    const int matId = 1, bc0 = -1, bc1 = -2, bc2=-3, bc3=-4; //MESMOS ids da malha geometrica
    const int dirichlet = 0, neumann = 1, mixed = 2; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita
    REAL visco=1.;
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    // Criando material
    TPZStokesMaterial *material = new TPZStokesMaterial(matId,dim,visco);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
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
////    
//    cmesh->InsertMaterialObject(BCond3);//insere material na malha
    
    //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}