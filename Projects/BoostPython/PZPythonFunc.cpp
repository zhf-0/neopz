#include "PZPythonFunc.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "TPZMatModelProblem.h"
#include "pzbndcond.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

using namespace std;

boost::python::object CreateGMesh(long nel, REAL elsize)
{
  // Inicializa objeto da classe TPZGeoMesh
  TPZGeoMesh * gmesh = new TPZGeoMesh;

  long nnodes = nel + 1; // Numero de nos do problema
  gmesh->NodeVec().Resize(nnodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
  int mat1d = 1; // Define id para um material(formulacao fraca)
  int bc0 = -1; // Define id para um material(cond contorno esq)
  int bc1 = -2; // Define id para um material(cond contorno dir)

  // Colocando nos na malha
  for (long i = 0 ; i < nnodes; i++)
  {
    const REAL pos = i * elsize;
    TPZVec <REAL> coord(3,0.);
    coord[0] = pos;
    gmesh->NodeVec()[i].SetCoord(coord); // Seta coordenada de um no no vetor de nos da malha
    gmesh->NodeVec()[i].SetNodeId(i); // Atribui identificacao para um no
  }

  // Criando Elementos
  TPZVec <long> topol(2); // Vetor que sera inicializado com o indice dos nos de um elemento unidimensional
  TPZVec <long> TopolPoint(1); // Vetor que sera inicializado com o indice do no de um elemento zero-dimensional
  long id; // Id do elemento que sera preenchido pelo metodo CreateGeoElement

  for (long iel = 0; iel < nel; iel++)
  {
    const long ino1 = iel;
    const long ino2 = iel + 1;
    topol[0] = ino1;
    topol[1] = ino2;
    gmesh->CreateGeoElement(EOned, topol, mat1d, id);// Cria elemento unidimensional
    gmesh->ElementVec()[id];
  }

  // Cond Contorno esquerda
  TopolPoint[0] = 0;
  gmesh->CreateGeoElement(EPoint, TopolPoint, bc0, id);

  // Cond Contorno Direita
  TopolPoint[0] = nnodes-1;
  gmesh->CreateGeoElement(EPoint, TopolPoint, bc1, id);

  gmesh->BuildConnectivity(); // Constroi a conectividade de vizinhanca da malha

  return boost::python::object(gmesh);
}

boost::python::object CMesh(TPZGeoMesh *gmesh, int pOrder)
{
  const int dim = 1; // Dimensao do problema
  const int matId = 1, bc0 = -1, bc1 = -2; // MESMOS ids da malha geometrica
  const int dirichlet = 0, neumann = 1, mixed = 2; // Tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita

  // Criando material
  TPZMatModelProblem *material = new TPZMatModelProblem(matId); // Criando material que implementa a formulacao fraca do problema modelo
  TPZMatModelProblem *material2 = new TPZMatModelProblem(2);// Criando material que implementa a formulacao fraca do problema modelo

  // Criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmesh->SetDimModel(dim);//seta dimensao do modelo

  // Inserir condicao de contorno esquerda
  TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
  TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);// Cria material que implementa a condicao de contorno da esquerda

  // Condicao de contorno da direita
  TPZMaterial * BCond1 = material->CreateBC(material, bc1, dirichlet, val1, val2); //cria material que implementa a condicao de contorno da direita

  cmesh->InsertMaterialObject(BCond0); //insere material na malha
  cmesh->InsertMaterialObject(BCond1); //insere material na malha
  // Inserindo material na malha
  cmesh->InsertMaterialObject(material);
  cmesh->InsertMaterialObject(material2);


  //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
  cmesh->AutoBuild();

  material->Print();

  TPZMaterial * result = 0;
  std::map<int, TPZMaterial * >::iterator mit;
  mit = cmesh->MaterialVec().find(1);
  if(mit != cmesh->MaterialVec().end())
  {
    result = mit->second;
    std::cout << "foundFunc!" << endl;
  }
  if(result)
  {
    result->Print();
  }
  else
  {
    std::cout << "null ptr Func" << endl;
    DebugStop();
  }

  return boost::python::object(cmesh);
}
