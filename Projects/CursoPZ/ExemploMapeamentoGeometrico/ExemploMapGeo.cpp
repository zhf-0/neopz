#include <iostream>


#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "tpzquadrilateral.h"
#include "TPZQuadSphere.h"
#include "tpzgeoblend.h"
#include "TPZWavyLine.h"
#include "tpztriangle.h"
#include "TPZVTKGeoMesh.h"
#include "tpzgeoelmapped.h"
#include <tpzarc3d.h>

void QuadWavy();
void QuadWavyMapped();
void TriCircleMapped();
void TriWavyMapped();
void CubeWavyMapped();
void QuadSphere();
void ExemploIni();

int main(int argc, char *argv[])
{
    TriCircleMapped();
	return 0;
}

void QuadWavy()
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nnodes = 4;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node0, node1, node2, node3;
    
    TPZVec<REAL> coord(3,0.);
    
    //no 0
    coord[0] = 0.;
    coord[1] = 0.;
    coord[2] = 0.;
    node0.SetNodeId(0);
    node0.SetCoord(coord);
    geomesh->NodeVec()[0] = node0;
    
    //no 1
    coord[0] = 1.;
    coord[1] = 0.;
    coord[2] = 0.;
    node1.SetNodeId(1);
    node1.SetCoord(coord);
    geomesh->NodeVec()[1] = node1;
    
    //no 2
    coord[0] = 1.;
    coord[1] = 2.;
    coord[2] = 0.;
    node2.SetNodeId(2);
    node2.SetCoord(coord);
    geomesh->NodeVec()[2] = node2;
    
    //no 3
    coord[0] = 0.;
    coord[1] = 1.;
    coord[2] = 0.;
    node3.SetNodeId(3);
    node3.SetCoord(coord);
    geomesh->NodeVec()[3] = node3;
    
    //INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
    TPZVec<long> topology(4);//Quadrilatero ira utilizar 4 nos
    topology[0] = 0;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 2;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
    
    int materialId = 1;
    
    TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > * quadrilatero =
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (topology,materialId,*geomesh);
    
    
    //INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
    TPZVec<long> topologyWavy(2);//Quadrilatero ira utilizar 4 nos
    topologyWavy[0] = 1;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topologyWavy[1] = 2;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    
    int materialIdW = -1;
    
    TPZGeoElRefPattern< pzgeom::TPZWavyLine > * wavy =
    new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy,materialIdW,*geomesh);
    
    TPZVec<REAL> wavedir(3,0.);
    wavedir[0] = 0.;
    wavedir[2] = 0.1;
    int numwaves = 5;
    wavy->Geom().SetData(wavedir,numwaves);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    ///EXEMPLOS DE MAPEAMENTO GEOMETRICO ATRAVES DO METODO X
    TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
    
    std::cout << "\n\n";
    
    ///Refinando o elemento quadrilateral em subelementos (chamados de "filhos")
    TPZVec<TPZGeoEl *> sons;
    
    const int nref = 7;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malha.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile,true);
    
}

void QuadWavyMapped()
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nnodes = 4;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node0, node1, node2, node3;
    
    TPZVec<REAL> coord(3,0.);
    
    //no 0
    coord[0] = 0.;
    coord[1] = 0.;
    coord[2] = 0.;
    node0.SetNodeId(0);
    node0.SetCoord(coord);
    geomesh->NodeVec()[0] = node0;
    
    //no 1
    coord[0] = 1.;
    coord[1] = 0.;
    coord[2] = 0.;
    node1.SetNodeId(1);
    node1.SetCoord(coord);
    geomesh->NodeVec()[1] = node1;
    
    //no 2
    coord[0] = 1.;
    coord[1] = 2.;
    coord[2] = 0.;
    node2.SetNodeId(2);
    node2.SetCoord(coord);
    geomesh->NodeVec()[2] = node2;
    
    //no 3
    coord[0] = 0.;
    coord[1] = 1.;
    coord[2] = 0.;
    node3.SetNodeId(3);
    node3.SetCoord(coord);
    geomesh->NodeVec()[3] = node3;
    
    //INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
    TPZVec<long> topology(4);//Quadrilatero ira utilizar 4 nos
    topology[0] = 0;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 2;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
    
    int materialId = 1;
    
    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoQuad> > * quadrilatero =
    new TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoQuad> > (topology,materialId,*geomesh);
    
    
    //INSTANCIACAO E INICIALIZACAO DO ELEMENTO SENOIDAL
    TPZVec<long> topologyWavy1(2);//Primeiro lado da senoidal
    topologyWavy1[0] = 0;//no local 0 da senoide corresponde ao no 1 da malha geometrica
    topologyWavy1[1] = 1;//no local 1 da senoide corresponde ao no 2 da malha geometrica
    
    TPZVec<long> topologyWavy2(2);//Segundo lado da senoidal
    topologyWavy2[0] = 1;// no local 0 da senoide corresponde ao no 0 da malha geometrica
    topologyWavy2[1] = 2;// no local 1 da senoide corresponde ao no 1 da malha geometrica
    
    TPZVec<long> topologyWavy3(2);//Segundo lado da senoidal
    topologyWavy3[0] = 2;// no local 0 da senoide corresponde ao no 0 da malha geometrica
    topologyWavy3[1] = 3;// no local 1 da senoide corresponde ao no 1 da malha geometrica
    
    TPZVec<long> topologyWavy4(2);//Segundo lado da senoidal
    topologyWavy4[0] = 3;// no local 0 da senoide corresponde ao no 0 da malha geometrica
    topologyWavy4[1] = 0;// no local 1 da senoide corresponde ao no 1 da malha geometrica

   
    int materialIdW = -1;
    
    TPZGeoElRefPattern< pzgeom::TPZWavyLine > *wavy1 =
    new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy1,materialIdW,*geomesh);
    
//    TPZGeoElRefPattern< pzgeom::TPZWavyLine > *wavy2 =
//    new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy2,materialIdW,*geomesh);
//    
//    TPZGeoElRefPattern< pzgeom::TPZWavyLine > *wavy3 =
//    new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy3,materialIdW,*geomesh);
//
//    TPZGeoElRefPattern< pzgeom::TPZWavyLine > *wavy4 =
//    new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy4,materialIdW,*geomesh);
    
    TPZVec<REAL> wavedir0(3,0.);
    TPZVec<REAL> wavedir1(3,0);
    TPZVec<REAL> wavedir2(3,0);

    
    wavedir0[0] = 0.;
    wavedir0[1] = 0.03;
    
    wavedir1[0] = 0.03;
    wavedir1[1] = 0;
    
    wavedir2[0] = 0;
    wavedir2[2] = 0.03;
    
    int numwaves = 1;
    
    wavy1->Geom().SetData(wavedir2,numwaves);
//    wavy2->Geom().SetData(wavedir1,numwaves);
//    wavy3->Geom().SetData(wavedir0,numwaves);
//    wavy4->Geom().SetData(wavedir1,numwaves);
    
    ///DEFININDO QSI ETA
    TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    ///EXEMPLOS DE MAPEAMENTO GEOMETRICO ATRAVES DO METODO X

    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    //quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
    //jac.Print("Jac:");
    //axes.Print("axes");
    //std::cout << "detjac = " << detjac << std::endl;
    
    ///Refinando o elemento quadrilateral em subelementos (chamados de "filhos")
    TPZVec<TPZGeoEl *> sons;
    
    const int nref = 7;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malha_TESTE.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile,true);
    
}

void TriCircleMapped()
{
    const int nNodes = 4 + 4 + 1;
    const int matId = 1; // define id para um material(formulacao fraca)
    const int bc0 = -1;  // define id para um material(cond contorno dirichlet)
    const REAL rDomain = 1.;

    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * gmesh = new TPZGeoMesh;
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


    TPZVec<TPZGeoEl *> sons;

    const int nref = 4;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }

    std::ofstream outfile("malha_circular.vtk");
  	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outfile,true);
}

void TriWavyMapped()
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nnodes = 3;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node0, node1, node2;
    
    TPZVec<REAL> coord(3,0.);
    
    //no 0
    coord[0] = 0.;
    coord[1] = 0.;
    coord[2] = 0.;
    node0.SetNodeId(0);
    node0.SetCoord(coord);
    geomesh->NodeVec()[0] = node0;
    
    //no 1
    coord[0] = 1.;
    coord[1] = 0.;
    coord[2] = 0.;
    node1.SetNodeId(1);
    node1.SetCoord(coord);
    geomesh->NodeVec()[1] = node1;
    
    //no 2
    coord[0] = 1.;
    coord[1] = 2.;
    coord[2] = 0.;
    node2.SetNodeId(2);
    node2.SetCoord(coord);
    geomesh->NodeVec()[2] = node2;
    
    
    //INSTANCIACAO E INICIALIZACAO DO ELEMENTO TRIANGULAR
    TPZVec<long> topology(3);//Quadrilatero ira utilizar 4 nos
    topology[0] = 0;//no local 0 do triangulo corresponde ao no 0 da malha geometrica
    topology[1] = 1;//no local 1 do triangulo corresponde ao no 1 da malha geometrica
    topology[2] = 2;//no local 2 do triangulo corresponde ao no 2 da malha geometrica
    
    int materialId = 1;


    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> > * triangulo =
    new TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > > (topology,materialId,*geomesh);

    //INSTANCIACAO E INICIALIZACAO DO ELEMENTO SENOIDAL
    TPZVec<long> topologyWavy0(2);//Primeiro lado da senoidal
    topologyWavy0[0] = 0;//no local 0 da senoide corresponde ao no 0 da malha geometrica
    topologyWavy0[1] = 1;//no local 1 da senoide corresponde ao no 1 da malha geometrica
    
    TPZVec<long> topologyWavy1(2);//Segundo lado da senoidal
    topologyWavy1[0] = 1;// no local 0 da senoide corresponde ao no 1 da malha geometrica
    topologyWavy1[1] = 2;// no local 1 da senoide corresponde ao no 2 da malha geometrica
    
    TPZVec<long> topologyWavy2(2);//Segundo lado da senoidal
    topologyWavy2[0] = 2;// no local 0 da senoide corresponde ao no 2 da malha geometrica
    topologyWavy2[1] = 0;// no local 1 da senoide corresponde ao no 0 da malha geometrica
   
    
    
    int materialIdW = -1;

    TPZGeoElRefPattern< pzgeom::TPZWavyLine > *wavy0 =
    new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy0,materialIdW,*geomesh);
    
//    TPZGeoElRefPattern< pzgeom::TPZWavyLine > *wavy1 =
//    new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy1,materialIdW,*geomesh);
//    
//    TPZGeoElRefPattern< pzgeom::TPZWavyLine > *wavy2 =
//    new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (topologyWavy2,materialIdW,*geomesh);
    
    
    TPZVec<REAL> wavedir0(3,0.);
    TPZVec<REAL> wavedir1(3,0);
    TPZVec<REAL> wavedir2(3,0);
    TPZVec<REAL> wavedir3(3,0);

    
    wavedir0[0] = 0.;
    wavedir0[1] = 0.05;
    
    wavedir1[0] = 0.05;
    wavedir1[1] = 0;
    
    wavedir2[0] = 0.05;
    wavedir2[1] = 0.05;
    
    wavedir3[0] = 0;
    wavedir3[2] = 0.05;

    
    int numwaves = 3;
    
    wavy0->Geom().SetData(wavedir3,numwaves);
//    wavy1->Geom().SetData(wavedir1,numwaves);
//    wavy2->Geom().SetData(wavedir2,numwaves);
    
    ///DEFININDO QSI ETA
    TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    ///EXEMPLOS DE MAPEAMENTO GEOMETRICO ATRAVES DO METODO X
    
    triangulo -> X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    //quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
    //jac.Print("Jac:");
    //axes.Print("axes");
    //std::cout << "detjac = " << detjac << std::endl;
    
    ///Refinando o elemento quadrilateral em subelementos (chamados de "filhos")
    TPZVec<TPZGeoEl *> sons;
    
    const int nref = 7;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malha_triang.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile,true);
    
}


void CubeWavyMapped()
{
    //Instanciacao da Malha Geometrica
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nnodes = 8; //Quantidade de nos da malha
    geomesh->NodeVec().Resize(nnodes);
    
    //Inicializacao dos nos da malha
    TPZGeoNode node0, node1, node2, node3, node4, node5, node6, node7;
    
    //Inicializacao da malha geometrica pela instanciacao dos nos
    
    TPZVec<REAL> coord(3,0);
    
    
    //no 0
    coord[0] = 0.;
    coord[1] = 0.;
    coord[2] = 0.;
    node0.SetCoord(coord);
    geomesh->NodeVec()[0] = node0;
    
    //no 1
    coord[0] = 1.;
    coord[1] = 0.;
    coord[2] = 0.;
    node1.SetCoord(coord);
    geomesh->NodeVec()[1] = node1;
    
    //no 2
    coord[0] = 1.;
    coord[1] = 1.;
    coord[2] = 0.;
    node2.SetCoord(coord);
    geomesh->NodeVec()[2] = node2;
    
    //no 3
    coord[0] = 0.;
    coord[1] = 1.;
    coord[2] = 0.;
    node3.SetCoord(coord);
    geomesh->NodeVec()[3] = node3;
    
    //no 4
    coord[0] = 0.;
    coord[1] = 0.;
    coord[2] = 1.;
    node4.SetCoord(coord);
    geomesh->NodeVec()[4] = node4;
    
    //no 5
    coord[0] = 1.;
    coord[1] = 0.;
    coord[2] = 1.;
    node5.SetCoord(coord);
    geomesh->NodeVec()[5] = node5;
    
    //no 6
    coord[0] = 1.;
    coord[1] = 1.;
    coord[2] = 1.;
    node6.SetCoord(coord);
    geomesh->NodeVec()[6] = node6;
    
    //no 7
    coord[0] = 0.;
    coord[1] = 1.;
    coord[2] = 1.;
    node7.SetCoord(coord);
    geomesh->NodeVec()[7] = node7;
    
    int materialId = 1;
    
    //Instanciacao e inicializacao do elemento cubico

//    TPZVec<long> topologyCube(8);
//    topologyCube[0] = 0; topologyCube[1] = 1; topologyCube[2] = 2; topologyCube[3] = 3; //no local 0,1,2,3 correspondem aos nos 0,1,2,3 na malha global
//    topologyCube[4] = 4; topologyCube[5] = 5; topologyCube[6] = 6; topologyCube[7] = 7; //no local 4,5,6,7 correspondem aos nos 4,5,6,7 na malha global
//
//    
//    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoCube > > *cube =
//    new TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoCube > >(topologyCube,materialId,*geomesh);
    
    //Instanciacao e inicializacao do elemento quadrilateral 0
    TPZVec<long> topologyQuad0(4);
    topologyQuad0[0] = 0; //no local 0 corresponde ao no 4 da malha global
    topologyQuad0[1] = 1; //no local 1 corresponde ao no 5 da malha global
    topologyQuad0[2] = 2; //no local 2 corresponde ao no 6 da malha global
    topologyQuad0[3] = 3; //no local 3 corresponde ao no 7 da malha global
    
    //Instanciacao e inicializacao do elemento quadrilateral 1
    TPZVec<long> topologyQuad1(4);
    topologyQuad1[0] = 0; //no local 0 corresponde ao no 4 da malha global
    topologyQuad1[1] = 1; //no local 1 corresponde ao no 5 da malha global
    topologyQuad1[2] = 5; //no local 2 corresponde ao no 6 da malha global
    topologyQuad1[3] = 4; //no local 3 corresponde ao no 7 da malha global
    
    //Instanciacao e inicializacao do elemento quadrilateral 2
    TPZVec<long> topologyQuad2(4);
    topologyQuad2[0] = 1; //no local 0 corresponde ao no 4 da malha global
    topologyQuad2[1] = 2; //no local 1 corresponde ao no 5 da malha global
    topologyQuad2[2] = 6; //no local 2 corresponde ao no 6 da malha global
    topologyQuad2[3] = 5; //no local 3 corresponde ao no 7 da malha global
    
    //Instanciacao e inicializacao do elemento quadrilateral 3
    TPZVec<long> topologyQuad3(4);
    topologyQuad3[0] = 5; //no local 0 corresponde ao no 4 da malha global
    topologyQuad3[1] = 4; //no local 1 corresponde ao no 5 da malha global
    topologyQuad3[2] = 7; //no local 2 corresponde ao no 6 da malha global
    topologyQuad3[3] = 6; //no local 3 corresponde ao no 7 da malha global
    
    //Instanciacao e inicializacao do elemento quadrilateral 4
    TPZVec<long> topologyQuad4(4);
    topologyQuad4[0] = 2; //no local 0 corresponde ao no 4 da malha global
    topologyQuad4[1] = 3; //no local 1 corresponde ao no 5 da malha global
    topologyQuad4[2] = 7; //no local 2 corresponde ao no 6 da malha global
    topologyQuad4[3] = 6; //no local 3 corresponde ao no 7 da malha global
    
    //Instanciacao e inicializacao do elemento quadrilateral 5
    TPZVec<long> topologyQuad5(4);
    topologyQuad5[0] = 3; //no local 0 corresponde ao no 4 da malha global
    topologyQuad5[1] = 0; //no local 1 corresponde ao no 5 da malha global
    topologyQuad5[2] = 4; //no local 2 corresponde ao no 6 da malha global
    topologyQuad5[3] = 7; //no local 3 corresponde ao no 7 da malha global
    
    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoQuad > > *quad0 =
    new TPZGeoElMapped<TPZGeoElRefPattern<pzgeom::TPZGeoQuad> > (topologyQuad0,materialId,*geomesh);

//    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoQuad > > *quad1 =
//    new TPZGeoElMapped<TPZGeoElRefPattern<pzgeom::TPZGeoQuad> > (topologyQuad1,materialId,*geomesh);
//
//    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoQuad > > *quad2 =
//    new TPZGeoElMapped<TPZGeoElRefPattern<pzgeom::TPZGeoQuad> > (topologyQuad2,materialId,*geomesh);
//    
//    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoQuad > > *quad3 =
//    new TPZGeoElMapped<TPZGeoElRefPattern<pzgeom::TPZGeoQuad> > (topologyQuad3,materialId,*geomesh);

    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoQuad > > *quad4 =
    new TPZGeoElMapped<TPZGeoElRefPattern<pzgeom::TPZGeoQuad> > (topologyQuad4,materialId,*geomesh);

//    TPZGeoElMapped<TPZGeoElRefPattern< pzgeom::TPZGeoQuad > > *quad5 =
//    new TPZGeoElMapped<TPZGeoElRefPattern<pzgeom::TPZGeoQuad> > (topologyQuad5,materialId,*geomesh);
    
    //Instanciacao e inicializacao do elemento senoidal
    TPZVec<long> topologyWavy(2);
    topologyWavy[0] = 2; //no local 0 corresponde ao no 5 da malha global
    topologyWavy[1] = 3; //no local 0 corresponde ao no 6 da malha global
    
    int materialIdW = -1;
    
    TPZGeoElRefPattern<pzgeom::TPZWavyLine> *wavy =
    new TPZGeoElRefPattern<pzgeom::TPZWavyLine> (topologyWavy,materialIdW,*geomesh);
    
    TPZVec<REAL> wavydir (3,0);
    wavydir[0] = 0;
    wavydir[1] = 0.1;
    wavydir[2] = 0;
    int numwaves;
    numwaves = 4;
    wavy->Geom().SetData(wavydir, numwaves);
    
    //Definindo qsi, eta e zeta para a malha
    
    TPZVec<REAL> qsi (3,0), xqsi (3,0);
    
    //Construindo a malha
    geomesh -> BuildConnectivity();
    
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    
    ///Refinando o elemento quadrilateral em subelementos (chamados de "filhos")
    TPZVec<TPZGeoEl *> sons;
    
    const int nref = 7;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malha_cube.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile,true);

    
    
}

void QuadSphere()
{
	///INSTANCIACAO DA MALHA GEOMETRICA
	TPZGeoMesh * geomesh = new TPZGeoMesh;
	
	int nnodes = 6;//quantidade de nos da malha geometrica
	geomesh->NodeVec().Resize(nnodes);
	
	///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
	TPZGeoNode node0, node1, node2, node3, node4, node5;	
	TPZVec<REAL> coord(3,0.);
	const REAL r = 2.;
	TPZManVector<REAL,3> xc(3,0.);
	xc[0] = 1.;
	xc[1] = 2.;
	xc[2] = 3.;
	
	//no 0
	coord[0] = xc[0] + r;
	coord[1] = xc[1] + 0.;
	coord[2] = xc[2] + 0.;
	node0.SetNodeId(0);
	node0.SetCoord(coord);
	geomesh->NodeVec()[0] = node0;
	
	//no 1
	coord[0] = xc[0] + 0.;
	coord[1] = xc[1] + 0.;
	coord[2] = xc[2] + -r;
	node1.SetNodeId(1);
	node1.SetCoord(coord);
	geomesh->NodeVec()[1] = node1;
	
	//no 2
	coord[0] = xc[0] + 0.;
	coord[1] = xc[1] + r;
	coord[2] = xc[2] + 0.;
	node2.SetNodeId(2);
	node2.SetCoord(coord);
	geomesh->NodeVec()[2] = node2;
	
	//no 3
	coord[0] = xc[0] + 0.;
	coord[1] = xc[1] + 0.;
	coord[2] = xc[2] + r;
	node3.SetNodeId(3);
	node3.SetCoord(coord);
	geomesh->NodeVec()[3] = node3;
	
	//no 4
	coord[0] = xc[0] + 0.;
	coord[1] = xc[1] + -r;
	coord[2] = xc[2] + 0.;
	node4.SetNodeId(4);
	node4.SetCoord(coord);
	geomesh->NodeVec()[4] = node4;
	
	//no 5
	coord[0] = xc[0] + -r;
	coord[1] = xc[1] + 0.;
	coord[2] = xc[2] + 0.;
	node5.SetNodeId(5);
	node5.SetCoord(coord);
	geomesh->NodeVec()[5] = node5;
	
	int materialId = 1;
	
	// El 0
	//INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
	TPZVec<long> topology(4);//Quadrilatero ira utilizar 4 nos
	topology[0] = 0;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 2;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > * quadrilatero =
	new TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > (topology,materialId,*geomesh);
	
	quadrilatero->Geom().SetData(r,xc);
	
	// El 1
	topology[0] = 2;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 5;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > * quadrilatero2 =
	new TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > (topology,materialId,*geomesh);
	
	quadrilatero2->Geom().SetData(r,xc);
	
	// El 2	
	topology[0] = 5;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 4;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > * quadrilatero3 =
	new TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > (topology,materialId,*geomesh);
	
	quadrilatero3->Geom().SetData(r,xc);
	
	// El 3
	topology[0] = 4;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
	topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
	topology[2] = 0;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
	topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
	
	TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > * quadrilatero4 =
	new TPZGeoElRefPattern< pzgeom::TPZQuadSphere <> > (topology,materialId,*geomesh);
	
	quadrilatero4->Geom().SetData(r,xc);
	
	//CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
	geomesh->BuildConnectivity();
	
	
	///EXEMPLOS DE MAPEAMENTO GEOMETRICO ATRAVES DO METODO X
	TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
	TPZFNMatrix<9,REAL> jac(2,2),axes(2,3),jacinv(2,2);
	REAL detjac;
	
	std::cout << "\n\n";
	
	qsi[0] = +0.23;
	qsi[1] = -0.17;
	quadrilatero->X(qsi,xqsi);
	std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
	quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
	jac.Print("Jac:");
	axes.Print("axes");
	std::cout << "detjac = " << detjac << std::endl;

	qsi[0] = -0.89;
	qsi[1] = +0.31;
	quadrilatero->X(qsi,xqsi);
	std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
	quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
	jac.Print("Jac:");
	axes.Print("axes");
	std::cout << "detjac = " << detjac << std::endl;		
	
	qsi[0] = +0.05;
	qsi[1] = +1.0;
	quadrilatero->X(qsi,xqsi);
	std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
	quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
	jac.Print("Jac:");
	axes.Print("axes");
	std::cout << "detjac = " << detjac << std::endl;
	
	qsi[0] = -0.73;
	qsi[1] = -0.35;
	quadrilatero->X(qsi,xqsi);
	std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
	quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
	jac.Print("Jac:");
	axes.Print("axes");
	std::cout << "detjac = " << detjac << std::endl;
	
	std::cout << "\n\n";
	
	
	TPZVec<TPZGeoEl *> sons;
	const int nref = 4;
	for (int iref = 0; iref < nref; iref++) {
		int nel = geomesh->NElements();
		for (int iel = 0; iel < nel; iel++) {
			TPZGeoEl *gel = geomesh->ElementVec()[iel];
			if(!gel->HasSubElement()) gel->Divide(sons);
		}
	}
	
	std::ofstream outfile("malhaesfera.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile);

	
}

/** Geracao de uma malha geometrica de 1 elemento quadrilateral apenas */
void ExemploIni()
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
	TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nnodes = 6;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node0, node1, node2, node3, node4, node5;
    
    TPZVec<REAL> coord(3,0.);
    
    //no 0
    coord[0] = 3.;
    coord[1] = 2.;
    coord[2] = 3.;
    node0.SetNodeId(0);
    node0.SetCoord(coord);
    geomesh->NodeVec()[0] = node0;
    
    //no 1
    coord[0] = 7.;
    coord[1] = 1.;
    coord[2] = 4.;
    node1.SetNodeId(1);
    node1.SetCoord(coord);
    geomesh->NodeVec()[1] = node1;
    
    //no 2
    coord[0] = 8.;
    coord[1] = 6.;
    coord[2] = 2.;
    node2.SetNodeId(2);
    node2.SetCoord(coord);
    geomesh->NodeVec()[2] = node2;
    
    //no 3
    coord[0] = 4.;
    coord[1] = 4.;
    coord[2] = 6.;
    node3.SetNodeId(3);
    node3.SetCoord(coord);
    geomesh->NodeVec()[3] = node3;
    
    //no 4
    coord[0] = 9.;
    coord[1] = 3.;
    coord[2] = 6.;
    node4.SetNodeId(4);
    node4.SetCoord(coord);
    geomesh->NodeVec()[4] = node4;
    
    //no 5
    coord[0] = 8.5;
    coord[1] = 7.;
    coord[2] = 6.;
    node5.SetNodeId(5);
    node5.SetCoord(coord);
    geomesh->NodeVec()[5] = node5;
    
    //INSTANCIACAO E INICIALIZACAO DO ELEMENTO QUADRILATERAL
    TPZVec<long> topology(4);//Quadrilatero ira utilizar 4 nos
    topology[0] = 0;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 1;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 2;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    topology[3] = 3;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
    
    int materialId = 1;
    
    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * quadrilatero =
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (topology,materialId,*geomesh);
    
    
    topology[0] = 1;//no local 0 do quadrilatero corresponde ao no 0 da malha geometrica
    topology[1] = 4;//no local 1 do quadrilatero corresponde ao no 1 da malha geometrica
    topology[2] = 5;//no local 2 do quadrilatero corresponde ao no 2 da malha geometrica
    topology[3] = 2;//no local 3 do quadrilatero corresponde ao no 3 da malha geometrica
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (topology,materialId,*geomesh);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    ///EXEMPLOS DE MAPEAMENTO GEOMETRICO ATRAVES DO METODO X
    TPZVec<REAL> qsi(2,0.), xqsi(3,0.);
    
    std::cout << "\n\n";
    
    qsi[0] = +0.23;
    qsi[1] = -0.17;
    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";

    TPZFNMatrix<9,REAL> jac(2,2),axes(2,3),jacinv(2,2);
    REAL detjac;
    quadrilatero->Jacobian(qsi, jac, axes, detjac, jacinv);
    qsi[0] = -0.89;
    qsi[1] = +0.31;
    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    
    
    qsi[0] = +0.05;
    qsi[1] = +1.0;
    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    
    qsi[0] = -0.73;
    qsi[1] = -0.35;
    quadrilatero->X(qsi,xqsi);
    std::cout << "qsi = { " << qsi[0] << " , " << qsi[1] << " }  ;  x(qsi) = { " << xqsi[0] << " , " << xqsi[1] << " , " << xqsi[2] << " }\n";
    
    std::cout << "\n\n";
    
    
    ///Refinando o elemento quadrilateral em subelementos (chamados de "filhos")
    TPZVec<TPZGeoEl *> sons0, sons1;
    quadrilatero->Divide(sons0);
    
    ///Refinando os filhos do elemento quadrilateral
    for(int s = 0; s < sons0.NElements(); s++)
    {
        sons1.Resize(0);
        TPZGeoEl * actSon = sons0[s];
        actSon->Divide(sons1);
    }
    
    std::ofstream outfile("malha.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile);
    
    geomesh->Print();
    
}
