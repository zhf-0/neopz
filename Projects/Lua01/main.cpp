/**
 * @file
 * @brief Implements the build of a computational mesh as tutorial example of the interpolation NeoPZ module
 */

#include "pzcmesh.h"
#include <TPZVTKGeoMesh.h>
#include "pzelast3d.h"
#include "pzstepsolver.h"
#include "pzanalysis.h"
#include "tpzautopointer.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include <pzvec.h>
#include <pzgmesh.h>
#include <pzcompel.h>
#include <pzgeoel.h>
#include <pzquad.h>
#include <pzmat2dlin.h>
#include <TPZGeoElement.h>
#include <pzskylstrmatrix.h>
#include <pzcmesh.h>

#include "pzgeoelbc.h"

#include "pzfilebuffer.h"
#include "pzmaterialid.h"
#include "pzmeshid.h"
#include "pzbfilestream.h"
#include <pzelast3d.h>
#include <pzplaca.h>
#include <pzvtkmesh.h>
#include <pzlog.h>


#include <iostream>
#include <fstream>

using namespace std;

// nx = number of nodes in x direction
// ny = number of nodes in y direction
TPZGeoMesh * GetMesh(int nx,int ny);
void InsertElasticity(TPZCompMesh *cmesh);

const REAL Pi=M_PI;

void Forcing(const TPZVec<REAL> &x, TPZVec<STATE> &f)
{
    double q=-5.0, E=20000000,v=0.3, h=0.1, lx=5., ly=6.;
    double D=0;
    D=(E*h*h*h)/(12*(1-v*v));
    
    double x1=x[0];
    double x2=x[1];
    
    f[2] = (q)*sin(Pi*x1/lx)*sin(Pi*x2/ly);
}

// definition of v analytic

void w_exact(TPZFMatrix<REAL> &axes, TPZVec<REAL> & x, TPZFMatrix<STATE> &uexact,TPZFMatrix<STATE> &duexact){
    
    axes.Resize(3, 3);
    
    axes(0,0)=1.0;
    axes(0,1)=0.0;
    axes(0,2)=0.0;
    
    axes(1,0)=0.0;
    axes(1,1)=1.0;
    axes(1,2)=0.0;
    
    axes(2,0)=0.0;
    axes(2,1)=0.0;
    axes(2,2)=1.0;
    
    //uexact.Resize(1,1);
    
    STATE xv = x[0];
    STATE yv = x[1];
    
    STATE w_x = (-22113.*sin((Pi*(-2.5 + xv))/5.)*sin((Pi*(-3. + yv))/6.))/(37210.*pow(Pi,4.));
    
    uexact(0,0) = w_x; // x direction
    
    STATE dw_x =(-22113.*cos((Pi*(-2.5 + xv))/5.)*sin((Pi*(-3. + yv))/6.))/(186050.*pow(Pi,3));
    STATE dw_y =(-7371.*cos((Pi*(-3. + yv))/6.)*sin((Pi*(-2.5 + xv))/5.))/(74420.*pow(Pi,3));

    //duexact.Resize(2,1);
    
    duexact(0,0)=dw_x;
    duexact(1,0)=dw_y;
}

/*
 
 void Forcing1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
	/*	double x = pt[0];
 double y = pt[1];
 disp[0]= -2.*(1.-x*x) -2.*(1.-y*y);
 */
/*  double x = pt[0];
 double y = pt[1];
 disp[0]= -2.*pow(Pi,(REAL)2.)*sin(Pi*x)*sin(Pi*y);
 return;
 }
 */

int main(){
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    //	TPZSaveable::Register(TPZSAVEABLEID,Restore<TPZSaveable>);
    //
    //    RegisterMeshClasses();
    //    RegisterMatrixClasses();
    //	  RegisterMaterialClasses();
    //
    int nx = 9, ny = 9;
    int order = 3;
    TPZCompEl::SetgOrder(order);
    
    //Creates the geometric mesh
    TPZGeoMesh *mesh = GetMesh(nx,ny);
    mesh->SetName("testing a space");
    ofstream out("all.dat");
    
    // Print geometrical mesh info
    ofstream geomeshout("geomesh.txt");
    mesh->Print(geomeshout);
    
    //Link between geometrical mesh and computational mesh
    TPZCompMesh *cmesh = new TPZCompMesh(mesh);
    
    //Inserting info (dimension, polinomial order, material properties) into computational mesh
    cmesh->SetDefaultOrder(order);
    cmesh->SetDimModel(2);
    InsertElasticity(cmesh);
    cmesh->AutoBuild();
    
    // Print computational mesh info
    ofstream cmeshout("../cmesh.txt");
    cmesh->Print(cmeshout);
    
    //Assembly and building solution
    
    bool optimizeBandwidth = true;
    TPZSkylineStructMatrix skylstruct(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    TPZAnalysis an(cmesh,optimizeBandwidth);
    an.SetStructuralMatrix(skylstruct);
    an.SetSolver(step);
    //an.Run();
    
    std::cout << "Assemble matrix with NDoF = " << cmesh->NEquations() << std::endl;
    
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
    }
#endif
    
    std::cout << "Solving Matrix " << std::endl;
    
    an.Solve();
    
    
    an.Solution().Print("Solucao");
    
    
    //
    //	TPZAutoPointer<TPZCompMesh> cmeshauto(cmesh);
    //	TPZMaterial * mat = cmeshauto->FindMaterial(1);
    //	int nstate = mat->NStateVariables();
    //
    
    //an.SetStep(8);
    //Defining plot properties
    int dimension = 2, resolution = 1;
    std::string plotfile("placaaf.vtk");
    
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Deslocz");
    vecnames.Push("Displacement");
    scalnames.Push("Mn1");
    //scalnames.Push("Mn2");
    scalnames.Push("Mn1n2");
    

    int postProcessResolution = 0; //  keep low as possible
    
    int dim = mesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    
    
    
    //
    //	TPZMaterial * mat = cmesh->FindMaterial(1);
    //	int nstate = mat->NStateVariables();
    //	int nscal = 0, nvec = 0, dim = 3;
    //	if(nstate ==1)
    //	{		nscal = 1;
    //	}
    //	else
    //	{
    //		nvec = 1;
    //	}
    //	TPZManVector<std::string> scalnames(nscal),vecnames(nvec);
    //	if(nscal == 1)
    //	{
    //		scalnames[0]="state";
    //	}
    //	else
    //	{
    //		vecnames[0] = "state";
    //	}
    //	std::string postprocessname("after.vtk");
    //	TPZVTKGraphMesh vtkmesh(cmesh,dim,mat,scalnames,vecnames);
    //	vtkmesh.SetFileName(postprocessname);
    //	vtkmesh.SetResolution(1);
    //	int numcases = 1;
    //
    //
    //	// Iteracoes de tempo
    //	int istep = 0, nsteps = 2;
    //	vtkmesh.DrawMesh(numcases);
    //	vtkmesh.DrawSolution(istep, 1.);
    //
    // // Print geometrical mesh info
    //	ofstream geomeshoutf("geomeshf.txt");
    //	mesh->Print(geomeshoutf);
    //
    //	// Print computatioinal mesh info
    //	ofstream cmeshoutf("cmeshf.txt");
    //	cmesh->Print(cmeshoutf);
    //
    delete cmesh;
    delete mesh;
    return 0;
    
}

TPZGeoMesh *GetMesh (int nx,int ny){
    int i,j;
    long id, index;
    
    //Let's try with an unitary domain
    REAL lx = 5.;
    REAL ly = 6.;
    
    //Creates the geometric mesh... The nodes and elements
    //will be inserted into mesh object during initilize process
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    //Auxiliar vector to store a coordinate
    TPZVec <REAL> coord (3,0.);
    
    //Nodes initialization
    for(i = 0; i < nx; i++){
        for(j = 0; j < ny; j++){
            id = i*ny + j;
            coord[0] = (i)*lx/(nx - 1);
            coord[1] = (j)*ly/(ny - 1);
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
            if (nodecoord[2] == 0. & nodecoord[1] == 6.)
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
            if (nodecoord[2] == 0. & nodecoord[0] == 5.)
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

void InsertElasticity(TPZCompMesh *mesh)
{
    mesh->SetDimModel(3);
    int nummat = 1, dirichlet = 0, misto=2;
    int neumann = 1;
    //int mixed = 2;
    int dirbott = -1, dirtop = -2, dirleft = -3, dirright = -4;
    TPZManVector<REAL> force(3,0.);
    
    force[2] = 1.;
    
    REAL h = 0.1;
    REAL f = 0;
    REAL E = 2.*10000000;
    REAL ni1 = 0.3;
    REAL ni2 = 0.3;
    REAL G = E/(2*(1+ni1));
    //	TPZFMatrix<REAL> naxes;
    TPZFMatrix<STATE> naxes(3,3,0.0);
    
    naxes(0,0)=1.0;
    naxes(0,1)=0.0;
    naxes(0,2)=0.0;
    
    naxes(1,0)=0.0;
    naxes(1,1)=1.0;
    naxes(1,2)=0.0;
    
    naxes(2,0)=0.0;
    naxes(2,1)=0.0;
    naxes(2,2)=1.0;
    
    ofstream file("axes.txt");
    naxes.Print(" axes = ",file, EMathematicaInput);
    
    //#ifdef log4cxx
    //	if (logdata->isDebugEnabled()) {
    //
    //		std::stringstream sout;
    //		naxes.Print(" axes ",sout, EMathematicaInput );
    //		LOGPZ_DEBUG(lodgata,sout.str());
    //	}
    //#endif
    
    TPZVec< STATE > xf(6,0.0);
    xf[2]=0;
    
    
    //TPZMatPoisson3d *p = new TPZMatPoisson3d(nummat, 2);//
    
    //TPZElasticity3D *elast = new TPZElasticity3D(nummat, Ela, poisson, force);//
    
    
    TPZPlaca *placa = new TPZPlaca(nummat, h, f, E, E, ni1, ni2, G, G, G, naxes, xf);
    
    TPZAutoPointer<TPZFunction<STATE> > minha = new TPZDummyFunction<STATE>(Forcing);
    
    
    //Print(std::ostream &out).Print("axes");
    
    /*
     TPZAutoPointer<TPZDummyFunction<STATE> > minha(Forcing);
     */
    
    /// brincar com execute de minha
    //minha->Execute(xf, TPZVec<double> &f, TPZFMatrix<double> &df);
    
    placa->SetForcingFunction(minha);
    
    placa->SetExactFunction(w_exact);
    
    TPZMaterial * elastauto(placa);
    //elastauto->SetForcingFunction(minha);
    
    mesh->InsertMaterialObject(elastauto);
    /*
     TPZFMatrix<STATE> val1bott(6,6,0.),val2bott(6,1,0.);
     TPZBndCond *bcbott = placa->CreateBC(elastauto, dirbott, dirichlet, val1bott, val2bott);
     TPZMaterial * bcbottauto(bcbott);
     mesh->InsertMaterialObject(bcbottauto);
     */
    
    TPZFMatrix<STATE> val1bott(6,6,0.),val2bott(6,1,0.);
    
    val1bott(0,0) = 1.e19;
    val1bott(1,1) = 1.e19;
    val1bott(2,2) = 1.e19;
    
    TPZBndCond *bcbott = placa->CreateBC(elastauto, dirbott, misto, val1bott, val2bott);
    TPZMaterial * bcbottauto(bcbott);
    mesh->InsertMaterialObject(bcbottauto);
    
    TPZFMatrix<STATE> val1top(6,6,0.),val2top(6,1,0.);
    
    val1top(0,0) = 1.e19;
    val1top(1,1) = 1.e19;
    val1top(2,2) = 1.e19;
    //    val2top(3.,0.) = 0;
    TPZBndCond *bctop = placa->CreateBC(elastauto, dirtop, misto, val1top, val2top);
    TPZMaterial * bctopauto(bctop);
    mesh->InsertMaterialObject(bctopauto);
    
    TPZFMatrix<STATE> val1left(6,6,0.),val2left(6,1,0.);
    
    val1left(0,0) = 1.e19;
    val1left(1,1) = 1.e19;
    val1left(2,2) = 1.e19;
    //    val2left(3.,0.) = 0;
    TPZBndCond *bcleft = placa->CreateBC(elastauto, dirleft, misto, val1left, val2left);
    TPZMaterial * bcleftauto(bcleft);
    mesh->InsertMaterialObject(bcleftauto);
    
    TPZFMatrix<STATE> val1right(6,6,0.),val2right(6,1,0.);
    
    val1right(0,0) = 1.e19;
    val1right(1,1) = 1.e19;
    val1right(2,2) = 1.e19;
    //    val2right(3.,0.) = 0;
    TPZBndCond *bcright = placa->CreateBC(elastauto, dirright, misto, val1right, val2right);
    TPZMaterial * bcrightauto(bcright);
    mesh->InsertMaterialObject(bcrightauto);
    
    
    
    /*
     // Dirichlet em 1 -1 -1 yz;
     val1(0,0) = 0.;
     val1(1,1) = 1.;
     val1(2,2) = 1.;
     TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
     TPZMaterial * bcauto2(bc2);
     mesh->InsertMaterialObject(bcauto2);
     
     // Dirichlet em 1 1 -1 z;
     val1(0,0) = 0.;
     val1(1,1) = 0.;
     val1(2,2) = 1.;
     TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
     TPZMaterial * bcauto3(bc3);
     mesh->InsertMaterialObject(bcauto3);
     */
}

