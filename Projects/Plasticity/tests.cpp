//$Id: poroelastoplastic.cpp,v 1.53 2010-06-11 22:13:02 diogo Exp $


/***************************************************************************
 *   Copyright (C) 2008 by Erick Slis   *
 *      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
/*
 using namespace std;
 
 int main(int argc, char *argv[])
 {
 cout << "Hello, world!" << endl;
 
 return EXIT_SUCCESS;
 }
 */

#include "pzelctemp.h" // TPZIntelGen<TSHAPE>
#include "pzshapecube.h" // TPZShapeCube
#include "pzcompelwithmem.h"
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZLadeKim.h"
#include "TPZSandlerDimaggio.h"
#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzporoanalysis.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"

#include "TPZTensor.h"

#include "pzelast3d.h"
#include "pzbndcond.h"

#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"

#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "pzbdstrmatrix.h"
#include "pzstepsolver.h"
#include "pzdiffmatrix.h"


#include "TPZPlasticStep.h"
#include "TPZElasticResponse.h"
#include "TPZLadeNelsonElasticResponse.h"
#include "TPZThermoForceA.h"
#include "TPZYCVonMises.h"
#include "TPZYCTresca.h"
#include "tpzyctrescaregularized.h"
#include "tpzycvonmisescombtresca.h"
#include "TPZYCLadeKim.h"
#include "TPZLadeKimThermoForceA.h"
#include "TPZLadeKim.h"
#include "TPZYCSandlerDimaggio.h"
#include "TPZSandlerDimaggio.h"
#include "TPZSandlerDimaggioThermoForceA.h"
#include "TPZMohrCoulomb.h"
//#include "checkconv.h"
#include <iostream>
#include <string>
#include "pzmatrix.h"
#include "TPZTensor.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include <math.h>
#include "pzelasmat.h"
#include "pzelast3d.h"
#include "TPZMohrCoulomb.h"

#include "TPZTensor.h"
#include "TPZYCDruckerPrager.h"
#include "pzelastoplasticanalysis.h"
#include "pzelastoplastic.h"
#include "TPZSandlerDimaggio.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"


#include "TPZDruckerPrager.h"


#include "TPZYCMohrCoulomb.h"
#include "TPZYCWillamWarnke.h"

#include "TPZSpStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZTimer.h"


#include "pzelastoplastic2d.h"


#include "pzvec.h"

#include "pzcmesh.h"

#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"

#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"

#include "pzsbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzadmchunk.h"
#include "gmres.h"
#include "pzbndcond.h"
#include "pzelast3d.h"
#include "pzblockdiag.h"
#include "pzvisualmatrix.h"
#include "cg.h"
#include "pzseqsolver.h"

#include "TPZYCVonMises.h"
#include "TPZVonMises.h"
#include "TPZParSkylineStructMatrix.h"


#include "TPZParFrontMatrix.h"
#include "TPZFrontMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"

#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"
#include "BrazilianTestGeoMesh.h"
#include "GeoMeshClass.h"




#include <sstream>

using namespace pzshape; // needed for TPZShapeCube and related classes
void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
                            int BCId,int BCId2, int nsteps, REAL PGRatio,
                            TPZFMatrix & val1Begin, TPZFMatrix & val1End,
                            TPZFMatrix & val2Begin, TPZFMatrix & val2End);

void SolveLinear(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    //TPZParSkylineStructMatrix full(fCmesh,2);
    TPZSkylineStructMatrix full(fCmesh);
    an.SetStructuralMatrix(full);
    TPZStepSolver step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
   
}

void CMeshPressuredCilinder(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
{
   
    TPZFMatrix kk(2,2,0.);
    TPZFMatrix ff(2,1,0.);
    ff(0,0)=1.;
    TPZAutoPointer<TPZMaterial> ContBC = mat->CreateBC(mat, -3, 5, kk, ff);
    CMESH->InsertMaterialObject(ContBC);
   
    TPZFMatrix k1(2,2,0.);
    TPZFMatrix f1(2,1,0.);
    //        k1(0,0)=BIG;
    //        k1(1,1)=BIG;
    f1(0,0)=1.;
    TPZAutoPointer<TPZMaterial> ContBC1 = mat->CreateBC(mat, -2,3, k1, f1);
    CMESH->InsertMaterialObject(ContBC1);
   
   
    TPZFMatrix k2(2,2,0.);
    TPZFMatrix f2(2,1,0.);
    //        k2(0,0)=BIG;
    //        k2(1,1)=BIG;
    f2(1,0)=1.;
    TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -1, 3, k2, f2);
    CMESH->InsertMaterialObject(ContBC2);
   
   
    CMESH->AutoBuild();
}




void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
                            int BCId,int BCId2, int nsteps, REAL PGRatio,
                            TPZFMatrix & val1Begin, TPZFMatrix & val1End,
                            TPZFMatrix & val2Begin, TPZFMatrix & val2End)
{
   
   
    if(!analysis.Mesh())return;
   
    // computing the initial value for the PG progression such that its sum equals one;
    REAL a0;
   
    if(fabs(PGRatio - 1.) < 1.e-3)
    {
        a0 = 1. / REAL(nsteps);
    }
   
    else
    {
        a0 = (PGRatio - 1) / (pow(PGRatio,nsteps) - 1.);
    }
    TPZFNMatrix<36> val1(6,6,0.), deltaVal1(6,6,0.);
    TPZFNMatrix< 6> val2(6,1,0.), deltaVal2(6,1,0.);
   
    deltaVal1 = val1End;
    deltaVal1.ZAXPY(-1., val1Begin);
    deltaVal2 = val2End;
    deltaVal2.ZAXPY(-1., val2Begin);   
   
    //-19
    TPZAutoPointer<TPZMaterial> mat = analysis.Mesh()->FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat.operator->());
    if(!pBC)return;
   
    int i;
    for(i = 0; i < nsteps; i++)
    {
        REAL stepLen;
        if(fabs(PGRatio - 1.) < 1.e-3)
        {
            stepLen = REAL(i+1) / REAL(nsteps);
        }
        else
        {
            stepLen = a0 * (pow(PGRatio,i+1) - 1) / (PGRatio - 1.);
        }
       
        val1 = val1Begin;
        val1.ZAXPY(stepLen, deltaVal1);
        val2 = val2Begin;
        val2.ZAXPY(stepLen, deltaVal2);
       
        pBC->Val1() = val1;
        pBC->Val2() = val2;
       
        TPZTimer time;
        time.start();
        analysis.IterativeProcess(out, tol, numiter);
        time.stop();
        REAL timeelapsed = time.seconds();
        cout << " \n TEMPO  GASTO NO ITERATIVE PROCESS "<<timeelapsed<< endl;
       
        time.reset();
        time.start();
        analysis.AcceptSolution();
        time.stop();
        REAL timeelapsed2 = time.seconds();
        cout << " \n TEMPO  GASTO NO AcceptSolution "<<timeelapsed2<< endl;
       
    }
   
}



void PressureCilinder()
{
   
    TPZFMatrix BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
    TPZFMatrix val1(3,1,0.);TPZFMatrix val2(3,1,0.);TPZFMatrix BeginForce(3,1,0.);TPZFMatrix EndForce(3,1,0.);
   
    int BC1,BC2,nsteps,taxa,nnewton;
    int h=3;
    int order = 2;
    REAL tol = 1.e-5;
    nnewton = 10;
    BC1 = -3;
    nsteps =1;
    taxa = 1;
    BeginForce(0,0) = 0.;
    EndForce(0,0) = 0.25;
   
   
    TPZGeoMesh * MESH = new TPZGeoMesh;
   
    int refdir = 3;
    MESH = BrazilianTestGeoMesh::MisesPressure(h,refdir);
   
    ofstream arg1("GeoMesh.txt");
    MESH->Print(arg1);
    TPZCompEl::SetgOrder(order);
    TPZCompMesh *CMESH = new TPZCompMesh(MESH);
   
    //TPZDruckerPrager * Pstep = new TPZDruckerPrager();
    TPZVonMises * Pstep = new TPZVonMises();
   
    REAL Yield = 0.24;//GPA
    Pstep->fTFA.SetUp(Yield,1000.);
    Pstep->fER.SetUp(/*young GPA */ 210., /*poisson*/ 0.30);
   
    TPZMatElastoPlastic2D<TPZVonMises> EPMat2(1,1);
    EPMat2.SetPlasticity(*Pstep);
    TPZAutoPointer<TPZMaterial> plastic(&EPMat2);
    plastic->Print(cout);
    CMESH->InsertMaterialObject(plastic);
   
    TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem();
   
    CMeshPressuredCilinder(CMESH,plastic);   
   
    ofstream arg("CMESHPLASTIC2D.txt");
    CMESH->Print(arg);
   
   
   
    TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
    TPZSkylineStructMatrix structmatrix(EPAnalysis.Mesh());
    structmatrix.SetNumThreads(2);
   
   
    SolveLinear(EPAnalysis,CMESH);
   
    //    TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
    //    TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
    //    PPAnalysis.SetStructuralMatrix(structmatrix);
    //    TPZVec<int> PostProcMatIds(1,1);
    //    TPZVec<std::string> PostProcVars, scalNames, vecNames;
    //    SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
    //    PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
   
    //    EPAnalysis.TransferSolution(PPAnalysis);
   
    cout << "\nDefining Graph Mesh\n";
    int dimension =2;
   
    //    PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"TuboNewPostProcess3.vtk");
   
    cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
   
    //    PPAnalysis.PostProcess(0/*pOrder*/);
   
    ManageIterativeProcess(EPAnalysis,cout,tol,nnewton,BC1,BC2,nsteps,taxa,BeginStress,EndStress,BeginForce,EndForce);//,&PPAnalysis,0);
   
    //cout << "\nInitial Solution Exported. Solving Problem\n";
    //        EPAnalysis.IterativeProcess(cout, tol, nnewton);
    //        cout << EPAnalysis.Solution() << endl;
    //        EPAnalysis.AcceptSolution();
    //        EPAnalysis.TransferSolution(PPAnalysis);
    //        PPAnalysis.PostProcess(0);
   
   
    //    PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"TuboNewPostProcess3.vtk");
    //    PPAnalysis.PostProcess(0);
    //    PPAnalysis.CloseGraphMesh();
   
   
}

int main()
{

    PressureCilinder();
    return 0;
}
