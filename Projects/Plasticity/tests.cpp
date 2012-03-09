// //$Id: poroelastoplastic.cpp,v 1.53 2010-06-11 22:13:02 diogo Exp $
// 
// 
// /***************************************************************************
//  *   Copyright (C) 2008 by Erick Slis   *
//  *      *
//  *                                                                         *
//  *   This program is free software; you can redistribute it and/or modify  *
//  *   it under the terms of the GNU General Public License as published by  *
//  *   the Free Software Foundation; either version 2 of the License, or     *
//  *   (at your option) any later version.                                   *
//  *                                                                         *
//  *   This program is distributed in the hope that it will be useful,       *
//  *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
//  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
//  *   GNU General Public License for more details.                          *
//  *                                                                         *
//  *   You should have received a copy of the GNU General Public License     *
//  *   along with this program; if not, write to the                         *
//  *   Free Software Foundation, Inc.,                                       *
//  *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
//  ***************************************************************************/
// 
// 
// #ifdef HAVE_CONFIG_H
// #include <config.h>
// #endif
// 
// #include <iostream>
// #include <cstdlib>
// /*
//  using namespace std;
//  
//  int main(int argc, char *argv[])
//  {
//  cout << "Hello, world!" << endl;
//  
//  return EXIT_SUCCESS;
//  }
//  */
// 
// #include "pzelctemp.h" // TPZIntelGen<TSHAPE>
// #include "pzshapecube.h" // TPZShapeCube
// #include "pzcompelwithmem.h"
// #include "pzelastoplastic.h"
// #include "pzporous.h"
// #include "TPZLadeKim.h"
// #include "TPZSandlerDimaggio.h"
// #include "TPZYCDruckerPrager.h"
// #include "TPZThermoForceA.h"
// #include "TPZElasticResponse.h"
// #include "pzelastoplasticanalysis.h"
// #include "pzporoanalysis.h"
// #include "pzanalysis.h"
// #include "pzfstrmatrix.h"
// #include "pzskylstrmatrix.h"
// #include "TPZParFrontStructMatrix.h"
// #include "TPZParFrontMatrix.h"
// #include "TPZFrontNonSym.h"
// #include "pzbdstrmatrix.h"
// #include "pzblockdiag.h"
// #include "TPZSpStructMatrix.h"
// 
// #include "TPZTensor.h"
// 
// #include "pzelast3d.h"
// #include "pzbndcond.h"
// 
// #include "pzcompelpostproc.h"
// #include "pzpostprocmat.h"
// #include "pzpostprocanalysis.h"
// 
// #include "pzblockdiag.h"
// #include "TPZSpStructMatrix.h"
// #include "pzbdstrmatrix.h"
// #include "pzstepsolver.h"
// #include "pzdiffmatrix.h"
// 
// 
// #include "TPZPlasticStep.h"
// #include "TPZElasticResponse.h"
// #include "TPZLadeNelsonElasticResponse.h"
// #include "TPZThermoForceA.h"
// #include "TPZYCVonMises.h"
// #include "TPZYCTresca.h"
// #include "tpzyctrescaregularized.h"
// #include "tpzycvonmisescombtresca.h"
// #include "TPZYCLadeKim.h"
// #include "TPZLadeKimThermoForceA.h"
// #include "TPZLadeKim.h"
// #include "TPZYCSandlerDimaggio.h"
// #include "TPZSandlerDimaggio.h"
// #include "TPZSandlerDimaggioThermoForceA.h"
// #include "TPZMohrCoulomb.h"
// //#include "checkconv.h"
// #include <iostream>
// #include <string>
// #include "pzmatrix.h"
// #include "TPZTensor.h"
// #include "pzvec.h"
// #include "pzfmatrix.h"
// #include "pzanalysis.h"
// #include "TPZParSkylineStructMatrix.h"
// #include "pzfstrmatrix.h"
// #include "pzstepsolver.h"
// #include "pzcmesh.h"
// #include "tpzcompmeshreferred.h"
// #include "pzpoisson3d.h"
// #include "pzbndcond.h"
// #include <math.h>
// #include "pzelasmat.h"
// #include "pzelast3d.h"
// #include "TPZMohrCoulomb.h"
// 
// #include "TPZTensor.h"
// #include "TPZYCDruckerPrager.h"
// #include "pzelastoplasticanalysis.h"
// #include "pzelastoplastic.h"
// #include "TPZSandlerDimaggio.h"
// #include "TPZThermoForceA.h"
// #include "TPZElasticResponse.h"
// 
// 
// #include "TPZDruckerPrager.h"
// 
// 
// #include "TPZYCMohrCoulomb.h"
// #include "TPZYCWillamWarnke.h"
// 
// #include "TPZSpStructMatrix.h"
// #include "TPZParSkylineStructMatrix.h"
// #include "pzstepsolver.h"
// #include "TPZTimer.h"
// 
// 
// #include "pzelastoplastic2D.h"
// 
// 
// #include "pzvec.h"
// 
// #include "pzcmesh.h"
// 
// #include "pzdebug.h"
// #include "pzcheckgeom.h"
// 
// #include "pzgeoel.h"
// #include "pzgnode.h"
// #include "pzgeoelside.h"
// 
// #include "pzintel.h"
// #include "pzcompel.h"
// 
// #include "pzmatrix.h"
// 
// #include "pzanalysis.h"
// #include "pzfstrmatrix.h"
// #include "pzskylstrmatrix.h"
// #include "TPZParFrontStructMatrix.h"
// #include "TPZParFrontMatrix.h"
// #include "TPZFrontNonSym.h"
// #include "pzbdstrmatrix.h"
// #include "pzblockdiag.h"
// #include "TPZSpStructMatrix.h"
// 
// #include "TPZCopySolve.h"
// #include "TPZStackEqnStorage.h"
// 
// #include "pzsbstrmatrix.h"
// #include "pzstepsolver.h"
// 
// #include "pzadmchunk.h"
// #include "gmres.h"
// #include "pzbndcond.h"
// #include "pzelast3d.h"
// #include "pzblockdiag.h"
// #include "pzvisualmatrix.h"
// #include "cg.h"
// #include "pzseqsolver.h"
// 
// #include "TPZYCVonMises.h"
// #include "TPZVonMises.h"
// #include "TPZParSkylineStructMatrix.h"
// 
// 
// #include "TPZParFrontMatrix.h"
// #include "TPZFrontMatrix.h"
// #include "TPZParFrontStructMatrix.h"
// #include "TPZFrontStructMatrix.h"
// 
// #include "TPZFrontSym.h"
// #include "TPZFrontNonSym.h"
// #include "BrazilianTestGeoMesh.h"
// #include "GeoMeshClass.h"
// 
// 
// 
// 
// #include <sstream>
// 
// 
// 
// using namespace pzshape; // needed for TPZShapeCube and related classes
// void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
//                             int BCId,int BCId2, int nsteps, REAL PGRatio,
//                             TPZFMatrix & val1Begin, TPZFMatrix & val1End,
//                             TPZFMatrix & val2Begin, TPZFMatrix & val2End);
// 
// void SolveLinear(TPZAnalysis &an, TPZCompMesh *fCmesh)
// {
//     //TPZParSkylineStructMatrix full(fCmesh,2);
//     TPZSkylineStructMatrix full(fCmesh);
//     an.SetStructuralMatrix(full);
//     TPZStepSolver step;
//     step.SetDirect(ELDLt);
//     an.SetSolver(step);
//    
// }
// 
// void CMeshPressuredCilinder(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
// {
//    
//     TPZFMatrix kk(2,2,0.);
//     TPZFMatrix ff(2,1,0.);
//     ff(0,0)=1.;
//     TPZAutoPointer<TPZMaterial> ContBC = mat->CreateBC(mat, -3, 5, kk, ff);
//     CMESH->InsertMaterialObject(ContBC);
//    
//     TPZFMatrix k1(2,2,0.);
//     TPZFMatrix f1(2,1,0.);
//     //        k1(0,0)=BIG;
//     //        k1(1,1)=BIG;
//     f1(0,0)=1.;
//     TPZAutoPointer<TPZMaterial> ContBC1 = mat->CreateBC(mat, -2,3, k1, f1);
//     CMESH->InsertMaterialObject(ContBC1);
//    
//    
//     TPZFMatrix k2(2,2,0.);
//     TPZFMatrix f2(2,1,0.);
//     //        k2(0,0)=BIG;
//     //        k2(1,1)=BIG;
//     f2(1,0)=1.;
//     TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -1, 3, k2, f2);
//     CMESH->InsertMaterialObject(ContBC2);
//    
//    
//     CMESH->AutoBuild();
// }
// 
// 
// 
// 
// void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
//                             int BCId,int BCId2, int nsteps, REAL PGRatio,
//                             TPZFMatrix & val1Begin, TPZFMatrix & val1End,
//                             TPZFMatrix & val2Begin, TPZFMatrix & val2End)
// {
//    
//    
//     if(!analysis.Mesh())return;
//    
//     // computing the initial value for the PG progression such that its sum equals one;
//     REAL a0;
//    
//     if(fabs(PGRatio - 1.) < 1.e-3)
//     {
//         a0 = 1. / REAL(nsteps);
//     }
//    
//     else
//     {
//         a0 = (PGRatio - 1) / (pow(PGRatio,nsteps) - 1.);
//     }
//     TPZFNMatrix<36> val1(6,6,0.), deltaVal1(6,6,0.);
//     TPZFNMatrix< 6> val2(6,1,0.), deltaVal2(6,1,0.);
//    
//     deltaVal1 = val1End;
//     deltaVal1.ZAXPY(-1., val1Begin);
//     deltaVal2 = val2End;
//     deltaVal2.ZAXPY(-1., val2Begin);   
//    
//     //-19
//     TPZAutoPointer<TPZMaterial> mat = analysis.Mesh()->FindMaterial(BCId);
//     TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat.operator->());
//     if(!pBC)return;
//    
//     int i;
//     for(i = 0; i < nsteps; i++)
//     {
//         REAL stepLen;
//         if(fabs(PGRatio - 1.) < 1.e-3)
//         {
//             stepLen = REAL(i+1) / REAL(nsteps);
//         }
//         else
//         {
//             stepLen = a0 * (pow(PGRatio,i+1) - 1) / (PGRatio - 1.);
//         }
//        
//         val1 = val1Begin;
//         val1.ZAXPY(stepLen, deltaVal1);
//         val2 = val2Begin;
//         val2.ZAXPY(stepLen, deltaVal2);
//        
//         pBC->Val1() = val1;
//         pBC->Val2() = val2;
//        
//         TPZTimer time;
//         time.start();
//         analysis.IterativeProcess(out, tol, numiter);
//         time.stop();
//         REAL timeelapsed = time.seconds();
//         cout << " \n TEMPO  GASTO NO ITERATIVE PROCESS "<<timeelapsed<< endl;
//        
//         time.reset();
//         time.start();
//         analysis.AcceptSolution();
//         time.stop();
//         REAL timeelapsed2 = time.seconds();
//         cout << " \n TEMPO  GASTO NO AcceptSolution "<<timeelapsed2<< endl;
//        
//     }
//    
// }
// 
// 
// 
// void PressureCilinder()
// {
//    
//     TPZFMatrix BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
//     TPZFMatrix val1(3,1,0.);TPZFMatrix val2(3,1,0.);TPZFMatrix BeginForce(3,1,0.);TPZFMatrix EndForce(3,1,0.);
//    
//     int BC1,BC2,nsteps,taxa,nnewton;
//     int h=3;
//     int order = 2;
//     REAL tol = 1.e-5;
//     nnewton = 10;
//     BC1 = -3;
//     nsteps =1;
//     taxa = 1;
//     BeginForce(0,0) = 0.;
//     EndForce(0,0) = 0.25;
//    
//    
//     TPZGeoMesh * MESH = new TPZGeoMesh;
//    
//     int refdir = 3;
//     MESH = BrazilianTestGeoMesh::MisesPressure(h,refdir);
//    
//     ofstream arg1("GeoMesh.txt");
//     MESH->Print(arg1);
//     TPZCompEl::SetgOrder(order);
//     TPZCompMesh *CMESH = new TPZCompMesh(MESH);
//     TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(CMESH);
//     //TPZDruckerPrager * Pstep = new TPZDruckerPrager();
//     TPZVonMises * Pstep = new TPZVonMises();
//    
//     REAL Yield = 0.24;//GPA
//     Pstep->fTFA.SetUp(Yield,1000.);
//     Pstep->fER.SetUp(/*young GPA */ 210., /*poisson*/ 0.30);
//    
//     TPZMatElastoPlastic2D<TPZVonMises> EPMat2(1,1);
//     EPMat2.SetPlasticity(*Pstep);
//     TPZAutoPointer<TPZMaterial> plastic(&EPMat2);
//     plastic->Print(cout);
//     CMESH->InsertMaterialObject(plastic);
//    
//     //TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(CMESH);
//    
//     CMeshPressuredCilinder(CMESH,plastic);   
//    
//     ofstream arg("CMESHPLASTIC2D.txt");
//     CMESH->Print(arg);
//    
//    
//    
//     TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
//     TPZFStructMatrix structmatrix(EPAnalysis.Mesh());
//     structmatrix.SetNumThreads(20);
//    
//    
//     SolveLinear(EPAnalysis,CMESH);
//    
//     //    TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
//     //    TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
//     //    PPAnalysis.SetStructuralMatrix(structmatrix);
//     //    TPZVec<int> PostProcMatIds(1,1);
//     //    TPZVec<std::string> PostProcVars, scalNames, vecNames;
//     //    SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
//     //    PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
//    
//     //    EPAnalysis.TransferSolution(PPAnalysis);
//    
//     cout << "\nDefining Graph Mesh\n";
//     int dimension =2;
//    
//     //    PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"TuboNewPostProcess3.vtk");
//    
//     cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
//    
//     //    PPAnalysis.PostProcess(0/*pOrder*/);
//    
//     ManageIterativeProcess(EPAnalysis,cout,tol,nnewton,BC1,BC2,nsteps,taxa,BeginStress,EndStress,BeginForce,EndForce);//,&PPAnalysis,0);
//    
//     //cout << "\nInitial Solution Exported. Solving Problem\n";
//     //        EPAnalysis.IterativeProcess(cout, tol, nnewton);
//     //        cout << EPAnalysis.Solution() << endl;
//     //        EPAnalysis.AcceptSolution();
//     //        EPAnalysis.TransferSolution(PPAnalysis);
//     //        PPAnalysis.PostProcess(0);
//    
//    
//     //    PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"TuboNewPostProcess3.vtk");
//     //    PPAnalysis.PostProcess(0);
//     //    PPAnalysis.CloseGraphMesh();
//    
//    
// }
// 
// void SandlerDimaggioIsotropicCompression()
// {
//     ofstream outfiletxt5("SDMcCormicRanchSand.txt");
//     TPZTensor<REAL> stress, strain, deltastress, deltastrain;
//     TPZFNMatrix<6*6> Dep(6,6,0.);
//    
//    
//     TPZSandlerDimaggio SD;
// 
//     cout << "\n Put the value of strain you want to add in each step of your loat test: ";
//     REAL straininput;
//     cin >> straininput;
//     deltastrain.XX() = -straininput;
//     deltastrain.XY() = 0.;
//     deltastrain.XZ() = 0.;
//     deltastrain.YY() = -straininput;
//     deltastrain.YZ() = 0.;
//     deltastrain.ZZ() = -straininput;
//     strain=deltastrain;   
//     cout << "Choose the material pareameters you want to set to SandlerDimaggio Test :";
//     cout << "\n0 - McCormicRanchSandMod";
//     cout << "\n1 - McCormicRanchSandMod2 ";
//     cout << "\n2 - UncDeepSandRes ";
//     cout << "\n3 - UncDeepSandResPSI";
//     cout << "\n4 - UncDeepSandResMPa";
//     cout << "\n5 - Put the material parameters you want ";
//     int choice;
//     cin >> choice;
//    
//     cout << "\n Put the numbers of steps you want: ";
//     int length;
//     cin >> length;
// 
//     switch (choice) {
//         case(0):
//         {
//             TPZSandlerDimaggio::McCormicRanchSandMod(SD);
//             std::ofstream outfiletxt("TPZSandlerDimaggioMcCormicRanchSandMod(SD).txt");
//            
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 SD.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//             }
//             break;
//         }
//         case(1):
//         {
//             TPZSandlerDimaggio::McCormicRanchSandMod2(SD);
//             std::ofstream outfiletxt("TPZSandlerDimaggioMcCormicRanchSandMod2.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 SD.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//         }
//         case(2):
//         {
//             TPZSandlerDimaggio::UncDeepSandRes(SD);
//             std::ofstream outfiletxt("TPZLadeKimUncDeepSandRes.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 SD.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//         }
//         case(3):
//         {
//             TPZSandlerDimaggio::UncDeepSandResPSI(SD);
//             std::ofstream outfiletxt("TPZLadeKimUncDeepSandResPSI.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 SD.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//         }
//         case(4):
//         {
//             TPZSandlerDimaggio::UncDeepSandResMPa(SD);
//             std::ofstream outfiletxt("TPZLadeKimUncDeepSandResMPa.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 SD.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//        
//             break;
//            
//         }   
//         case(5):
//         {
//             cout<< "\n Young Modulus ";
//             REAL E;
//             cin >> E;
//            
//             cout<< "\n poisson ";
//             REAL poisson;
//             cin >> poisson;
//            
//            
//             SD.fER.SetUp(E, poisson);
//            
//             cout<< "\n A ";
//             REAL A;
//             cin >> A;
//            
//             cout << "\n B ";
//             REAL B;
//             cin >> B;
//            
//             cout << "\n C ";
//             REAL C;
//             cin >> C;
//            
//             cout << "\n D ";
//             REAL D;
//             cin >> D;
//            
//             cout << "\n R ";
//             REAL R;
//             cin >> R;
//            
//             cout << "\n W ";
//             REAL W;
//             cin >> W;
//            
//             SD.fYC.SetUp(A, B, C, D, R, W);
//    
//             std::ofstream outfiletxt("SandlerDimaggioYOURMODEL.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 SD.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//         }
//         default:
//         {
//             cout << "Unknown Test Type. Exiting...";
//             break;
//         }
//     }
//    
//    
//    
// }
// 
// void LKFineSilicaLoadTest()
// {
//     ofstream outfiletxt1("FineSilica.txt");
//     TPZLadeKim LK;
//     TPZLadeKim::FineSilicaSand(LK);
//     TPZTensor<REAL> stress, strain, deltastress, deltastrain;
//     TPZFNMatrix<6*6> Dep(6,6,0.);
//     deltastress.XX() = -4.;
//     deltastress.XY() = 0.;
//     deltastress.XZ() = 0.;
//     deltastress.YY() = -4.;
//     deltastress.YZ() = 0.;
//     deltastress.ZZ() = -4.;
//     stress = deltastress;
//    
//    
//     int length =30;
//     for(int step=0;step<length;step++)
//     {
//         cout << "\nstep "<< step;
//         if(step == 69)deltastress *=-1.;
//         LK.ApplyLoad(stress,strain);           
//         outfiletxt1 << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//         stress +=  deltastress;
//         cout << "strain = "<<strain <<"\n";
//         cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
//        
//     }
//    
//    
//    
//    
// }
// 
// void LKIsotropicCompression()
// {
//    
//     TPZTensor<REAL> stress, strain, deltastress, deltastrain;
//     TPZFNMatrix<6*6> Dep(6,6,0.);
//    
//     cout << "\n Put the value of strain you want to add in each step of your loat test: ";
//     REAL straininput;
//     cin >> straininput;
//    
//     deltastrain.XX() = -straininput;
//     deltastrain.XY() = 0.;
//     deltastrain.XZ() = 0.;
//     deltastrain.YY() = -straininput;
//     deltastrain.YZ() = 0.;
//     deltastrain.ZZ() = -straininput;
//     strain=deltastrain;
//    
//     cout << "Choose the material pareameters you want to set to Lade Kim Test :";
//     cout << "\n0 - Plain Concrete ";
//     cout << "\n1 - Loose Sacramento River Sand ";
//     cout << "\n2 - Dense Sacramento River Sand ";
//     cout << "\n3 - Fine Silica Sand";
//     cout << "\n4 - Put the material parameters you want";
//     int choice;
//     cin >> choice;
//    
//     cout << "\n Put the numbers of steps you want: ";
//     int length;
//     cin >> length;
// 
//     TPZLadeKim LK2;
//     switch (choice) {
//         case(0):
//         {
//             TPZLadeKim::PlainConcrete(LK2);
//             std::ofstream outfiletxt("TPZLadeKim::PlainConcrete.txt");
// 
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 LK2.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//             }
//             break;
//         }
//         case(1):
//         {
//             TPZLadeKim::LooseSacrRiverSand(LK2);
//             std::ofstream outfiletxt("TPZLadeKim::LooseSacrRiverSand.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 LK2.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//         }
//         case(2):
//         {
//             TPZLadeKim::DenseSacrRiverSand(LK2);
//             std::ofstream outfiletxt("TPZLadeKim::DenseSacrRiverSand.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 LK2.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//         }
//         case(3):
//         {
//             TPZLadeKim::FineSilicaSand(LK2);
//             std::ofstream outfiletxt("TPZLadeKim::FineSilicaSand.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 LK2.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//         }
//         case(4):
//         {
//             cout << "\n poisson ";
//             REAL poisson;// = 0.18;
//             cin>>poisson;
//            
//             cout << "\n M ";
//             REAL M;//       = 361800.;
//             cin >> M;
//            
//             cout << "\nlambda";
//             REAL lambda;//  = 0.;
//             cin >> lambda;
//            
//             cout << "\n a";
//             REAL a;//       = 28.5;
//             cin >> a;
//            
//             cout << "\n m";
//             REAL m;//       = 1.113;
//             cin >> m;
//            
//             cout << "\n neta1";
//             REAL neta1;//   = 159800.;
//             cin >> neta1;
//            
//             cout << "\n ksi2";
//             REAL ksi2; //   = -2.92;
//             cin >> ksi2;
//            
//             cout << "\n mu ";
//             REAL mu;//     = 5.06;
//             cin >> mu;
//            
//             cout << "\n C";
//             REAL C;//       = 0.712E-12;
//             cin >> C;
//            
//             cout << "\n p ";
//             REAL p;//       = 3.8;
//             cin >> p;
//            
//             cout <<"\n h";
//             REAL h;//       = 1.990;
//             cin >> h;
//            
//             cout << "\n alpha";
//             REAL alpha;//   = 0.75;
//             cin >> alpha;
//            
//             cout << "\n pa";
//             REAL pa;//      = 14.7;
//             cin >> pa;
//            
//             REAL restol;
//             cout << "\n Tolerance";
//             cin >> restol;
//            
//             LK2.fResTol = restol;
//            
//             LK2.SetUp(poisson, M, lambda,
//                       a, m, neta1,
//                       ksi2, mu,
//                       C, p,
//                       h, alpha,
//                       pa);
//             std::ofstream outfiletxt("TPZLadeKim::YOURMODEL.txt");
//             for(int step=0;step<length;step++)
//             {
//                 cout << "\nstep "<< step;
//                 LK2.ApplyStrainComputeDep(strain, stress,Dep);
//                 outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//                 strain += deltastrain;
//                
//             }
//             break;
//         }
//         default:
//         {
//             cout << "Unknown Test Type. Exiting...";
//             break;
//         }
//     }
// 
// }
// 
// 
// void LKKoCompressionLoadTest()
// {
//     ofstream outfiletxt1("LKKoCompressionLoadTest.txt");
//     TPZTensor<REAL> stress, strain, deltastress, deltastrain;
//     TPZFNMatrix<6*6> Dep(6,6,0.);
//    
//     deltastress.XX()=-1.;
//     deltastress.YY()=-0.5;
//     deltastress.ZZ()=-0.5;
//     stress=deltastress;
//    
//    
//     TPZLadeKim LK2;
//     TPZLadeKim::FineSilicaSand(LK2);
//    
//     int length2 =20;
//     for(int step=0;step<length2;step++)
//     {
//         cout << "\nstep "<< step;
//         LK2.ApplyLoad(stress,strain);
//         LK2.ApplyStrainComputeDep(strain, stress, Dep);
//         outfiletxt1 << fabs(strain.I1()) << " " << fabs(stress.XX()/stress.ZZ()) << "\n";
//         stress += deltastress;
//         cout << "strain = "<<strain <<"\n";
//         cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
//        
//     }
//    
// }
// 
// 
// void LKLoadingTest()
// {
//    
//     ofstream outfiletxt1("FineSilica.txt");
//     ofstream outfiletxt2("PlainConcretee1.txt");
//     ofstream outfiletxt3("PlainConcretee2.txt");
//     ofstream outfiletxt4("PlainConcretee3.txt");
//     ofstream outfiletxt5("SDMcCormicRanchSand.txt");
//     TPZTensor<REAL> stress, strain, deltastress, deltastrain;
//     TPZFNMatrix<6*6> Dep(6,6,0.);
//    
//    
//    
//    
//     deltastress.XX() = -0.001;
//     deltastress.XY() = 0.;
//     deltastress.XZ() = 0.;
//     deltastress.YY() = -0.001;
//     deltastress.YZ() = 0.;
//     deltastress.ZZ() = -0.001;
//     stress = deltastress;
//    
//     TPZLadeKim LK;
//    
//     TPZLadeKim::FineSilicaSand(LK);
//     //    LK.ApplyLoad(stress,deltastrain);
//    
//     deltastress.XX() = -4.;
//     deltastress.XY() = 0.;
//     deltastress.XZ() = 0.;
//     deltastress.YY() = -4.;
//     deltastress.YZ() = 0.;
//     deltastress.ZZ() = -4.;
//     stress = deltastress;
//    
//    
//     //    deltastrain.XX() = -0.0001;
//     //    deltastrain.XY() = 0.;
//     //    deltastrain.XZ() = 0.;
//     //    deltastrain.YY() = -0.0001;
//     //    deltastrain.YZ() = 0.;
//     //    deltastrain.ZZ() = -0.0001;
//     //    strain=deltastrain;
//    
//     int length =72;
//     for(int step=0;step<length;step++)
//     {
//         cout << "\nstep "<< step;
//         if(step == 69)deltastress *=-1.;
//         LK.ApplyLoad(stress,strain);
//         //        if(step == 16)deltastrain *=-1.;
//         //        LK.ApplyStrainComputeDep(strain, stress,Dep);
//         if(step==0)
//         {
//            
//             //outfiletxt1 << 0. << " " << 0. << "\n";
//            
//         }
//        
//         else
//         {
//            
//             outfiletxt1 << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//            
//         }
//        
//         //        deltastress.Multiply(1.1, 1.);
//         stress +=  deltastress;
//         //        strain +=deltastrain;
//         cout << "strain = "<<strain <<"\n";
//         cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
//        
//     }
//    
//    
//     deltastrain.XX() = -0.0001;
//     deltastrain.XY() = 0.;
//     deltastrain.XZ() = 0.;
//     deltastrain.YY() = -0.;
//     deltastrain.YZ() = 0.;
//     deltastrain.ZZ() = -0.;
//     strain=deltastrain;
//     TPZLadeKim LK2;
//     TPZLadeKim::PlainConcrete(LK2);
//    
//     int length2 =30;
//     for(int step=0;step<length2;step++)
//     {
//         cout << "\nstep "<< step;
//         //    if(step == 10 || step == 16|| step == 40 || step==51 || step==80)deltastrain *=-1.;
//         LK2.ApplyStrainComputeDep(strain, stress,Dep);
//         outfiletxt2 << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//         outfiletxt3 << fabs(strain.YY()) << " " << fabs(stress.XX()) << "\n";
//         outfiletxt4 << fabs(strain.ZZ()) << " " << fabs(stress.XX()) << "\n";
//         strain += deltastrain;
//         cout << "strain = "<<strain <<"\n";
//         cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
//        
//     }
//    
//     TPZSandlerDimaggio SD;
//     TPZSandlerDimaggio::McCormicRanchSand(SD);
//    
//     deltastrain.XX() = -0.005;
//     deltastrain.XY() = 0.;
//     deltastrain.XZ() = 0.;
//     deltastrain.YY() = -0.;
//     deltastrain.YZ() = 0.;
//     deltastrain.ZZ() = -0.;
//     strain=deltastrain;
//    
//     int length3 =23;
//     for(int step=0;step<length3;step++)
//     {
//         cout << "\nstep "<< step;
//         if(step == 14 )deltastrain *=-1.;
//         SD.ApplyStrainComputeDep(strain, stress,Dep);
//         outfiletxt5 << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//         strain += deltastrain;
//         cout << "strain = "<<strain <<"\n";
//         cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
//        
//     }   
//    
//    
// }
// 
// 
// 
// //int main()
// //{
// //
// //   
// //   
// //   
// //    PressureCilinder();
// //    return 0;
// //}
// 
// 
// void DruckerIsotropicCompression()
// {
//     TPZDruckerPrager DP;
//     ofstream outfiletxt("DruckerPragerIsotropicCompression.txt");
//     cout << "\n Put the value of strain you want to add in each step of your loat test: ";
//     REAL straininput;
//     cin >> straininput;
//     TPZFNMatrix<6*6> Dep(6,6,0.);
//     TPZTensor<REAL> deltastrain,strain,stress,deltastress;
//     deltastrain.XX() = -straininput;
//     deltastrain.XY() = 0.;
//     deltastrain.XZ() = 0.;
//     deltastrain.YY() = -straininput;
//     deltastrain.YZ() = 0.;
//     deltastrain.ZZ() = -straininput;
//     strain=deltastrain;
//    
//     cout << "\n4 - Put the material parameters you want";
//    
//     cout << "\n Young modulus ";
//     REAL E;
//     cin >> E;
//    
//     cout << "\n Poisson";
//     REAL poisson;
//     cin >> poisson;
//    
//     int mcfit;
//     cout << "\n choose 0 for Iner Morh-Coulomb fit or 1 for outer Morh-Coulomb Fit ";
//     cin >> mcfit;
//    
//     if(mcfit!= 0 || mcfit!= 1)
//     {
//         cout << "\n wrong choice in Morh-Coulomb fit tipe 0 or 1";
//         return;
//     }
//    
//     REAL phi;
//     cout << "\n Type the internal frictional angle";
//     cin >> phi;
//    
//     REAL c;
//     cout << "\n Type the material coesion ";
//     cin >> c;
//    
//     REAL h;
//     cout << "\n Type the material hardening modulus ";
//     cin >> h;
//    
//     DP.fYC.SetUp(phi/180. *M_PI ,mcfit);
//     DP.fTFA.SetUp(c,h);
//     DP.fER.SetUp(E,poisson);
//    
//     int length =30;
//     for(int step=0;step<length;step++)
//     {
//         cout << "\nstep "<< step;   
//         DP.ApplyStrainComputeDep(strain, stress, Dep);
//         strain += deltastrain;
//         outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//     }
//    
// }
// 
// 
// void DruckerTest()
// {
//    
//     ofstream outfiletxt1("e1LK.txt");
//     ofstream outfiletxt2("e2LK.txt");
//     ofstream outfiletxt3("e3LK.txt");
//     ofstream outfiletxt4("VolLK.txt");    
//     TPZTensor<REAL> stress, strain, deltastress, deltastrain;
//     //
//     //    deltastress.XX() = -0.5;
//     //    deltastress.XY() = -0.001;
//     //    deltastress.XZ() = -0.001;
//     //    deltastress.YY() = -0.001;
//     //    deltastress.YZ() = -0.001;
//     //    deltastress.ZZ() = -0.001;
//    
//     deltastress.XX() = -0.1;
//     deltastress.XY() = -0.001;
//     deltastress.XZ() = -0.003;
//     deltastress.YY() = -0.15;
//     deltastress.YZ() = -0.0015;
//     deltastress.ZZ() = -0.17;
//    
//    
//    
//     //    deltastrain.XX() = -0.00001;
//     //    deltastrain.XY() = -0.0000001;
//     //    deltastrain.XZ() = -0.0000003;
//     //    deltastrain.YY() = -0.0000015;
//     //    deltastrain.YZ() = -0.00000015;
//     //    deltastrain.ZZ() = -0.0000017;
//     //    strain=deltastrain;
//    
//     //    TPZSandlerDimaggio SD;
//     //    TPZSandlerDimaggio::McCormicRanchSand(SD);
//    
//     TPZLadeKim LK;
//     //TPZLadeKim::FineSilicaSand(LK);
//     //TPZLadeKim::DenseSacrRiverSand(LK);
//     TPZLadeKim::PlainConcrete(LK);
//    
//     typedef TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> TPZDruckerPrager;
//     TPZDruckerPrager DP;
//     REAL pi = M_PI;
//     /*innerMCFit = 0*/
//     /*OuterMCFit = 1*/
//     DP.fYC.SetUp(/*phi=20*/ 20./180. * pi ,/*MCFit*/0);
//     REAL coesao = 9.2376;
//     DP.fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ coesao, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
//     DP.fER.SetUp(/*young*/ 20000., /*poisson*/ 0.2);
//    
//     //    LK.Print(cout);
//     LK.ApplyLoad(stress, strain);
//    
//     //deltastress.XX() = -0.8;
//     //deltastress.XY() = -0.001;
//     //deltastress.XZ() = -0.001;
//     //deltastress.YY() = -0.001;
//     //deltastress.YZ() = -0.001;
//     //deltastress.ZZ() = -0.001;
//     deltastress.XX() = -1000.;
//     deltastress.XY() = -0.001;
//     deltastress.XZ() = -0.003;
//     deltastress.YY() = -0.15;
//     deltastress.YZ() = -0.0015;
//     deltastress.ZZ() = -0.17;
//     stress = deltastress;
//    
//     int length =30;
//     for(int step=0;step<length;step++)
//     {
//         cout << "\nstep "<< step;   
//         LK.ApplyLoad(stress,strain);
//         REAL pa = stress.I1()/3.;
//         outfiletxt1 << strain.XX() << " " << fabs(stress.XX()) << "\n";
//         outfiletxt2 << strain.YY() << " " << fabs(stress.XX()) << "\n";
//         outfiletxt3 << strain.ZZ() << " " << fabs(stress.XX()) << "\n";
//         outfiletxt4 << strain.I1() << " " << fabs(stress.XX()) << "\n";
//         stress += deltastress;
//         cout << "strain = "<<strain <<"\n";
//         cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
//        
//     }
//    
// }
// 
// 
// int main_USANDO_OUTRO_main()
// {
// /*   
//     cout << "\nChoose Plasticity test:";
//     cout << "\n0 - Isotropic compression ";
//     cout << "\n1 - = ";
//     cout << "\n2 - Uniaxial traction ";
//     cout << "\n";
//     int choice;
//     cin >> choice;
//    
//     switch(choice)
//     {
//         case(0):
//             cout << "\n Choose the Plastic model tou need to run Isotropic compression: ";
//             cout << "\n0 - Lade - Kim ";
//             cout << "\n1 - Sandler Dimaggio ";
//             cout << "\n2 - Drucker Prager ";
//             cout << "\n";
//             int choice2;
//             cin >> choice2;
//         switch(choice2)
//         {
//             case(0):
//                 LKIsotropicCompression();
//                 break;
//             case(1):
//                 SandlerDimaggioIsotropicCompression();
//                 break;
//             case(2):
//                 DruckerIsotropicCompression();
//                 break;
//         }
//            
//            
//             break;
//         case(1):
//             cout << "NOT IMPLEMENTED YET";
//             //LadeKim_SimpleTest();
//             break;
//         case(2):
//             cout << "NOT IMPLEMENTED YET";
//         //    LadeKim_ReversalTest();
//             break;
//         default:
//             cout << "Unknown Test Type. Exiting...";
//     }
//     */
//    
//     PressureCilinder();
//    
//     return 0;
//    
// }