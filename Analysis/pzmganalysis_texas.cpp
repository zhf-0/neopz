// -*- c++ -*-
/* Generated by Together */

#include <stdlib.h>
#include "pzmganalysis.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pztransfer.h"
#include "pzadmchunk.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzskylmat.h"
#include "pzskylstrmatrix.h"

#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontStructMatrix.h"
#include "pzmgsolver.h"
#include "pzseqsolver.h"
#include "pzstepsolver.h"

#include "pzquad.h"

#include "pzelcq2d.h"
#include "pzelct2d.h"
#include "pzmat2dlin.h"

#include "pzonedref.h"

class TPZTransfer;
static ofstream deduce("deduce.txt");

TPZMGAnalysis::TPZMGAnalysis(TPZCompMesh *cmesh) : TPZAnalysis(cmesh) {
  fMeshes.Push(cmesh);
  fExact = 0;
//   fIterative = 0;
//   fPrecond = 0;
}

TPZMGAnalysis::~TPZMGAnalysis() {
  while (fMeshes.NElements()) delete fMeshes.Pop();
  while(fSolutions.NElements()) delete fSolutions.Pop();
  while(fSolvers.NElements()) delete fSolvers.Pop();
  while(fPrecondition.NElements()) delete fPrecondition.Pop();
//   delete fIterative;
//   delete fPrecond;
}


void TPZMGAnalysis::AppendMesh(TPZCompMesh * mesh){
  if(fMeshes.NElements() != fSolvers.NElements() || fMeshes.NElements() != fSolutions.NElements() ||
     fPrecondition.NElements() != fMeshes.NElements()) {
    cout << "TPZMGAnalysis::AppendMesh can only be called after solving the coarse mesh\n";
    return;
  }
  fMeshes.Push(mesh);
  TPZCompMesh *coarse = fCompMesh;

  fCompMesh = mesh;
  SetBlockNumber();
  int numeq = mesh->NEquations();
  //TPZSkylineStructMatrix skstr(mesh);
  TPZFrontStructMatrix<TPZFrontSym> sfstr(mesh);
  //  SetStructuralMatrix(skstr);
  SetStructuralMatrix(sfstr);
  int mg = 0;
  if(numeq > 15000 && mg) {
    int nmeshes = fSolvers.NElements();
    TPZTransfer *tr = new TPZTransfer;
    mesh->BuildTransferMatrix(*coarse,*tr);
    TPZFMatrix sol(mesh->NEquations(),1);
    tr->TransferSolution(fSolution,sol);
    fSolution = sol;
    mesh->LoadSolution(fSolution);
    TPZBlockDiagonalStructMatrix bdstr(mesh);
    TPZBlockDiagonal *bd = (TPZBlockDiagonal *) bdstr.Create();
    bdstr.AssembleBlockDiagonal(*bd);
    TPZStepSolver s1(bd);
    s1.SetDirect(ELU);
      //    TPZSkylMatrix *skyl = (TPZSkylMatrix *) skstr.CreateAssemble(fRhs);
    int nvar = mesh->MaterialVec()[0]->NStateVariables();

    TPZMatrixSolver *prec;
    prec = fPrecondition[nmeshes-1];
    if(!prec) prec = fSolvers[nmeshes-1];
    TPZMGSolver s2(tr,*prec,nvar);
    TPZSequenceSolver s3;
    s3.ShareMatrix(s2);
    s3.AppendSolver(s1);
//     s3.AppendSolver(s1);
    s3.AppendSolver(s2);
//     s3.AppendSolver(s1);
    s3.AppendSolver(s1);
//     fPrecond = (TPZMatrixSolver *) s3.Clone();
    TPZStepSolver s4;
    s4.ShareMatrix(s3);
    fPrecondition.Push((TPZMatrixSolver *)s3.Clone());
    s4.SetCG(200,s3,1.e-6,1);
//     fIterative = new TPZStepSolver(s4);
    SetSolver(s4);
    fSolvers.Push((TPZMatrixSolver *) s4.Clone());
  } else if(!mg && numeq > 15000) {
    TPZTransfer tr;
    //isso modifica a forma de apresenta��o da
    // matriz solu��o do coarse ou do mesh
    mesh->BuildTransferMatrix(*coarse,tr);
    TPZFMatrix sol(mesh->NEquations(),1);
    tr.TransferSolution(fSolution,sol);
    fSolution = sol;
    mesh->LoadSolution(fSolution);
    TPZBlockDiagonalStructMatrix bdstr(mesh);
    TPZBlockDiagonal *bd = (TPZBlockDiagonal *) bdstr.Create();
    bdstr.AssembleBlockDiagonal(*bd);
    TPZStepSolver s1(bd);
    s1.SetDirect(ELU);
      //    TPZSkylMatrix *skyl = (TPZSkylMatrix *) skstr.CreateAssemble(fRhs);
    TPZStepSolver s4;
    fPrecondition.Push((TPZMatrixSolver *)s1.Clone());
    s4.SetCG(200,s1,1.e-6,1);
//     fIterative = new TPZStepSolver(s4);
    SetSolver(s4);
    fSolvers.Push((TPZMatrixSolver *) s4.Clone());
  } else {
    TPZStepSolver s4;
    s4.SetDirect(ECholesky);
    fSolvers.Push((TPZMatrixSolver *) s4.Clone());
    fPrecondition.Push(0);
    SetSolver(s4);
  }
}

TPZCompMesh *TPZMGAnalysis::PopMesh() {
  if(fMeshes.NElements() == 1) {
    cout << "TPZMGAnalysis cannot delete the root mesh, sorry\n";
    return 0;
  }
//   if(!fPrecond) {
//     cout << "TPZMGAnalysis::PopMesh I don't understand\n";
//   }
  if(fSolutions.NElements() == fMeshes.NElements()) delete fSolutions.Pop();
//   delete fPrecond;
//   delete fIterative;
  delete fSolvers.Pop();
  delete fPrecondition.Pop();
//   TPZMatrixSolver *Precond = dynamic_cast<TPZMatrixSolver *>(fSolvers[fSolvers.NElements()-1]->Clone());
//   if(!Precond) {
//     cout << "TPZMGAnalysis::PopMesh dynamic cast failed?";
//   }
  SetSolver(*fSolvers[fSolvers.NElements()-1]);
  fCompMesh = fMeshes[fMeshes.NElements()-2];
  fSolution = *fSolutions[fSolutions.NElements()-1];
  return fMeshes.Pop();
}


int TPZMGAnalysis::main() {

  TPZCompMesh *cmesh = TPZCompElQ2d::CreateMesh();
  cmesh->CleanUpUnconnectedNodes();
  TPZFMatrix sol(cmesh->NEquations(),1);
  ofstream out("output.txt");
  int row = sol.Rows();
  int r;
   for(r=0; r<row; r++) {
//     //    sol(r,0) = rand()/(RAND_MAX+1.);
     sol(r,0) = 1.;
    
   }
  TPZMGAnalysis mgan(cmesh);
  TPZSkylineStructMatrix strskyl(cmesh);
  mgan.SetStructuralMatrix(strskyl);
  TPZStepSolver direct;
  direct.SetDirect(ELDLt);
  mgan.SetSolver(direct);
  mgan.Run();
  TPZGeoMesh *gmesh = cmesh->Reference();
  int nel = gmesh->ElementVec().NElements();
  int el;
  TPZVec<TPZGeoEl *> sub;
  for(el=0; el<nel; el++) {
    TPZGeoEl *gel = gmesh->ElementVec()[el];
    if(!gel) continue;
    //     if(!gel->HasSubElement(0)) {
       gel->Divide(sub);
       //       break;
       //     }
  }
  gmesh->ResetReference();
  TPZCompMesh *cmesh2 = new TPZCompMesh(gmesh);
  //  cmesh2.LoadReferences();
  TPZAdmChunkVector<TPZMaterial *> &matvec2 = cmesh2->MaterialVec();
  TPZAdmChunkVector<TPZMaterial *> &matvec = cmesh->MaterialVec();
  int nmat = matvec.NElements();
  TPZMaterial *mat;
  int im;
  for(im=0; im<nmat; im++) {
    mat = matvec[im];
    if(!mat) continue;
    mat->Clone(matvec2);
  }
//   TPZFMatrix xk(1,1,1.),xc(1,1,1.),xf(1,1,1.);
//    //   xk(0,1) = xk(1,0) = xc(0,1) = xc(1,0) = 0.;
//   mat->SetMaterial(xk,xc,xf);
//   TPZFMatrix val1(1,1,0.),val2(1,1,0.);
//   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
//   cmesh2->InsertMaterialObject(mat);
//   cmesh2->InsertMaterialObject(bc);

  //  gmesh->Print();
  cmesh2->AutoBuild();
  cmesh2->AdjustBoundaryElements();
  cmesh2->CleanUpUnconnectedNodes();
  mgan.AppendMesh(cmesh2);

  mgan.Run();
  // TPZMGAnalysis an(cmesh2);
  //  TPZTransfer *trf = TPZMGAnalysis::BuildTransferMatrix(cmesh2,cmesh);
  cmesh->LoadSolution(sol);
  TPZTransfer trf;
  cmesh2->BuildTransferMatrix(*cmesh,trf);
  trf.Print("Transfer Matrix",out);
  TPZFMatrix sol2(cmesh2->NEquations(),1,0.);
  trf.TransferSolution(sol,sol2);
  cmesh2->LoadSolution(sol2);
  gmesh->Print(out);
  cmesh->Print(out);
  cmesh2->Print(out);
  cmesh->Solution().Print("Coarse mesh solution",out);
  cmesh2->Solution().Print("Fine mesh solution",out);

  TPZVec<REAL> ervec,truervec;
  TPZMGAnalysis::MeshError(cmesh2,cmesh,ervec,mgan.fExact,truervec);
  int i;
  cout << "TPZMGAnalysis the error between both meshes is equal to \n";
  for(i=0; i<ervec.NElements(); i++) cout << ervec[i] << ' ';
  cout << endl;
  return 1;
}


void TPZMGAnalysis::MeshError(TPZCompMesh *fine, TPZCompMesh *coarse, TPZVec<REAL> &ervec,
			      void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix &deriv),TPZVec<REAL> &truervec){
  coarse->Reference()->ResetReference();
  coarse->LoadReferences();
  int dim = fine->MaterialVec()[0]->Dimension();
  ervec.Resize(coarse->NElements());
  if(f) {
    truervec.Resize(coarse->NElements());
    truervec.Fill(0.,0);
  }
  ervec.Fill(0.,0);
  TPZCompEl *cel;
  TPZAdmChunkVector<TPZCompEl *> &elementvec = fine->ElementVec();
  int numel = elementvec.NElements();
  int el;
  for(el=0; el<numel; el++) {
    cel = elementvec[el];
    if(!cel || !cel->IsInterpolated()) continue;
    TPZInterpolatedElement *cint = (TPZInterpolatedElement *) cel;
    int ncon = cint->NConnects();
    TPZGeoElSide gelside(cint->Reference(),ncon-1);
    if(gelside.Dimension() != dim) continue;
    TPZGeoElSide gellarge(gelside);
    while(!gellarge.Reference().Exists() && gellarge.Father2().Exists()) gellarge = gellarge.Father2();
    if(!gellarge.Reference().Exists()) {
      cout << "TPZMGAnalsysis::BuildTransferMatrix element " << el << " found no corresponding element\n";
      continue;
    }
    TPZCompElSide cellargeside = gellarge.Reference();
    TPZCompEl *cellarge = cellargeside.Element();
    TPZInterpolatedElement *cintlarge = (TPZInterpolatedElement *) cellarge;
    TPZTransform transform(gelside.Dimension(),gellarge.Dimension());
    gelside.SideTransform3(gellarge,transform);
    int index = cellarge->Index();
    REAL truerror = 0.;
    ervec[index] += ElementError(cint,cintlarge,transform,f,truerror);
    if(f) truervec[index]  += truerror;
  }
}


REAL TPZMGAnalysis::ElementError(TPZInterpolatedElement *fine, TPZInterpolatedElement *coarse, TPZTransform &tr,
void (*f)(TPZVec<REAL> &loc, TPZVec<REAL> &val, TPZFMatrix &deriv),REAL &truerror){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t
  int locnod = fine->NConnects();
  int cornod = coarse->NConnects();
  int locmatsize = fine->NShapeF();
  int cormatsize = coarse->NShapeF();
  REAL error = 0.;
  truerror = 0.;

  REAL loclocmatstore[500] = {0.},loccormatstore[500] = {0.};
  TPZFMatrix loclocmat(locmatsize,locmatsize,loclocmatstore,500);
  TPZFMatrix loccormat(locmatsize,cormatsize,loccormatstore,500);

  TPZIntPoints &intrule = fine->GetIntegrationRule();
  int dimension = fine->Dimension();
  int numdof = fine->Material()->NStateVariables();
  TPZBlock &locblock = fine->Mesh()->Block();
  TPZFMatrix &locsolmesh = fine->Mesh()->Solution();

  TPZBlock &corblock = coarse->Mesh()->Block();
  TPZFMatrix &corsolmesh = coarse->Mesh()->Solution();

  TPZVec<REAL> locsol(numdof);
  TPZFMatrix locdsol(dimension,numdof);

  TPZVec<REAL> corsol(numdof);
  TPZFMatrix cordsol(dimension,numdof);

  TPZManVector<int> prevorder(dimension),
    order(dimension);
  intrule.GetOrder(prevorder);

  TPZManVector<int> interpolation(dimension);
  fine->GetInterpolationOrder(interpolation);

// compute the interpolation order of the shapefunctions squared
  int dim;
  int maxorder = interpolation[0];
  for(dim=0; dim<interpolation.NElements(); dim++) {
    maxorder = interpolation[dim] < maxorder ? maxorder : interpolation[dim];
  }
  for(dim=0; dim<dimension; dim++) {
	  order[dim] = 20;
  }
  intrule.SetOrder(order);


  REAL locphistore[50]={0.},locdphistore[150]={0.};
  TPZFMatrix locphi(locmatsize,1,locphistore,50);
  TPZFMatrix locdphi(dimension,locmatsize,locdphistore,150),locdphix(dimension,locmatsize);
  // derivative of the shape function
  // in the master domain

  REAL corphistore[50]={0.},cordphistore[150]={0.};
  TPZFMatrix corphi(cormatsize,1,corphistore,50);
  TPZFMatrix cordphi(dimension,cormatsize,cordphistore,150), cordphix(dimension,cormatsize);
   				// derivative of the shape function
  // in the master domain

  REAL jacobianstore[9],
    axesstore[9];
  TPZManVector<REAL> int_point(dimension),
    coarse_int_point(dimension);
  TPZFMatrix jacfine(dimension,dimension,jacobianstore,9),jacinvfine(dimension,dimension);
  TPZFMatrix axesfine(3,3,axesstore,9);
  TPZManVector<REAL> xfine(3);
  TPZFMatrix jaccoarse(dimension,dimension),jacinvcoarse(dimension,dimension);
  TPZFMatrix axescoarse(3,3), axesinner(3,3);
  TPZManVector<REAL> xcoarse(3);

  REAL jacdetcoarse;
  int numintpoints = intrule.NPoints();
  REAL weight;
  int i,j,k;

  TPZVec<REAL> truesol(numdof);
  TPZFMatrix truedsol(dimension,numdof);
  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule.Point(int_ind,int_point,weight);
    REAL jacdetfine;
    fine->Reference()->Jacobian( int_point, jacfine , axesfine, jacdetfine, jacinvfine);
    fine->Reference()->X(int_point, xfine);
    if(f) f(xfine,truesol,truedsol);
    fine->Shape(int_point,locphi,locdphi);
    tr.Apply(int_point,coarse_int_point);
    coarse->Shape(coarse_int_point,corphi,cordphi);
    coarse->Reference()->Jacobian( coarse_int_point,jaccoarse,axescoarse, jacdetcoarse, jacinvcoarse);
    coarse->Reference()->X(coarse_int_point,xcoarse);
    REAL dist = (xfine[0]-xcoarse[0])*(xfine[0]-xcoarse[0])+(xfine[1]-xcoarse[1])*(xfine[1]-xcoarse[1])+(xfine[2]-xcoarse[2])*(xfine[2]-xcoarse[2]);
    if(dist > 1.e-6) cout << "TPZMGAnalysis::ElementError transformation between fine and coarse is wrong\n";
    for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
	axesinner(i,j) = 0.;
	for(k=0; k<3; k++) axesinner(i,j) += axesfine(i,k)*axescoarse(j,k);
      }
    }
    if(fabs(axesinner(0,0)-1.) > 1.e-6 || fabs(axesinner(1,1)-1.) > 1.e-6 || fabs(axesinner(0,1)) > 1.e-6 || fabs(axesinner(1,0)) > 1.e-6) {
      cout << "TPZMGAnalysis axesinner is not identify?\n";
    }
    weight *= fabs(jacdetfine);

    locdphix.Zero();

    int ieq,d;
    switch(dim) {
    case 0:
      //dphix.Redim(1,1);
      //dphix(0,0) = dphi(0,0);
      break;
    case 1:
      for(d=0; d<dimension; d++) {
	for(ieq=0; ieq<locmatsize; ieq++) locdphix(d,ieq) = locdphi(d,ieq)*(1./jacdetfine);
	for(ieq=0; ieq<cormatsize; ieq++) cordphix(d,ieq) = cordphi(d,ieq)*(axesinner(0,0)/jacdetcoarse);
      }
      break;
    case 2:
	    for(ieq = 0; ieq < locmatsize; ieq++) {
	        locdphix(0,ieq) = jacinvfine(0,0)*locdphi(0,ieq) + jacinvfine(1,0)*locdphi(1,ieq);
	        locdphix(1,ieq) = jacinvfine(0,1)*locdphi(0,ieq) + jacinvfine(1,1)*locdphi(1,ieq);
		REAL tmp[2];
		tmp[0] = locdphix(0,ieq)*axesfine(0,0)+locdphix(1,ieq)*axesfine(1,0);
		tmp[1] = locdphix(0,ieq)*axesfine(0,1)+locdphix(1,ieq)*axesfine(1,1);
		locdphix(0,ieq) = tmp[0];
		locdphix(1,ieq) = tmp[1];
	    }
	    for(ieq = 0; ieq < cormatsize; ieq++) {
	        cordphix(0,ieq) = jacinvcoarse(0,0)*cordphi(0,ieq) + jacinvcoarse(1,0)*cordphi(1,ieq);
	        cordphix(1,ieq) = jacinvcoarse(0,1)*cordphi(0,ieq) + jacinvcoarse(1,1)*cordphi(1,ieq);
		REAL tmp[2];
		tmp[0] = cordphix(0,ieq)*axescoarse(0,0)+cordphix(1,ieq)*axescoarse(1,0);
		tmp[1] = cordphix(0,ieq)*axescoarse(0,1)+cordphix(1,ieq)*axescoarse(1,1);
		cordphix(0,ieq) = tmp[0];
		cordphix(1,ieq) = tmp[1];
	    }
      break;
    case 3:
	    for(ieq = 0; ieq < locmatsize; ieq++) {
	        locdphix(0,ieq) = jacinvfine(0,0)*locdphi(0,ieq) + jacinvfine(0,1)*locdphi(1,ieq) + jacinvfine(0,2)*locdphi(2,ieq);
	        locdphix(1,ieq) = jacinvfine(1,0)*locdphi(0,ieq) + jacinvfine(1,1)*locdphi(1,ieq) + jacinvfine(1,2)*locdphi(2,ieq);
		locdphix(2,ieq) = jacinvfine(2,0)*locdphi(0,ieq) + jacinvfine(2,1)*locdphi(1,ieq) + jacinvfine(2,2)*locdphi(2,ieq);
	    }
	    for(ieq = 0; ieq < cormatsize; ieq++) {
	        cordphix(0,ieq) = jacinvcoarse(0,0)*cordphi(0,ieq) + jacinvcoarse(0,1)*cordphi(1,ieq) + jacinvcoarse(0,2)*cordphi(2,ieq);
	        cordphix(1,ieq) = jacinvcoarse(1,0)*cordphi(0,ieq) + jacinvcoarse(1,1)*cordphi(1,ieq) + jacinvcoarse(1,2)*cordphi(2,ieq);
		cordphix(2,ieq) = jacinvcoarse(2,0)*cordphi(0,ieq) + jacinvcoarse(2,1)*cordphi(1,ieq) + jacinvcoarse(2,2)*cordphi(2,ieq);
		REAL tmp[3];
		tmp[0] = cordphix(0,ieq)*axesinner(0,0)+cordphix(1,ieq)*axesinner(0,1)+cordphix(2,ieq)*axesinner(0,2);
		tmp[1] = cordphix(0,ieq)*axesinner(1,0)+cordphix(1,ieq)*axesinner(1,1)+cordphix(2,ieq)*axesinner(1,2);
		tmp[2] = cordphix(0,ieq)*axesinner(2,0)+cordphix(1,ieq)*axesinner(2,1)+cordphix(2,ieq)*axesinner(2,2);
		cordphix(0,ieq) = tmp[0];
		cordphix(1,ieq) = tmp[1];
		cordphix(2,ieq) = tmp[2];
	    }
      break;
    default:
      PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      PZError.flush();
    }
    int iv=0;
    locsol.Fill(0.);
    locdsol.Zero();
    iv=0;
	int in;
    for(in=0; in<locnod; in++) {
        TPZConnect *df = &fine->Connect(in);
        int dfseq = df->SequenceNumber();
        int dfvar = locblock.Size(dfseq);
	int pos = locblock.Position(dfseq);

        for(int jn=0; jn<dfvar; jn++) {
          locsol[iv%numdof] += locphi(iv/numdof,0)*locsolmesh(pos+jn,0);
          for(d=0; d<dim; d++)
             locdsol(d,iv%numdof) += locdphix(d,iv/numdof)*locsolmesh(pos+jn,0);
	  iv++;
        }
    }
    corsol.Fill(0.);
    cordsol.Zero();
    iv=0;
    for(in=0; in<cornod; in++) {
        TPZConnect *df = &coarse->Connect(in);
        int dfseq = df->SequenceNumber();
        int dfvar = corblock.Size(dfseq);
	int pos = corblock.Position(dfseq);
        for(int jn=0; jn<dfvar; jn++) {
          corsol[iv%numdof] += corphi(iv/numdof,0)*corsolmesh(pos+jn,0);
          for(d=0; d<dim; d++)
             cordsol(d,iv%numdof) += cordphix(d,iv/numdof)*corsolmesh(pos+jn,0);
	  iv++;
        }
    }
    int jn;
    for(jn=0; jn<numdof; jn++) {
//       error += (locsol[jn]-corsol[jn])*(locsol[jn]-corsol[jn])*weight;
//       if(f) truerror += (corsol[jn]-truesol[jn])*(corsol[jn]-truesol[jn])*weight;
      for(d=0; d<dim; d++) {
	error += (locdsol(d,jn)-cordsol(d,jn))*(locdsol(d,jn)-cordsol(d,jn))*weight;
	if(f) truerror += (cordsol(d,jn)-truedsol(d,jn))*(cordsol(d,jn)-truedsol(d,jn))*weight;
      }
    }
  }
  intrule.SetOrder(prevorder);
  return error;
}

void TPZMGAnalysis::Solve() {
  if(fMeshes.NElements() == 1) {
    TPZAnalysis::Solve();
    if(fSolvers.NElements() == 0) {
      fSolvers.Push((TPZMatrixSolver *) fSolver->Clone());
    }
    if(fPrecondition.NElements() == 0) {
      fPrecondition.Push(0);
    }
    if(fSolutions.NElements() == 0) {
      fSolutions.Push(new TPZFMatrix(fSolution));
    } else {
      int nsol = fSolutions.NElements();
      *(fSolutions[nsol-1]) = fSolution;
    }
    return;
  }
  int numeq = fCompMesh->NEquations();
  if(fRhs.Rows() != numeq ) return;
  int nsolvers = fSolvers.NElements();
  
  TPZFMatrix residual(fRhs);
  TPZFMatrix delu(numeq,1);
  TPZMatrixSolver *solve = dynamic_cast<TPZMatrixSolver *> (fSolvers[nsolvers-1]);
  if(fSolution.Rows() != numeq) {
    fSolution.Redim(numeq,1);
  } else {
    solve->Matrix()->Residual(fSolution,fRhs,residual);
  }
  
  REAL normrhs = Norm(fRhs);
  REAL normres  = Norm(residual);
  if(normrhs*1.e-6 >= normres) {
    cout << "TPZMGAnalysis::Solve no need for iterations normrhs = " << normrhs << " normres = " << normres << endl;
    if(fSolutions.NElements() < fMeshes.NElements()) {
      fSolutions.Push(new TPZFMatrix(fSolution));
    } else {
      int nsol = fSolutions.NElements();
      *(fSolutions[nsol-1]) = fSolution;
    }
    return ;
  }
//   REAL tol = 1.e-6*normrhs/normres;
//   if(numeq > 1500) {
//     fIterative->SetCG(200,*fPrecond,tol,0);
//   } else {
//     fIterative->SetDirect(ELDLt);
//   }
  TPZStepSolver *stepsolve = dynamic_cast<TPZStepSolver *> (solve);
  if(stepsolve) stepsolve->SetTolerance(1.e-6*normrhs/normres);
  cout << "TPZMGAnalysis::Run res : " << Norm(residual) << " neq " << numeq << endl;
  solve->Solve(residual, delu);
  fSolution += delu;
	
  fCompMesh->LoadSolution(fSolution);
  if(fSolutions.NElements() < fMeshes.NElements()) {
    fSolutions.Push(new TPZFMatrix(fSolution));
  } else {
    int nsol = fSolutions.NElements();
    *(fSolutions[nsol-1]) = fSolution;
  }
}

TPZCompMesh  *TPZMGAnalysis::UniformlyRefineMesh(TPZCompMesh *mesh) {

  TPZGeoMesh *gmesh = mesh->Reference();
  if(!gmesh) {
    cout << "TPZMGAnalysis::UniformlyRefineMesh encountered no geometric mesh\n";
    return 0;
  }
  gmesh->ResetReference();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  int nmat = mesh->MaterialVec().NElements();
  int m;
  for(m=0; m<nmat; m++) {
    TPZMaterial *mat = mesh->MaterialVec()[m];
    if(!mat) continue;
    mat->Clone(cmesh->MaterialVec());
  }
  TPZAdmChunkVector<TPZCompEl *> &elementvec = mesh->ElementVec();
  int el,nelem = elementvec.NElements();
  for(el=0; el<nelem; el++) {
    TPZCompEl *cel = elementvec[el];
    if(!cel) continue;
    if(!cel->IsInterpolated()) {
      cout << "TPZMGAnalysis::UniformlyRefineMesh encountered a non interpolated element\n";
      continue;
    }
    TPZInterpolatedElement *cint = (TPZInterpolatedElement *) cel;
    int ncon = cint->NConnects();
    int porder = cint->PreferredSideOrder(ncon-1);

    TPZGeoEl *gel = cint->Reference();
    if(!gel) {
      cout << "TPZMGAnalysis::UniformlyRefineMesh encountered an element without geometric reference\n";
      continue;
    }
    TPZVec<TPZGeoEl *> sub;
    gel->Divide(sub);
    int nsub = sub.NElements();
    int isub;
    int celindex;
    for(isub=0; isub<nsub; isub++) {
      TPZInterpolatedElement *csint;
      csint = (TPZInterpolatedElement *) sub[isub]->CreateCompEl(*cmesh,celindex);
      csint->PRefine(porder+1);
    }
  }
  return cmesh;
}

void TPZMGAnalysis::ComputeError(TPZVec<REAL> &error) {
  int nmesh = fMeshes.NElements();
  int nsol = fSolutions.NElements();
  if(nsol != nmesh || nsol < 2) {
    cout << "TPZMGAnalysis cannot compute the errors\n";
    return;
  }
  fMeshes[nsol-2]->LoadSolution(*fSolutions[nsol-2]);
  fMeshes[nsol-1]->LoadSolution(*fSolutions[nsol-1]);
  TPZVec<REAL> truerror;
  MeshError(fMeshes[nsol-1],fMeshes[nsol-2],error,0,truerror);
  //  MeshError(fMeshes[nsol-1],fMeshes[nsol-2],error);
}

void TPZMGAnalysis::Sort(TPZVec<REAL> &vec, TPZVec<int> &perm) {

  int i,j;
  int imin = 0;
  int imax = vec.NElements();
  for(i=imin; i<imax; i++) {
    for(j=i+1; j<imax; j++) {
      if(vec[perm[i]] < vec[perm[j]]) {
	int kp = perm[i];
	perm[i] = perm[j];
	perm[j] = kp;
      }
    }
  }
}

TPZCompMesh *TPZMGAnalysis::RefinementPattern(TPZCompMesh *fine, TPZCompMesh *coarse,
					      REAL &totalerror, REAL &totaltruerror, TPZVec<REAL> &effect) {

  TPZVec<REAL> ervec,truervec;
  TPZStack<TPZGeoEl *> gelstack;
  TPZStack<int> porders;

  MeshError(fine,coarse,ervec,fExact,truervec);
  coarse->SetElementSolution(0,ervec);
  if(fExact) coarse->SetElementSolution(1,truervec);
  int nstate;
  TPZMaterial *mat = fine->MaterialVec()[0];
  nstate = mat->NStateVariables();

  TPZGeoMesh *gmesh = fine->Reference();
  gmesh->ResetReference();
  fine->LoadReferences();
  int nelem = ervec.NElements();
  TPZVec<int> perm(nelem);
  int i;
  for(i=0; i<nelem; i++) perm[i] = i;
  Sort(ervec,perm);
  totalerror = 0.;
  totaltruerror = 0.;
  for(i=0; i<nelem; i++) totalerror += ervec[i];
  if(fExact) {
    for(i=0; i<nelem; i++) totaltruerror += truervec[i];
  }
  effect.Resize(truervec.NElements());
  effect.Fill(0.);
  if(fExact) {
    for(i=0; i<nelem; i++) {
      if(truervec[i] >= 1.e-4*totaltruerror) effect[i] = ervec[i]/truervec[i];
    }
    coarse->SetElementSolution(2,effect);
  }
  REAL ninetyfivepercent = 0.65*totalerror;
  REAL r = 0.;
  i=0;
//   for(i=0; i<nelem; i++) cout << ervec[perm[i]] << " ";
//   cout << endl;
  int ilast=0;
  while(r<ninetyfivepercent && ilast<nelem) {
    r+= ervec[perm[ilast++]];
  }
  REAL lasterror = ervec[perm[ilast-1]];
  while(ilast < nelem && fabs(ervec[perm[ilast]]-lasterror)/lasterror < 1.e-3) ilast++;
//   cout << "elements which are acountable for 95% of the error\n";
//   for(i=0; i<ilast; i++) {
//     cout << "element " << perm[i] << " error " << ervec[perm[i]] << ' ';
//     int j;
//     TPZGeoEl *ref = coarse->ElementVec()[perm[i]]->Reference();
//     cout << " elid " << ref->Id() << ' ';
//     for(j=0; j<ref->NCornerNodes(); j++) cout << ref->NodeIndex(j) << ' ';
//     cout << endl;
//   }
  for(i=ilast; i<nelem; i++) {
    TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (coarse->ElementVec()[perm[i]]);
    //  cint->Print();
    if(!cint) {
      if(coarse->ElementVec()[perm[i]]) {
	cout << "TPZMGAnalysis RefinementPattern only for interpolated meshes - I don't understand\n";
	cout << " i = " << i << " perm[i] = " << perm[i] << " nelem = " << nelem << endl;
      }
      continue;
    }
    int pord = cint->PreferredSideOrder(cint->NConnects()-1);
    gelstack.Push(cint->Reference());
    porders.Push(pord);
  }
  TPZOneDRef f(nstate);
  for(i=0; i<ilast; i++) {
    TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *>( coarse->ElementVec()[perm[i]]);
//     cout << "analysing element " << perm[i] << endl;
    AnalyseElement(f,cint,gelstack,porders);
  }
  return CreateCompMesh(coarse,gelstack,porders);  
}

void TPZMGAnalysis::AnalyseElement(TPZOneDRef &f, TPZInterpolatedElement *cint, TPZStack<TPZGeoEl *> &subels, TPZStack<int> &porders) {

  TPZGeoEl *gel = cint->Reference();
  int ncon = cint->NConnects();
  int intorder = cint->SideOrder(ncon-1);
  deduce << "Internal order = " << intorder << endl;
  TPZGeoMesh *gmesh = gel->Mesh();
  int ncorners = gel->NCornerNodes();
  TPZVec<int> cornerids(ncorners);
  TPZVec<int> localporders(ncorners);
  int id;
  for(id=0; id<ncorners; id++) cornerids[id] = gel->NodeIndex(id);
  int n1dsides = 0;
  int nsides = gel->NSides();
  int side;
  for(side=0; side<nsides; side++) if(gel->SideDimension(side) == 1) n1dsides++;
  TPZVec<TPZRefPattern> refpattern(n1dsides);
  n1dsides = 0;
  for(side=0; side<nsides; side++) {
    int sidedimension = gel->SideDimension(side);
    if(sidedimension != 1) continue;
//    int level = gel->Level();
    TPZStack<TPZCompElSide> subelsides;
    TPZGeoElSide gelside(gel,side);
//     if(gelside.Neighbour().Exists()) gelside = gelside.Neighbour();
//     else {
//       cout << "TPZAnalyseElement coarse element without neighbour\n";
//       continue;
//     }
    gelside.HigherLevelCompElementList2(subelsides,1,1);
    if(subelsides.NElements() != 2) {
      cout << "A one dimensional side with more than one subelement\n";
      gelside.HigherLevelCompElementList2(subelsides,1,1);
      continue;
    }
    TPZGeoElSide gels1 = subelsides[0].Reference();
    TPZGeoElSide gels2 = subelsides[1].Reference();
    if(gels1.SideNodeIndex(1) != gels2.SideNodeIndex(0)) {
      TPZCompElSide temp = subelsides[0];
      subelsides[0]=subelsides[1];
      subelsides[1] = temp;
      gels1 = subelsides[0].Reference();
      gels2 = subelsides[1].Reference();
    }
    if(gels1.SideNodeIndex(1) != gels2.SideNodeIndex(0)) {
      cout << "Unexpected situation\n";
      cout << gels1.SideNodeIndex(0) << ' ' << gels1.SideNodeIndex(1) << ' ' << gels2.SideNodeIndex(0) << ' ' << gels2.SideNodeIndex(1) <<endl;
      gels1.Element()->Print();
      gels2.Element()->Print();
      continue;
    }
    TPZVec<int> id(3);
    id[0] = gels1.SideNodeIndex(0);
    id[1] = gels1.SideNodeIndex(1);
    id[2] = gels2.SideNodeIndex(1);
//     cout << "side = " << side << " indexes " << id[0] << ' ' << id[1] << ' ' << id[2] << endl;
    REAL del[3];
    int i;
    for(i=0; i<3; i++) del[i] = gmesh->NodeVec()[id[1]].Coord(i)-gmesh->NodeVec()[id[0]].Coord(i);

    REAL delx = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);
    TPZCompMesh *cmesh = subelsides[0].Element()->Mesh();
    TPZConnect *connects[5];
    TPZInterpolatedElement *c1 = (TPZInterpolatedElement *) subelsides[0].Element();
    int s1 = subelsides[0].Side();
    TPZInterpolatedElement *c2 = (TPZInterpolatedElement *) subelsides[1].Element();
    int s2 = subelsides[1].Side();
    connects[0] = c1->SideConnect(0,s1);
    connects[2] = c1->SideConnect(1,s1);
    connects[1] = c1->SideConnect(2,s1);
    connects[4] = c2->SideConnect(1,s2);
    connects[3] = c2->SideConnect(2,s2);
    int dof = 0;
    for(i=0; i<5; i++) dof += connects[i]->NDof(*cmesh);
    TPZFMatrix U(dof,1);
    int p1 = c1->SideOrder(s1);
    int p2 = c2->SideOrder(s2);
    deduce << "side = " << side << " order = " << p1 << ' ' << p2 << endl;
    int c,counter=0;
    for(c=0; c<5; c++) {
      int nd = connects[c]->NDof(*cmesh);
      int seq = connects[c]->SequenceNumber();
      int d;
      for(d=0; d<nd; d++) U(counter++,0) = (*cmesh).Block()(seq,0,d,0);
    }
    int hp1, hp2;
    REAL hperror;
    REAL error = f.BestPattern(U,id,p1,p2,hp1, hp2, hperror,delx);
    TPZRefPattern optimal = {id[0],id[1],id[2],p1,p2,hp1,hp2,hperror,error};
    refpattern[n1dsides] = optimal;
    n1dsides++;
  }
  refpattern.Resize(n1dsides);
  DeduceRefPattern(refpattern,cornerids,localporders,intorder);
  if(localporders[1] == -1) {
    subels.Push(gel);
    porders.Push(localporders[0]);
//     cout << "prefine order " << localporders[0] << endl;
    return;
  }
  TPZStack<TPZGeoElSide> gelsides;
  int dimension = gel->Dimension();

  gel->GetSubElements2(nsides-1,gelsides,dimension);
  for(id=0; id<ncorners; id++) {
    subels.Push(gelsides[id].Element());
    porders.Push(localporders[id]);
  }
  // hook for triangular, tetrahedral, prism and pyramid elements
  if(gelsides.NElements() > ncorners) {
    // compute the maximum interpolation order of all "corner" elements"
    int maxorder = 0;
    for(id=0; id<ncorners; id++) maxorder = (maxorder<localporders[id]) ? localporders[id] : maxorder;
    // Assign the maximum order to all remaining elements
    for(id=ncorners; id<gelsides.NElements(); id++) {
      subels.Push(gelsides[id].Element());
      porders.Push(maxorder);
    }
  }
}

void TPZMGAnalysis::HeapSort(TPZVec<REAL> &sol, TPZVec<int> &perm){

  int nelem = perm.NElements();
  int i,j;
  for(i=0; i<nelem; i++) perm[i] = i;

  if(nelem == 1) return;
  int l, ir,ind;
  REAL q;
  l= nelem/2;
  ir = nelem-1;
  while(l>0 && ir>0) {
    if(l> 0) {
      l--;
      ind = perm[l];
      q=sol[ind];
    } else {
      ind = perm[ir];
      q = sol[ind];
      perm[ir] = perm[0];
      ir--;
    }
    i=l;
    j=l+l+1;
    while(j<=ir) {
      if(j<ir && sol[perm[j]] < sol[perm[j+1]]) j++;
      if(q < sol[perm[j]]) {
	perm[i] = perm[j];
	i=j;
	j= i+i+1;
      } else {
	break;
      }
    }
    perm[i] = ind;
  }
}


void TPZMGAnalysis::DeduceRefPattern(TPZVec<TPZRefPattern> &refpat, TPZVec<int> &cornerids, TPZVec<int> &porders, int originalp) {

  // Eliminate the refinement pattern suggestion if the error is smaller than 10% of the total error
  int nref = refpat.NElements();
  REAL totalerror = 0.;
  int ir;
  for(ir=0; ir<nref; ir++) {
    totalerror += refpat[ir].fError;
  }

  int ncorners = cornerids.NElements();
  //Print the incoming refpattern to the log file
  for(ir=0; ir<ncorners; ir++) deduce << cornerids[ir] << ' ';
  deduce << endl;
  for(ir=0; ir<nref; ir++) {
    deduce << refpat[ir].fId[0] << ' ' << refpat[ir].fId[1] << ' ' << refpat[ir].fId[2] << ' ' << refpat[ir].fp[0] <<' ' << refpat[ir].fp[1] <<" error " << refpat[ir].fError;
    deduce << ' ' << refpat[ir].fh[0] << ' ' << refpat[ir].fh[1] << ' ' << refpat[ir].fhError << endl;
  }
  for(ir=0; ir<nref; ir++) {
    if(refpat[ir].fError < totalerror*1.e-3) {
      refpat[ir].fp[1] = -1;
      refpat[ir].fp[0] = originalp+1;
    }
  }

  deduce << "originalp " << endl;
  for(ir=0; ir<nref; ir++) {
    deduce << refpat[ir].fId[0] << ' ' << refpat[ir].fId[1] << ' ' << refpat[ir].fId[2] << ' ' << refpat[ir].fp[0] <<' ' << refpat[ir].fp[1] <<" error " << refpat[ir].fError;
    deduce << ' ' << refpat[ir].fh[0] << ' ' << refpat[ir].fh[1] << ' ' << refpat[ir].fhError << endl;
  }
  // for each corner id, identify the edges which are connected to it
  // if any edge suggests an h-refinement, use the h-refinement
  // if no edge suggests an h-refinement return a unique parameter p in porders
  int pref = 1;
  for(ir=0; ir<nref; ir++) if(refpat[ir].fp[1] != -1) pref = 0;
  // if pref is still equal 1, we will use the prefinement
  // determine the maximum p-order suggested and return it.
  if(pref) {
    int maxp = 0;
    for(ir=0; ir<nref; ir++) {
      if(refpat[ir].fp[0] > maxp) maxp = refpat[ir].fp[0];
      if(refpat[ir].fp[1] > maxp) maxp = refpat[ir].fp[1];
    }
    porders.Fill(-1);
    porders[0] = maxp;
    deduce << "prefinement order " << maxp << endl;
    return;
  }
  TPZVec<int> perm(refpat.NElements());
  TPZVec<REAL> error(refpat.NElements());
  totalerror = 0.;
  for(ir=0; ir<nref; ir++) {
    error[ir] = refpat[ir].fhError;
    totalerror += refpat[ir].fhError;
    perm[ir] = ir;
  }
  Sort(error,perm);

  // h-refinement will be used
  // determine the order of interpolation of the sub elements
  int ic;
  deduce << "errors in their order ";
  for(ir=0; ir<nref; ir++) deduce << perm[ir] << ' ' << refpat[perm[ir]].fhError << ' ';
  deduce << endl;
  deduce << "h-refinement order of the elements ";
  for(ic=0; ic<ncorners; ic++) porders[ic]=0;
  for(ir=0; ir<nref; ir++) {
    if(refpat[perm[ir]].fhError < 1.e-3*totalerror) continue;
    for(ic=0; ic<ncorners; ic++) {
      int id = cornerids[ic];
      if(refpat[perm[ir]].fId[0] == id) {
	if(porders[ic] == 0) porders[ic] = refpat[perm[ir]].fh[0];
      }
      if(refpat[perm[ir]].fId[2] == id) {
	if(porders[ic] == 0) porders[ic] = refpat[perm[ir]].fh[1];
      }
    }
  }
  for(ic=0; ic<ncorners; ic++) {
    if(porders[ic] == 0) porders[ic] = originalp/2+1;
    deduce << porders[ic] << ' ';
  }
  deduce << endl;
}

TPZCompMesh *TPZMGAnalysis::CreateCompMesh(TPZCompMesh *mesh, TPZVec<TPZGeoEl *> &gelstack, TPZVec<int> &porders) {

  TPZGeoMesh *gmesh = mesh->Reference();
  if(!gmesh) {
    cout << "TPZMGAnalysis::UniformlyRefineMesh encountered no geometric mesh\n";
    return 0;
  }
  gmesh->ResetReference();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  int nmat = mesh->MaterialVec().NElements();
  int m;
  for(m=0; m<nmat; m++) {
    TPZMaterial *mat = mesh->MaterialVec()[m];
    if(!mat) continue;
    mat->Clone(cmesh->MaterialVec());
  }
  int el,nelem = gelstack.NElements();
  for(el=0; el<nelem; el++) {

    TPZGeoEl *gel = gelstack[el];
    if(!gel) {
      cout << "TPZMGAnalysis::CreateCompMesh encountered an null element\n";
      continue;
    }
    int celindex;
    TPZInterpolatedElement *csint;
    csint = dynamic_cast<TPZInterpolatedElement *> (gel->CreateCompEl(*cmesh,celindex));
    if(!csint) continue;
    csint->PRefine(porders[el]);
  }
  cmesh->AdjustBoundaryElements();
  return cmesh;

}


#include "pzvec.h"

template class TPZVec<TPZMGAnalysis::TPZRefPattern>;
