#include "pzlog.h"
#include "TPZFracAnalysis.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

TPZFracAnalysis::TPZFracAnalysis(TPZAutoPointer<TPZFracData> Data)
{
  fData = Data;
  fmeshvec.Resize(2);
  fgmesh = NULL;
  fcmeshMixed = NULL;
  for (int i = 0; i < 2; i++) {
    fmeshvec[i] = NULL;
  }
  fLastStepRhs.Redim(0, 0);
  fMatFrac = NULL;
  fmustStop = false;
}

TPZFracAnalysis::~TPZFracAnalysis()
{
  for (int i = 0; i < 2; i++) {
    delete fmeshvec[i];
  }
  delete fcmeshMixed;
  delete fgmesh;
  /*
  if (fMatFrac) {
    delete fMatFrac;
  }
  */
}

void TPZFracAnalysis::Run()
{
  
  // Creating empty gmesh
  fgmesh = new TPZGeoMesh;
  
  REAL vl = this->RunUntilOpen();

  const REAL lFrac = fData->ElSize();
  fData->SetLfrac(lFrac);
  
  // Criando primeira malha geometrica com um elemento
  this->CreateFirstGeoElWithBC();
  
  // Malhas computacionais - FluxoH1, PressaoL2, Multifisica para sistema misto quente
  fmeshvec[0] = CreateCMeshFluxH1();
  fmeshvec[1] = CreateCMeshPressureL2();
  TPZFMatrix<> vlMatrix(1,1,vl);
  fcmeshMixed = CreateCMeshMixed(vlMatrix);
  
  // Analysis
  bool mustOptimizeBandwidth = false;
  TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
  TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  an->SetSolver(step);
  an->SetStructuralMatrix(skyl);
  
  // Plot of first solution
  this->PostProcessVTK(an);
  
  // Solving transiente system
  while (fmustStop == false) {
    bool propagate = SolveSistTransient(an);
    if (propagate) {
      
      // Novo comprimento de fratura
      REAL newLfrac = fData->Lfrac() + fData->ElSize();
      fData->SetLfrac(newLfrac);
      std::cout << "Lfrac = " << newLfrac << std::endl;
      
      // Apagando elemento de contorno que impoe pressao
      TPZGeoEl *gel = this->FindPressureBCElement();
      fgmesh->ResetReference();
      fmeshvec[0]->LoadReferences();
      delete gel->Reference();
      fgmesh->ResetReference();
      fmeshvec[1]->LoadReferences();
      delete gel->Reference();
      fgmesh->ResetReference();
      fcmeshMixed->LoadReferences();
      gel->Reference()->SetFreeIntPtIndices();
      delete gel->Reference();
      gel->RemoveConnectivities();
      delete gel;
      
      fmeshvec[0]->CleanUpUnconnectedNodes();
      fmeshvec[1]->CleanUpUnconnectedNodes();
      fcmeshMixed->CleanUpUnconnectedNodes();

      // Criando novo noh
      TPZAdmChunkVector<TPZGeoNode> &nodevec = fgmesh->NodeVec();
      const int nnodes = nodevec.NElements() + 1;
      nodevec.Resize(nnodes);
      const int lastnodeid = nnodes - 1;
      nodevec[lastnodeid].SetNodeId(lastnodeid);
      TPZVec<REAL> x(3,0.);
      x[0] = newLfrac;
      nodevec[lastnodeid].SetCoord(x);
      
      // Criando novo elemento geometrico
      TPZVec<long> TopolLine(2);
      TopolLine[0] = nnodes - 2;
      TopolLine[1] = nnodes - 1;
      long index;
      const int matid = 1;
      TPZGeoEl *newGel = fgmesh->CreateGeoElement(EOned, TopolLine, matid, index);
      fgmesh->BuildConnectivity();

      // Colocano a cond de contorno de pressao na malha geometrica
      const int bcpressureid = -2;
      TPZGeoEl *gelBCPress = newGel->CreateBCGeoEl(1, bcpressureid);
      fgmesh->BuildConnectivity();
      
      // Criando novo elemento computacional de fratura
      fgmesh->ResetReference();
      fmeshvec[0]->LoadReferences();
      fmeshvec[0]->CreateCompEl(newGel, index);
      fgmesh->ResetReference();
      fmeshvec[1]->LoadReferences();
      TPZCompEl *celPress = fmeshvec[1]->CreateCompEl(newGel, index);
      fgmesh->ResetReference();
      fcmeshMixed->LoadReferences();
      TPZCompEl *celMixed = fcmeshMixed->CreateCompEl(newGel, index);
      
      // Criando novo elemento computacional de BC de pressao
      fgmesh->ResetReference();
      fmeshvec[0]->LoadReferences();
      fmeshvec[0]->CreateCompEl(gelBCPress, index);
      fgmesh->ResetReference();
      fmeshvec[1]->LoadReferences();
      fmeshvec[1]->CreateCompEl(gelBCPress, index);
      fgmesh->ResetReference();
      fcmeshMixed->LoadReferences();
      TPZCompEl *celBCMixed = fcmeshMixed->CreateCompEl(gelBCPress, index);
      
      // Ajustando a estrutura de dados
      fmeshvec[0]->ComputeNodElCon();
      fmeshvec[1]->ComputeNodElCon();
      fcmeshMixed->ComputeNodElCon();
      fmeshvec[0]->CleanUpUnconnectedNodes();
      fmeshvec[1]->CleanUpUnconnectedNodes();
      fcmeshMixed->CleanUpUnconnectedNodes();
      fmeshvec[0]->ExpandSolution();
      fmeshvec[1]->ExpandSolution();
      fcmeshMixed->ExpandSolution();
      
      // Setando a pressao inicial no elemento novo
      TPZConnect &c = celPress->Connect(0);
      TPZBlock<STATE> &Block = fmeshvec[1]->Block();
      int sq = c.SequenceNumber();
      int sz = Block.Size(sq);
      if (sz > 1) {
        DebugStop(); // Only works for p = 0 in pressure space
      }
      int pos = Block.Position(sq);
      fmeshvec[1]->Solution()(pos,0) = fData->SigmaConf();
      
      // Transferindo para a multifisica
      TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshMixed);
      TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshMixed);
      TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
      
#ifdef DEBUG
      std::ofstream out2("meshes2.txt");
      fgmesh->Print(out2);
      fmeshvec[0]->Print(out2);
      fmeshvec[1]->Print(out2);
      fcmeshMixed->Print(out2);
#endif
    
      // Preparando os index dos pontos de integracao.
      TPZFMatrix<> Vl(1,1,fData->AccumVl());
      fMatFrac->SetDefaultMem(Vl);
      celMixed->PrepareIntPtIndices();
      celBCMixed->PrepareIntPtIndices();
      // Vl eh resetado depois de inicializar o chute inicial de newton
      
#ifdef DEBUG
      std::ofstream outvtk("newfrac.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(fgmesh, outvtk,true);
#endif
    }
    
    TPZEquationFilter newEquationFilter(fcmeshMixed->NEquations());
    an->StructMatrix()->EquationFilter() = newEquationFilter; //AQUINATHAN
    
    std::cout << "\nSolucao apos a propagacao:" << std::endl;
    fcmeshMixed->Solution().Print("sol");
    
    an->LoadSolution(fcmeshMixed->Solution());
    this->PostProcessVTK(an);
  }
  
  //fData->PrintDebugMapForMathematica("DebugMapQl.nb");
  
  delete an;
}

void TPZFracAnalysis::RunTest()
{
  
  DebugStop();
  
  // Malha geometrica
  fgmesh = CreateGMesh();
  
  // Malhas computacionais - FluxoH1, PressaoL2, Multifisica para sistema misto quente
  fmeshvec[0] = CreateCMeshFluxH1();
  fmeshvec[1] = CreateCMeshPressureL2();
  TPZFMatrix<> vlZero(1,1,0.);
  fcmeshMixed = CreateCMeshMixed(vlZero);
  
  // Analysis
  bool mustOptimizeBandwidth = false;
  TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
  TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  an->SetSolver(step);
  an->SetStructuralMatrix(skyl);
  
  // Plot of first solution
  this->PostProcessVTK(an);
  
  // Solving transiente system
  SolveSistTransient(an);
  
  //fData->PrintDebugMapForMathematica("DebugMapQl.nb");
  
  delete an;
}

TPZGeoMesh * TPZFracAnalysis::CreateGMesh()
{
  const int nel = fData->NelFrac();
  const REAL lfrac = fData->Lfrac();
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  
  // Nos
  const int nnodes = nel+1;
  const REAL elsize = lfrac / nel;
  gmesh->NodeVec().Resize(nnodes);
  TPZVec<REAL> coord(3,0.);
  for (int in = 0; in < nnodes; in++) {
    coord[0] = in * elsize;
    gmesh->NodeVec()[in].SetNodeId(in);
    gmesh->NodeVec()[in].SetCoord(coord);
  }
  
  // Elementos
  TPZVec<long> TopolLine(2,0);
  long index = 0;
  for (int iel = 0; iel < nel; iel++) {
    TopolLine[0] = iel;
    TopolLine[1] = iel+1;
    gmesh->CreateGeoElement(EOned, TopolLine, matIdFrac, index);
  }
  
  gmesh->BuildConnectivity();
  
  // Left
  TPZVec<long> TopolPoint(1,0);
  gmesh->CreateGeoElement(EPoint, TopolPoint, bcLeftId, index);
  
  // Right
  TopolPoint[0] = nnodes-1;
  gmesh->CreateGeoElement(EPoint,TopolPoint,bcRightId,index);
  
  gmesh->BuildConnectivity();
  
#ifdef DEBUG
  std::ofstream out("GeoMesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
  
  return gmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshFluxH1()
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  const int fluxorder = fData->PorderFlow();
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  cmesh->InsertMaterialObject(mat);
  
  // Condicao de contorno na esquerda
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  // Condicao de contorno na direita
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando H1
  cmesh->SetDimModel(1);
  cmesh->SetDefaultOrder(fluxorder);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();
  
#ifdef DEBUG
  std::ofstream out("CMeshFluxH1.txt");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshPressureL2()
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  const int pressureorder = fData->PorderPressure();
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  cmesh->InsertMaterialObject(mat);
  
  // Condicao de contorno na esquerda
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  
  // Condicao de contorno na direita
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando L2
  cmesh->SetDimModel(1);
  cmesh->SetDefaultOrder(pressureorder);
  cmesh->SetAllCreateFunctionsDiscontinuous();
  
  cmesh->AutoBuild();
  
//  int ncon = cmesh->NConnects();
//  for(int i=0; i<ncon; i++)
//  {
//    TPZConnect &newnod = cmesh->ConnectVec()[i];
//    newnod.SetLagrangeMultiplier(1);
//  }
  
  for (int i = 0; i < cmesh->Solution().Rows(); i++) {
    cmesh->Solution()(i,0) = fData->SigmaConf();
  }
  
#ifdef DEBUG
  std::ofstream out("CMeshPressureL2.txt");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshMixed(TPZFMatrix<STATE> vlMatrix)
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  fMatFrac = new TPZMatfrac1dhdiv(matIdFrac);
  fMatFrac->SetSimulationData(fData);
  cmesh->InsertMaterialObject(fMatFrac);
  fMatFrac->SetDefaultMem(vlMatrix);
  
  // Condicao de contorno na esquerda
  val2(0,0) = fData->Q();
  TPZBndCond * bcLeft = fMatFrac->CreateBC(fMatFrac, bcLeftId, typeFlux, val1, val2);
  cmesh->InsertMaterialObject(bcLeft);
  
  // Condicao de contorno na direita
  val2(0,0) = fData->SigmaConf();
  TPZBndCond * bcRight = fMatFrac->CreateBC(fMatFrac, bcRightId, typePressure, val1, val2);
  cmesh->InsertMaterialObject(bcRight);
  
  // Setando Multifisico
  cmesh->SetDimModel(1);
  cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
  cmesh->AutoBuild();
  
  // Transferindo para a multifisica
	TPZBuildMultiphysicsMesh::AddElements(fmeshvec, cmesh);
	TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, cmesh);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, cmesh);
  
  // Preparando os index dos pontos de integracao.
	long nel = cmesh->NElements();
	for (long iel = 0; iel < nel; iel++) {
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		cel->PrepareIntPtIndices();
	}
  
#ifdef DEBUG
  std::ofstream out("CMeshMultiPhysics.txt");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

void TPZFracAnalysis::IterativeProcess(TPZAnalysis *an, std::ostream &out)
{
	int iter = 0;
	REAL error = 1.e10, NormResLambdaLast = 1.e10;
  const REAL tol = 1.e-8;
  const int numiter = 50;
  
  fData->SetCurrentState();
	int numeq = an->Mesh()->NEquations();
	
	TPZFMatrix<STATE> prevsol(an->Solution());
  TPZFMatrix<STATE> SoliterK(prevsol);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
  
  an->Assemble();
  an->Rhs() += fLastStepRhs;
  TPZAutoPointer< TPZMatrix<REAL> > matK;
  
	while(error > tol && iter < numiter) {
		
		an->Solve(); // o an->Solution() eh o deltaU aqui
    SoliterK = prevsol - an->Solution();
		REAL normDeltaSol = Norm(an->Solution());
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
      std::stringstream sout;
      matK=an->Solver().Matrix();
      matK->Print("matK = ", sout,EMathematicaInput);
      an->Solution().Print("DeltaX = ", sout,EMathematicaInput);
      SoliterK.Print("Xk = ", sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    an->LoadSolution(SoliterK); // Aqui o an->Solution() eh o U
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
    an->Assemble();
    an->Rhs() += fLastStepRhs;
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
      std::stringstream sout;
      an->Rhs().Print("Res = ", sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    double NormResLambda = Norm(an->Rhs());
		double norm = normDeltaSol;
		out << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << std::endl;
    
    //SoliterK.Print("Sol");
    //an->Rhs().Print("Rhs:");
    
		if(norm < tol) {
			out << "\nTolerancia do DELTAU atingida na iteracao : " << (iter+1) << std::endl;
			out << "\n\nNorma da solucao |Delta(Un)|  : " << norm << std::endl << std::endl;
      
		}
    else if( (NormResLambda - NormResLambdaLast) > 1.e-9 ) {
      out << "\nDivergent Method\n";
    }
    
    NormResLambdaLast = NormResLambda;
		error = norm;
		iter++;
    prevsol = SoliterK;
		out.flush();
	}
  
  if (error > tol) {
    DebugStop(); // Metodo nao convergiu!!
  }
  
}

void TPZFracAnalysis::AssembleLastStep(TPZAnalysis *an)
{
  fData->SetLastState();
  
  an->Assemble();
  fLastStepRhs = an->Rhs();
}

bool TPZFracAnalysis::SolveSistTransient(TPZAnalysis *an)
{
  
  bool propagate = false;
  static int itglob = 0;
  int it = 0;
  while (fmustStop == false && propagate == false) {
    
    AssembleLastStep(an);
    TPZFMatrix<> lastSol = an->Solution();
    
    if (it == 0) {

      const long pSolSize = fmeshvec[1]->Solution().Rows();
      if(itglob == 0){
        fmeshvec[0]->Solution()(0,0) = fData->Q();
        
        const REAL pfrac = 0.;
        const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
        const REAL vlnext = fData->VlFtau(pfrac, tstar+fData->TimeStep());
        const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
        const REAL totalLeakOffPrev = 2. * fData->ElSize() * fData->AccumVl();
        const REAL ql = (totalLeakOff - totalLeakOffPrev)/fData->TimeStep();
        REAL qout = fData->Q() - ql;

        fmeshvec[0]->Solution()(1,0) = qout;
        const REAL dwdp = fData->GetDwDp();
        const REAL pini = fData->SigmaConf() + pow(12. * fData->Viscosity() * qout * fData->ElSize() / (dwdp*dwdp*dwdp),1./4.);
        for (int i = 0; i < pSolSize; i++) {
          fmeshvec[1]->Solution()(i,0) = pini;
        }
        
      }
      else{
        fmeshvec[1]->Solution()(pSolSize-1,0) = fmeshvec[1]->Solution()(pSolSize-2,0);
      }
      
      fData->SetAccumVl(0.); // Zerando accumvl para os proximos
      TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
      an->LoadSolution(fcmeshMixed->Solution());
      an->Solution().Print("anini: ");
    }

    IterativeProcess(an, std::cout);
    
    if (0)
    {
      an->Solver().Matrix()->Print("Jac");
      an->Rhs().Print("Rhs");
    }
    
    const REAL qtip = this->Qtip();
    fcmeshMixed->Solution().Print("meshsol");
    if (qtip < 0.) {
      propagate = false;
    }
    else{
      propagate = this->VerifyIfPropagate(qtip);
    }
    
    if (propagate) {
      an->Solution() = lastSol;
      an->LoadSolution();
      TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
      std::cout << "\n******************************** FRACTURE PROPAGATED ********************************" << std::endl;
    }
    else{
      fData->SetNextTime();
      if (qtip > 0.) {
        AcceptSolution(an); // updates leak off
      }
      this->PostProcessVTK(an);
    }
    
    REAL peteleco = 1.e-8;
    if( fData->Time() > (fData->TotalTime() - peteleco) )
    {
      fmustStop = true;
    }
    it++;
    itglob++;
  }
  return propagate;
}

void TPZFracAnalysis::AcceptSolution(TPZAnalysis *an)
{
  fMatFrac->SetUpdateMem();
  an->AssembleResidual();
  fMatFrac->SetUpdateMem(false);
}

void TPZFracAnalysis::PostProcessVTK(TPZAnalysis *an)
{
  const int dim = 1;
  TPZStack<std::string> scalnames, vecnames;
  scalnames.Push("Pressure");
  scalnames.Push("Flow");
  scalnames.Push("Opening");
  an->DefineGraphMesh(dim, scalnames, vecnames, fData->PostProcessFileName());
  an->PostProcess(0);
}

REAL TPZFracAnalysis::Qtip()
{
  fgmesh->ResetReference();
  fcmeshMixed->LoadReferences();
  const int bcRightId = -2;
  const long nel = fcmeshMixed->NElements();
  TPZCompEl *cel = NULL;
  for (long iel = 0; iel < nel; iel++) {
    cel = fcmeshMixed->Element(iel);
    if (!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if (gel->Dimension() != 1) {
      continue;
    }
    TPZGeoElSide gelside(gel,1);
    TPZGeoElSide neigh = gelside.Neighbour();
    while (neigh != gelside) {
      if (neigh.Element()->MaterialId() == bcRightId) {
        break;
      }
      neigh = neigh.Neighbour();
    }
    if (neigh == gelside) {
      continue;
    }
    break;
  }
  
  TPZVec<REAL> qsi(3,1.), sol(1,0.);
  const int varQ = fMatFrac->VariableIndex("Flow");
  cel->Solution(qsi, varQ, sol);
  const REAL qTip = sol[0];
  std::cout << "\nqtip = " << qTip << std::endl;

  return qTip;
}



TPZGeoEl* TPZFracAnalysis::CreateFirstGeoElWithBC()
{
  const int nnodes = 2;
  const int nel = nnodes - 1;
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  
  // Nos
  const REAL elsize = fData->ElSize();
  fgmesh->NodeVec().Resize(nnodes);
  TPZVec<REAL> coord(3,0.);
  for (int in = 0; in < nnodes; in++) {
    coord[0] = in * elsize;
    fgmesh->NodeVec()[in].SetNodeId(in);
    fgmesh->NodeVec()[in].SetCoord(coord);
  }
  
  // Elementos
  TPZVec<long> TopolLine(2,0);
  long index = 0;
  TPZGeoEl *gel = NULL;
  for (int iel = 0; iel < nel; iel++) {
    TopolLine[0] = iel;
    TopolLine[1] = iel+1;
    gel = fgmesh->CreateGeoElement(EOned, TopolLine, matIdFrac, index);
  }
  
  fgmesh->BuildConnectivity();
  
  // Left
  TPZVec<long> TopolPoint(1,0);
  fgmesh->CreateGeoElement(EPoint, TopolPoint, bcLeftId, index);
  
  // Right
  TopolPoint[0] = nnodes-1;
  fgmesh->CreateGeoElement(EPoint,TopolPoint,bcRightId,index);
  
  fgmesh->BuildConnectivity();
  
#ifdef DEBUG
  std::ofstream out("GeoMesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(fgmesh, out, true);
#endif

  return gel;
}

bool TPZFracAnalysis::VerifyIfPropagate(REAL qtip)
{
  const REAL dt = fData->TimeStep();
  const REAL AccumVolThroughTip = fData->AccumVl() * fData->ElSize() * 2.;
  const REAL volThroughTip = qtip * dt + AccumVolThroughTip;
  const REAL pfrac = 0.;
  const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
  REAL vl = fData->VlFtau(pfrac, tstar+dt);
  const REAL totalLeakOff = 2. * fData->ElSize() * vl;
  const REAL totalLeakOffprev = 2. * fData->ElSize() * fData->AccumVl();
  const REAL ql = (totalLeakOff - totalLeakOffprev)/dt;
  if (qtip > 2.5 * ql) { // AQUINATHAN
    return true;
  }
  else{
    vl = fData->AccumVl() + qtip*dt/fData->ElSize()/2.;
    fData->SetAccumVl(vl);
    return false;
  }

}

REAL TPZFracAnalysis::RunUntilOpen()
{
  const int maxinitialit = 10000;
  const REAL qtip = fData->Q();
  REAL vlForThisTime = 0.;
  int it = 0;
  for (it = 0; it < maxinitialit; it++) {
    const REAL dt = fData->TimeStep();
    const REAL AccumVolThroughTip = fData->AccumVl() * fData->ElSize() * 2.;
    const REAL volThroughTip = qtip * dt + AccumVolThroughTip;
    const REAL pfrac = 0.;
    const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
    REAL vlnext = fData->VlFtau(pfrac, tstar+dt);
    const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
    const REAL totalLeakOffPrev = 2. * fData->ElSize() * fData->AccumVl();
    const REAL ql = (totalLeakOff - totalLeakOffPrev)/dt;
    if (qtip > 2.5 * ql) { //AQUINATHAN
      break;
    }
    
    vlnext = fData->AccumVl() + qtip*dt/fData->ElSize()/2.;
    fData->SetAccumVl(vlnext);
    fData->SetNextTime();
  }
  if (it == maxinitialit) {
    DebugStop();
  }
  
  std::cout << "#################### Opening of the fracture occured at time t = " << fData->Time() << " s ####################" << std::endl;
  std::cout << "Total vol injected = " << qtip *fData->Time() << std::endl;
  std::cout << "\nStarting First Simulation" << std::endl;
  
  return fData->AccumVl();
}

TPZGeoEl * TPZFracAnalysis::FindPressureBCElement()
{
  TPZGeoEl *gel = NULL;
  const int bcpressureid = -2;
  for (long iel = 0; iel < fgmesh->NElements(); iel++) {
    gel = fgmesh->ElementVec()[iel];
    if (gel->MaterialId() == bcpressureid) {
      break;
    }
  }
#ifdef DEBUG
  if (gel == NULL) {
    DebugStop();
  }
#endif
  return gel;
}