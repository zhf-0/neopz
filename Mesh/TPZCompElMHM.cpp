/*
 *  TPZCompElMHM.cpp
 *  PZ
 *
 *  Created by Frederico on 30/01/13.
 *  Copyright 2011 LNCC. All rights reserved.
 *
 */

#include <ostream>

#include "TPZCompElMHM.h"

using namespace pzshape;
using namespace pzgeom;

template class TPZGeoElRefPattern<TPZGeoTriangle>;
template class TPZGeoElRefPattern<TPZGeoPoint>;
template class TPZGeoElRefPattern<TPZGeoLinear>;
template class TPZGeoElRefPattern<TPZGeoQuad>;
template class TPZGeoElRefPattern<TPZGeoTetrahedra>;
template class TPZGeoElRefPattern<TPZGeoPrism>;
template class TPZGeoElRefPattern<TPZGeoPyramid>;
template class TPZGeoElRefPattern<TPZGeoCube>;

int block_desired_size = 0;

template <class TSHAPE, class TGEO>
TPZCompElMHM<TSHAPE, TGEO>::TPZCompElMHM() :
TPZIntelGen<TSHAPE>(),
LocalBasisComputed(false),
LocalMaterial(51),
LocalInterpolationOrder(2),
FaceInterpolationOrder(1),
LocalRefinement(1), // must be 1 at least!
NumberOfLocalProblems(0),
compElFaces(0),
geoElFaces(0)
{
    // Nothing to do!
}

template <class TSHAPE, class TGEO>
TPZCompElMHM<TSHAPE,TGEO>::TPZCompElMHM(TPZCompMesh &cmesh, TPZGeoEl *ref, int &index) :
TPZIntelGen<TSHAPE>(cmesh, ref, index),
LocalBasisComputed(false),
LocalMaterial(51),
LocalInterpolationOrder(2),
FaceInterpolationOrder(1),
LocalRefinement(1), // must be 1 at least!
NumberOfLocalProblems(0),
compElFaces(0),
geoElFaces(0)
{
    for(int i = 0; i < TSHAPE::NSides; i++)
    {
        if (TSHAPE::SideDimension(i) != cmesh.Dimension() )
            continue;

        cmesh.ConnectVec()[ this->ConnectIndex(i) ].SetNState(0);
    }
}

template <class TSHAPE, class TGEO>
TPZCompElMHM<TSHAPE,TGEO>::~TPZCompElMHM()
{
    // Nothing to do!
}

template <class TSHAPE, class TGEO>
void TPZCompElMHM<TSHAPE,TGEO>::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix &jacobian, TPZFMatrix &axes,
                  REAL &detjac, TPZFMatrix &jacinv, TPZFMatrix &phi, TPZFMatrix &dphix)
{
    // for all computational elements
    TPZInterpolationSpace::ComputeShape(intpoint, X, jacobian, axes, detjac, jacinv, phi, dphix);

    if ( TSHAPE::Dimension != this->Mesh()->Dimension() )
        return;

    // for computational elements which evaluate local basis.
    if( LocalBasisComputed == false )
    {
        LocalBasisComputed = true;

        TPZCompMesh* cmesh = this->Mesh();

        // Current geometric element.
        TPZGeoEl* ref = this->Reference();

        TPZCreateApproximationSpace fCreate;
        int meshdim = cmesh->Dimension();

        // ====================
        // Local geometric mesh
        // ====================

        int id = 0;
        int index, nsides;
        LocalGeoMesh = new TPZGeoMesh();
        TPZVec<int> element( ref->NNodes(), 0 );

        for(int i = 0; i < ref->NNodes(); i++)
        {
            TPZVec<REAL> coord(3, 0.0);
            ref->NodePtr(i)->GetCoordinates(coord);
            int nodeId = LocalGeoMesh->NodeVec().AllocateNewElement();
            LocalGeoMesh->NodeVec()[nodeId].Initialize(coord, *LocalGeoMesh);
            element[i] = nodeId;
        }

        new TPZGeoElRefPattern<TGEO>(id, element, LocalMaterial, *LocalGeoMesh);

        nsides = LocalGeoMesh->ElementVec()[0]->NSides();
        for(int i = 0; i < nsides; i++)
        {
            int sidedim = LocalGeoMesh->ElementVec()[0]->SideDimension(i);

            if (sidedim != meshdim-1)
                continue;

            LocalGeoMesh->ElementVec()[0]->CreateBCGeoEl(i, -LocalMaterial);
        }

        LocalGeoMesh->BuildConnectivity();

        for( int nref = 0; nref < LocalRefinement; nref++ )
        {
            TPZVec<TPZGeoEl *> filhos;
            int h = LocalGeoMesh->NElements();
            for ( int i = 0; i < h; i++ )
            {
                TPZGeoEl * gel = LocalGeoMesh->ElementVec()[i];
                gel->Divide ( filhos );
            }
        }

        std::ostringstream name_in;
        name_in << this->Index();
        std::string name = "localGeoMesh" + name_in.str() + ".txt";
        //LocalGeoMesh->SetName( name.c_str() );
        //std::ofstream gm( name.c_str() );
        //LocalGeoMesh->Print(gm);

        // ==============================
        // Computational element on faces
        // ==============================

        int is;
        nsides = ref->NSides();
        int cmeshOrder = cmesh->GetDefaultOrder();
        cmesh->SetDefaultOrder(FaceInterpolationOrder);
        for (is = 0; is < nsides; ++is)
        {
            // =======================
            // Geometric element faces: over the dimension-1 sides.
            // =======================

            if (ref->SideDimension(is) != meshdim-1)
                continue;

            TPZGeoEl* geo_temp = ref->CreateBCGeoEl(is, LocalMaterial);
            if (!geo_temp)
                continue;

            geo_temp->SetMaterialId(LocalMaterial);
            geoElFaces.push_back( geo_temp );

            // ==============================
            // Computational element on faces
            // ==============================

            //TPZCompEl::SetgOrder(FaceInterpolationOrder);
            TPZCompEl* comp_temp = fCreate.CreateCompEl( geo_temp, *cmesh, index);
            compElFaces.push_back( comp_temp );
        }
        cmesh->SetDefaultOrder(cmeshOrder);

        // ========================
        // Number of local problems
        // ========================

        int dofPerFace = 0;
        int NumberOfFaces = compElFaces.size();
        for(int j = 0; j < NumberOfFaces; j++)
        {
            dofPerFace = 0;
            int nconnects = compElFaces[j]->NConnects();
            for(int i = 0; i < nconnects; i++)
            {
                dofPerFace += compElFaces[j]->Connect(i).NDof();
                NumberOfLocalProblems += compElFaces[j]->Connect(i).NDof();
            }
        }
        ++NumberOfLocalProblems;

		std::cout << "Number of local problems: " << NumberOfLocalProblems << std::endl;
        LocalCompMesh.Resize(NumberOfLocalProblems);
        
        // ======================
        // 1st computational mesh (Continuous elements)
        // ======================

        for(int loc = 0; loc < NumberOfLocalProblems; loc++)
        {
            LocalGeoMesh->ResetReference();
            LocalCompMesh[loc] = new TPZCompMesh(LocalGeoMesh);
            LocalCompMesh[loc]->SetDefaultOrder(LocalInterpolationOrder);
            LocalCompMesh[loc]->SetDimModel(meshdim);
            LocalCompMesh[loc]->SetAllCreateFunctionsContinuous();
            
            LocalMHM* mat_cmesh1 = new LocalMHM(LocalMaterial, meshdim);    // fake...
            TPZAutoPointer<TPZMaterial> auto_mat_cmesh1(mat_cmesh1);
            LocalCompMesh[loc]->InsertMaterialObject(mat_cmesh1);
            
            int nstates = mat_cmesh1->NStateVariables();
            TPZFMatrix val1(nstates, nstates, 0.), val2(nstates, 1, 0.);
            TPZMaterial * BCondN1 = mat_cmesh1->CreateBC(auto_mat_cmesh1, -LocalMaterial, 1 /*Neumann*/, val1, val2);
            LocalCompMesh[loc]->InsertMaterialObject(BCondN1);
            
            name = "cmesh1_" + name_in.str() + ".txt";
            
            LocalCompMesh[loc]->SetName(name);
            LocalCompMesh[loc]->AutoBuild();
            LocalCompMesh[loc]->AdjustBoundaryElements();
            
            //std::ofstream cm1( name.c_str() );
            //LocalCompMesh[loc]->Print( cm1 );
        }

        // ======================
        // 2nd computational mesh (Lagrange multiplier, once we are dealing with pure Neumann boundary conditions)
        // ======================

        LocalGeoMesh->ResetReference();
        TPZCompEl::SetgOrder(0);
        TPZCompMesh * cmesh2 = new TPZCompMesh(LocalGeoMesh);
        cmesh2->SetDimModel( meshdim );
        cmesh2->SetAllCreateFunctionsDiscontinuous();
        LocalMHM* mat_cmesh2 = new LocalMHM(LocalMaterial, meshdim);  // fake...
        cmesh2->InsertMaterialObject( mat_cmesh2 );
        int nstates = mat_cmesh2->NStateVariables();
        TPZFMatrix val12(nstates, nstates, 0.), val22(nstates, 1, 0.);
        TPZAutoPointer<TPZMaterial> auto_mat_cmesh2 (mat_cmesh2);
        TPZBndCond *fake_bnd2 = mat_cmesh2->CreateBC( auto_mat_cmesh2, -LocalMaterial, 1 /*Neumann*/, val12, val22 );
        cmesh2->InsertMaterialObject(fake_bnd2);

        name = "cmesh2_" + name_in.str() + ".txt";

        cmesh2->SetName(name);
        cmesh2->AutoBuild();
        RemoveInterfaces( cmesh2 );
        cmesh2->AdjustBoundaryElements();
        
        //std::ofstream cm2( name.c_str() );
        //cmesh2->Print( cm2 );

        // ============================
        // Definition of Local Problems
        // ============================

        LocalMHM* localProblem = new LocalMHM(LocalMaterial, meshdim);
        localProblem->SetFaces(compElFaces, geoElFaces);
        nstates = localProblem->NStateVariables();
        TPZAutoPointer<TPZMaterial> autoLocalProblem (localProblem);
        TPZFMatrix val13(nstates, nstates, 0.), val23(nstates, 1, 0.);
        TPZBndCond *bnd3 = localProblem->CreateBC( autoLocalProblem, -LocalMaterial, 1 /*Neumann*/, val13, val23 );

        // ==========================
        // Solution of local problems (p^lambda)
        // ==========================

        for(int loc = 0; loc < NumberOfLocalProblems; loc++)
        {
            // =================
            // Multiphysics mesh (it is required to deal with Lagrange multipliers)
            // =================

            TPZVec<TPZCompMesh*> meshvec(2);
            meshvec[0] = LocalCompMesh[loc];
            meshvec[1] = cmesh2;

            LocalGeoMesh->ResetReference();
            TPZCompMesh* mphysics = new TPZCompMesh(LocalGeoMesh);
            mphysics->SetDimModel(meshdim);
            mphysics->SetDefaultOrder(LocalInterpolationOrder);
            mphysics->SetAllCreateFunctionsMultiphysicElem();

            mphysics->InsertMaterialObject( localProblem );
            mphysics->InsertMaterialObject(bnd3);

            mphysics->AutoBuild();
            mphysics->AdjustBoundaryElements();

            ConfigureLangrangeCompMesh(mphysics, meshvec);

            //name = "localCompMesh" + name_in.str() + ".txt";
            //std::ofstream cm( name.c_str() );
            //mphysics->Print( cm );

            if(loc == (NumberOfLocalProblems-1))
                localProblem->SetLambdaProblem(false);
            else
                localProblem->SetLambdaProblem(true);

            int curDof = loc % dofPerFace;
            int curEdge = (loc - curDof) / dofPerFace;

            TPZAnalysis localAn( mphysics );
            TPZFStructMatrix localMatrix( mphysics );
            TPZStepSolver localSolver;
            localSolver.SetDirect( ELU );
            localAn.SetStructuralMatrix( localMatrix );
            localAn.SetSolver( localSolver );

            localProblem->SetProblemParams(curEdge, curDof);

            localAn.Assemble();

            ResizeAndSolveSystem( localAn );
            mphysics->LoadSolution( localAn.Solution() );

            TPZVec<std::string> scalnames(2), vecnames(0);
            scalnames[0] = "Solution1";
            scalnames[1] = "Solution2";

            std::string solname;
            std::ostringstream sol_in;
            sol_in << curDof;
            std::ostringstream faceOrder;
            faceOrder << curEdge;
            
            if(loc == (NumberOfLocalProblems-1))
                solname = "Solution" + name_in.str() + "_F.vtk";
            else
                solname = "Solution" + name_in.str() + "_Lambda_Face" + faceOrder.str() + "_Dof" + sol_in.str() + ".vtk";

            std::cout << solname << std::endl;

            // Printing informations
            //std::string printname;
            //printname = "Rigidez" + name_in.str() + "_" + sol_in.str() + ".nb";
            //std::ofstream file_stiff2( printname.c_str() );
            //localAn.Solver().Matrix()->Print("Rigidez = ", file_stiff2, EMathematicaInput);

            //printname = "Rhs" + name_in.str() + "_" + sol_in.str() + ".nb";
            //std::ofstream file_rhs2( printname.c_str() );
            //localAn.Rhs().Print("Rhs = ", file_rhs2, EMathematicaInput);

            //printname = "SolutionVec" + name_in.str() + "_" + sol_in.str() + ".nb";
            //std::ofstream file_sol( printname.c_str() );
            //localAn.Solution().Print("Solution = ", file_sol, EMathematicaInput);

            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);

            // Post-processing
            localAn.DefineGraphMesh(meshdim, scalnames, vecnames, solname);
            localAn.PostProcess(0, meshdim);

			std::cout << "=====================================================================" << std::endl;

            // armazenar solução de cada problema local...
        }
    }

	std::cout << "\n -------------------------------------------------------------------- \n" << std::endl;

    TPZSolVec sol;
    TPZGradSolVec dsol;

    this->ComputeSolution(intpoint, sol, dsol, axes);

    // sol e dsol devem ser inseridos em phi e dphi.
    phi.Resize(NumberOfLocalProblems, 1);
    dphix.Resize(this->Mesh()->Dimension(), NumberOfLocalProblems);

    for(int loc = 0; loc < NumberOfLocalProblems; loc++)
    {
        phi(loc,0) = sol[loc][0];
        for(int d = 0; d < this->Mesh()->Dimension(); d++)
            dphix(d, loc) = dsol[loc](d, 0);
    }
}

template <class TSHAPE, class TGEO>
void TPZCompElMHM<TSHAPE,TGEO>::ComputeSolution(TPZVec<REAL> &intpoint, TPZSolVec &sol, TPZGradSolVec &dsol, TPZFMatrix& axes)
{
    sol.Resize(NumberOfLocalProblems);
    dsol.Resize(NumberOfLocalProblems);

    TPZVec<REAL> real(3, 0.0);
    TPZGeoEl* ref = this->Reference();
    ref->X(intpoint, real);

    // std::cout << "X: " << intpoint[0] << ", " << intpoint[1] << std::endl;
    // std::cout << "real: " << real[0] << ", " << real[1] << std::endl;

    for(int loc = 0; loc < NumberOfLocalProblems; loc++)
    {
        sol[loc].Resize(1);
        sol[loc].Fill(0.0);
        dsol[loc].Resize(LocalCompMesh[loc]->Dimension(), 1);
        dsol[loc].Zero();

        LocalCompMesh[loc]->LoadReferences();
        int nel = LocalCompMesh[loc]->NElements();

        TPZSolVec temp_sol;
        TPZGradSolVec temp_dsol;
        temp_sol.Resize(1);
        temp_sol[0].Resize(1, 0.0);
        temp_dsol.Resize(1);
        temp_dsol[0].Resize(LocalCompMesh[loc]->Dimension(), 1);

        for(int i = 0; i < nel; i++)
        {
            TPZCompEl* loc_ref = LocalCompMesh[loc]->ElementVec()[i];

            // this part is tailored for integration on faces.
            if( loc_ref->Dimension() == LocalCompMesh[loc]->Dimension() )
                continue;

            if ( IsXOnTheFace(loc_ref->Reference(), real) )
            {
                loc_ref->ComputeSolution(intpoint, temp_sol, temp_dsol, axes);
            }
        }

        sol[loc][0] = temp_sol[0][0];
        for(int d = 0; d < LocalCompMesh[loc]->Dimension(); d++)
            dsol[loc](d,0) = temp_dsol[0](d,0);
    }
}

template <class TSHAPE, class TGEO>
void TPZCompElMHM<TSHAPE,TGEO>::ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi)
{
    // Under construction...
}

template <class TSHAPE, class TGEO>
void TPZCompElMHM<TSHAPE,TGEO>::RemoveInterfaces(TPZCompMesh* cmesh)
{
    int n = cmesh->ElementVec().NElements();
    for(int i = 0; i < n; i++)
    {
        TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc*>( cmesh->ElementVec()[i] );
        if(!disc)
            continue;
        disc->RemoveInterfaces();
    }
}

template <class TSHAPE, class TGEO>
void TPZCompElMHM<TSHAPE,TGEO>::ConfigureLangrangeCompMesh(TPZCompMesh* MFMesh, TPZVec<TPZCompMesh*> cmeshVec)
{
    // ============
    // Add Elements
    // ============
    
	TPZGeoMesh *gmesh = MFMesh->Reference();
	gmesh->ResetReference();
	int nMFEl = MFMesh->NElements();
	int nmesh = cmeshVec.size();
	int imesh;
	for(imesh = 0; imesh < nmesh; imesh++)
	{
		cmeshVec[imesh]->LoadReferences();
		int iel;
		for(iel = 0; iel < nMFEl; iel++)
		{
            TPZCompEl *cel = MFMesh->ElementVec()[iel];
			TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *> (cel);
			if(mfcel)
			{
                int found = 0;
                TPZGeoEl *gel = mfcel->Reference();
                //std::cout << "imesh: " << imesh << ", iel: " << iel << ", sides: " << gel->NSides() << ", matId: " << gel->MaterialId() << std::endl;
                TPZStack<TPZCompElSide> celstack;
                TPZGeoElSide gelside(gel, gel->NSides()-1);
                // if the geometric element has a reference, it is an obvious candidate
                if (gel->Reference()) 
                    celstack.Push(gelside.Reference());
                
                // put all large and small elements on the stack
                gelside.ConnectedCompElementList(celstack, 0, 1);
                while (celstack.size())
                {
                    int ncel = celstack.size();
                    // the last element on the stack
                    TPZGeoElSide gelside = celstack[ncel-1].Reference();
                    TPZStack<TPZCompElSide> celstack2;
                    // put te last element on the new stack
                    celstack2.Push(celstack[ncel-1]);
                    celstack.Pop();
                    // put all equal level elements on the new stack
                    gelside.EqualLevelCompElementList(celstack2, 0, 1);
                    
                    while (celstack2.size()) 
                    {
                        int ncel2 = celstack2.size();
                        TPZGeoElSide gelside2 = celstack2[ncel2-1].Reference();
                        
                        // put all elements in the stack - if there is one element, stop the search
                        // if(gelside2.Element()->Dimension()==gel->Dimension())
                        if( gelside2.Element()->Reference() )
                        {
                            //std::cout << "iel: " << iel << std::endl;
                            mfcel->AddElement(gelside2.Element()->Reference(), imesh);
                            found = 1;
                            celstack2.Resize(0);
                            celstack.Resize(0);
                            break;
                        }
                        
                        celstack2.Pop();
                    }
                }

                if (!found)
                {
                    std::cout << "NOT_FOUND - imesh: " << imesh << ", iel: " << iel << ", sides: " << gel->NSides() << ", matId: " << gel->MaterialId() << std::endl;
                    mfcel->AddElement(0, imesh);
                }
            }
		}
		gmesh->ResetReference();
	}

    // ============
    // Add Connects
    // ============

	int nmeshes = cmeshVec.size();
	TPZVec<int> FirstConnect(nmeshes,0);
	int nconnects = 0;
	for (imesh = 0; imesh < nmeshes; imesh++)
	{
		FirstConnect[imesh] = nconnects;
        nconnects += cmeshVec[imesh]->ConnectVec().NElements();
	}
	MFMesh->ConnectVec().Resize(nconnects);
	MFMesh->Block().SetNBlocks(nconnects);
	int counter = 0;
    block_desired_size = 0;
    int seqnum2 = 0;
	for (imesh = 0; imesh < nmeshes; imesh++)
	{
		int ic;
		int nc = cmeshVec[imesh]->ConnectVec().NElements();

		for (ic = 0; ic < nc; ic++)
		{
			TPZConnect &refcon =  cmeshVec[imesh]->ConnectVec()[ic];
			MFMesh->ConnectVec()[counter] = refcon;

			if (refcon.SequenceNumber() >= 0)
            {
				MFMesh->ConnectVec()[counter].SetSequenceNumber(block_desired_size);
				MFMesh->ConnectVec()[counter].SetNState(refcon.NState());
				MFMesh->ConnectVec()[counter].SetNShape(refcon.NShape());
				int ndof = refcon.NDof(*cmeshVec[imesh]);
				MFMesh->Block().Set(seqnum2,ndof);
                if( imesh == 0)
                    block_desired_size++;
                seqnum2++;
			}
			counter++;
		}

		///ajustar as dependencias
		for (ic = 0; ic < nc; ic++)
		{
			TPZConnect &cn = MFMesh->ConnectVec()[FirstConnect[imesh]+ic];
			if (cn.HasDependency()) 
			{
				TPZConnect::TPZDepend *dep = cn.FirstDepend();
				while (dep)
                {
					dep->fDepConnectIndex = dep->fDepConnectIndex+FirstConnect[imesh];
					dep = dep->fNext;
				}
			}
		}
	}

	MFMesh->Block().SetNBlocks(seqnum2);    
	MFMesh->ExpandSolution();
	int iel;
	int nelem = MFMesh->NElements();
	for (iel = 0; iel < nelem; iel++) 
	{
		TPZMultiphysicsElement *cel = dynamic_cast<TPZMultiphysicsElement *> (MFMesh->ElementVec()[iel]);
		if (!cel)
            continue;

		TPZStack<int> connectindexes;
		int imesh;
		for (imesh=0; imesh < nmeshes; imesh++)
        {
			TPZCompEl *celref = cel->ReferredElement(imesh);
            if (!celref) 
                continue;
            
			int ncon = celref->NConnects();
			int ic;
			for (ic=0; ic<ncon; ic++) {
				connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[imesh]);
			}
		}
		cel->SetConnectIndexes(connectindexes);
	}

    // ==================
    // TransferFromMeshes
    // ==================

    nmeshes = cmeshVec.size();
    TPZManVector<int> FirstConnectIndex(nmeshes+1, 0);
    for (imesh = 0; imesh < nmeshes; imesh++)
    {
        FirstConnectIndex[imesh+1] = FirstConnectIndex[imesh] + cmeshVec[imesh]->NConnects();
    }

    TPZBlock &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++)
    {
        int ncon = cmeshVec[imesh]->NConnects();
        
		TPZBlock &block = cmeshVec[imesh]->Block();
		int ic;
		for (ic = 0; ic < ncon; ic++)
        {
			TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
			int seqnum = con.SequenceNumber();

			if(seqnum < 0)  // Whether connect was deleted by previous refined process
                continue;

			int blsize = block.Size(seqnum);
			TPZConnect &conMF = MFMesh->ConnectVec()[FirstConnectIndex[imesh]+ic];
			int seqnumMF = conMF.SequenceNumber();
            // std::cout << "imesh: " << imesh << ", ic: " << ic << ", seqnum: " << seqnum << ", blsize: " << blsize << ", seqnumMF: " << seqnumMF << "\n";
			int idf;
			for (idf = 0; idf < blsize; idf++)
            {
				blockMF.Put(seqnumMF, idf, 0, block.Get(seqnum, idf, 0) );
			}
		}
	}
}

template <class TSHAPE, class TGEO>
void TPZCompElMHM<TSHAPE,TGEO>::ResizeAndSolveSystem(TPZAnalysis& an)
{
    // ======
    // Resize
    // ======

    int nelem = an.Mesh()->NElements();
    int nelem_dim = 0;
    int elem_dim = 0;
    int mesh_dim = an.Mesh()->Dimension();
    for(int i = 0; i < nelem; i++)
    {
        elem_dim = an.Mesh()->ElementVec()[i]->Dimension();
        
        if( elem_dim == mesh_dim )
            nelem_dim++;
    }

	int old_size = an.Mesh()->NEquations();
    int new_size = old_size - nelem_dim + 1;

    an.Solver().Matrix()->Resize(new_size, new_size);
    an.Solution().Redim(new_size, 1);
    an.Rhs().Resize(new_size, 1);

    // =====
    // Solve
    // =====

    TPZFMatrix residual( an.Rhs() );
    TPZFMatrix delu( an.Solution() );

    an.Solver().Solve(residual, delu);

    an.LoadSolution(delu);
}

template <class TSHAPE, class TGEO>
bool TPZCompElMHM<TSHAPE,TGEO>::IsXOnTheFace( TPZGeoEl* geo, TPZVec<REAL>& X )
{
    REAL tol = 10E-6;
    TPZGeoNode *np;
    int nnodes = geo->NNodes();
    TPZFMatrix coord(3, nnodes, 0.0);
    
    for(int i = 0; i < nnodes; i++)
    {
        np = geo->NodePtr(i);
        for(int j = 0; j < 3; j++)
            coord(j,i) = np->Coord(j);
    }
    
    int dim = geo->Dimension();
    REAL line = sqrt( pow(coord(0,1) - coord(0,0),2) + pow(coord(1,1) - coord(1,0),2) );
    REAL dis = 0;
    switch (dim)
    {
        case (1) :
            // distance between point and line
            dis = fabs( (coord(0,1) - coord(0,0))*(coord(1,0) - X[1]) - (coord(0,0) - X[0])*(coord(1,1) - coord(1,0)) ) / line;
            break;
        case (2) :
            // distance between point and plane.
            dis = 1.0;
            break;
        default:
            DebugStop();
            break;
    }
    
    if( fabs(dis) < tol )
        return true;
    else 
        return false;
}

// ============================================================================

template class TPZCompElMHM<TPZShapeTriang, TPZGeoTriangle>;
template class TPZCompElMHM<TPZShapePoint, TPZGeoPoint>;
template class TPZCompElMHM<TPZShapeLinear, TPZGeoLinear>;
template class TPZCompElMHM<TPZShapeQuad, TPZGeoQuad>;
template class TPZCompElMHM<TPZShapeTetra, TPZGeoTetrahedra>;
template class TPZCompElMHM<TPZShapePrism, TPZGeoPrism>;
template class TPZCompElMHM<TPZShapePiram, TPZGeoPyramid>;
template class TPZCompElMHM<TPZShapeCube, TPZGeoCube>;

TPZCompEl * CreateMHMPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZCompElMHM<TPZShapePoint, TPZGeoPoint>(mesh, gel, index);
}

TPZCompEl * CreateMHMLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZCompElMHM<TPZShapeLinear, TPZGeoLinear>(mesh, gel, index);
}

TPZCompEl * CreateMHMTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZCompElMHM<TPZShapeTriang, TPZGeoTriangle>(mesh, gel, index);
}

TPZCompEl * CreateMHMQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZCompElMHM<TPZShapeQuad, TPZGeoQuad>(mesh, gel, index);
}

TPZCompEl * CreateMHMCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZCompElMHM<TPZShapeCube, TPZGeoCube>(mesh, gel, index);
}

TPZCompEl * CreateMHMPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZCompElMHM<TPZShapePrism, TPZGeoPrism>(mesh, gel, index);
}

TPZCompEl * CreateMHMTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZCompElMHM<TPZShapeTetra, TPZGeoTetrahedra>(mesh, gel, index);
}

TPZCompEl * CreateMHMPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZCompElMHM<TPZShapePiram, TPZGeoPyramid>(mesh, gel, index);
}
