/*
 *  TPZMultiphyisicsCompElMHM.cpp
 *  PZ
 *
 *  Created by Frederico on 01/02/2013.
 *  Copyright 2013 LNCC. All rights reserved.
 *
 */

#include "TPZMultiphysicsCompElMHM.h"

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pztrnsform.h"
#include "pzmaterial.h"
#include "tpzautopointer.h"
#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeopyramid.h"
#include "pzmaterial.h"
#include "pzelmat.h"
#include "pzconnect.h"
#include "pzmaterialdata.h"
#include "pzinterpolationspace.h"
#include "pzlog.h"

#include <pzshapelinear.h>
#include <pzshapequad.h>
#include <pzshapetriang.h>
#include <pzshapetetra.h>
#include <pzshapecube.h>
#include <pzshapepiram.h>
#include <pzshapepoint.h>
#include <pzshapeprism.h>

#include "TPZCompElMHM.h"

#include <set>

using namespace pzgeom;

template <class TSHAPE, class TGeometry>
TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::TPZMultiphysicsCompElMHM() : TPZMultiphysicsCompEl<TGeometry>()
{

}

template <class TSHAPE, class TGeometry>
TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::TPZMultiphysicsCompElMHM(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) : TPZMultiphysicsCompEl<TGeometry>(mesh, ref, index)
{

}

template <class TSHAPE, class TGeometry>
TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::~TPZMultiphysicsCompElMHM()
{

}

template <class TSHAPE, class TGeometry>
void TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::Solution(TPZVec<REAL> &qsi, int var,TPZVec<REAL> &sol)
{
	if(var >= 100) {
		TPZCompEl::Solution(qsi,var,sol);
		return;
	}
	int dim = this->Dimension();
	
	TPZAutoPointer<TPZMaterial> material = this->Material();
	if(!material){
		sol.Resize(0);
		return;
	}
	
	TPZManVector<TPZTransform> trvec;
	this->AffineTransform(trvec);
	
	//int neqsi= qsi.size();
	TPZVec<REAL> myqsi;
	myqsi.resize(qsi.size());
	
	const int numdof = material->NStateVariables();
	int nref = this->fElementVec.size();
	TPZVec<TPZMaterialData> datavec;
	datavec.resize(nref);
	
	for (int iref = 0; iref<nref; iref++)
	{		
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->fElementVec[iref]);
		datavec[iref].p = msp->MaxOrder();
        datavec[iref].sol.resize(1);
		datavec[iref].sol[0].Resize(numdof);
		datavec[iref].sol.Fill(0.);
        datavec[iref].dsol.resize(1);
		datavec[iref].dsol[0].Redim(dim,numdof);
		datavec[iref].dsol[0].Zero();
		datavec[iref].axes.Redim(dim,3);
		datavec[iref].axes.Zero();
		
		trvec[iref].Apply(qsi, myqsi);
		msp->ComputeSolution(myqsi, datavec[iref].sol, datavec[iref].dsol, datavec[iref].axes);
		
		datavec[iref].x.Resize(3);
		msp->Reference()->X(myqsi, datavec[iref].x);
	}	
		int solSize = material->NSolutionVariables(var);
		sol.Resize(solSize);
		sol.Fill(0.);
		material->Solution(datavec, var, sol);
}

template <class TSHAPE, class TGeometry>
void TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix &axes)
{
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TSHAPE, class TGeometry>
void TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::ComputeSolution(TPZVec<REAL> &qsi,
					 TPZVec<REAL> &normal,
					 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix &leftaxes,
					 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix &rightaxes)
{
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TSHAPE, class TGeometry>
void TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
					 const TPZFMatrix &axes, TPZSolVec &sol, TPZGradSolVec &dsol)
{
	PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
	DebugStop();
}

template <class TSHAPE, class TGeometry>
void TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
	const int ncon = this->NConnects();
	int numeq = 0;
	int ic;

	for(ic=0; ic<ncon; ic++)
	{
		numeq += this->Connect(ic).NDof( *this->Mesh() );
	}

	int nstate = 0;

    int nref = this->fElementVec.size();

	for (int iref = 0; iref < nref; iref++) 
	{
		TPZCompEl *msp  = dynamic_cast < TPZCompEl* >(this->fElementVec[iref]);
        nstate += msp->Material()->NStateVariables();
	}

	const int numstate = nstate;
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,1);
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = numstate;
	ef.fNumStateVars = numstate;

	int i;
	for(i=0; i<ncon; i++){
        int ndof = this->Connect(i).NDof( *this->Mesh() );
#ifdef DEBUG
        TPZConnect &c = this->Connect(i);
        if (c.NShape()*c.NState() != ndof) {
            DebugStop();
        }
#endif
		ek.fBlock.Set(i,ndof);
		ef.fBlock.Set(i,ndof);
	}

	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ek.fConnect)[i] = this->ConnectIndex(i);
		(ef.fConnect)[i] = this->ConnectIndex(i);
	}
}//void

template <class TSHAPE, class TGeometry>
void TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::InitMaterialData(TPZVec<TPZMaterialData > &dataVec)
{
	this->Material()->FillDataRequirements(dataVec);
	const int dim = this->Dimension();
	int nref = this->fElementVec.size();
	
#ifdef DEBUG
	if (nref != dataVec.size())
    {
		PZError << "Error at " << __PRETTY_FUNCTION__ << " The number of materials can not be different from the size of the fElementVec !\n";
		DebugStop();
	}
#endif

	TPZVec<int> nshape(nref);
	for (int iref = 0; iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->fElementVec[iref]);

        const int nstate = msp->Material()->NStateVariables();
            nshape[iref] =  msp->NShapeF();

		dataVec[iref].phi.Redim(nshape[iref],1);
		dataVec[iref].dphix.Redim(dim,nshape[iref]);
		dataVec[iref].axes.Redim(dim,3);
		dataVec[iref].jacobian.Redim(dim,dim);
		dataVec[iref].jacinv.Redim(dim,dim);
		dataVec[iref].x.Resize(3);

		if (dataVec[iref].fNeedsSol)
		{
            dataVec[iref].sol.resize(1);
            dataVec[iref].dsol.resize(1);
			dataVec[iref].sol[0].Resize(nstate);
			dataVec[iref].dsol[0].Redim(dim,nstate);
		}
	}
}

template <class TSHAPE, class TGeometry>
void TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
	TPZAutoPointer<TPZMaterial> material = this->Material();

	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
		ek.Reset();
		ef.Reset();
		return;
	}

	InitializeElementMatrix(ek,ef);

	if (this->NConnects() == 0) return; // boundary discontinuous elements have this characteristic

	TPZVec<TPZMaterialData> datavec;
	const int nref = this->fElementVec.NElements();
	datavec.resize(nref);
	InitMaterialData(datavec);

	TPZManVector<TPZTransform> trvec;
	this->AffineTransform(trvec);

	int dim = this->Dimension();
	TPZAutoPointer<TPZIntPoints> intrule;

	TPZManVector<REAL,3> intpoint(dim, 0.), intpointtemp(dim, 0.);
	REAL weight = 0.;
	
	TPZVec<int> ordervec;
	ordervec.resize(nref);
	for (int iref = 0; iref < nref; iref++) 
	{
		TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->fElementVec[iref]);

        datavec[iref].p = msp->MaxOrder();
        ordervec[iref] = datavec[iref].p;
	}

	int order = material->IntegrationRuleOrder(ordervec);

	TPZGeoEl *ref = this->Reference();
	intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);

	TPZManVector<int,3> intorder(dim,order);
	intrule->SetOrder(intorder);
	int intrulepoints = intrule->NPoints();

	TPZFMatrix jac, axe, jacInv;
	REAL detJac;
	for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
	{
		intrule->Point(int_ind, intpointtemp, weight);
		ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
		weight *= fabs(detJac);

        // =============
        // Local problem: computation of new basis functions
        // =============

		TPZCompElMHM<TSHAPE, TGeometry> *mspMHM  = dynamic_cast < TPZCompElMHM<TSHAPE, TGeometry> *>(this->fElementVec[0]);
		trvec[0].Apply(intpointtemp, intpoint);

        mspMHM->ComputeShape(intpoint, datavec[0].x, datavec[0].jacobian, datavec[0].axes, datavec[0].detjac, datavec[0].jacinv, datavec[0].phi, datavec[0].dphix);

		datavec[0].intPtIndex = int_ind;

        mspMHM->ComputeRequiredData(datavec[0], intpoint);

        // ==============
        // Global problem: basis functions for one-dimensional problem.
        // ==============

        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->fElementVec[1]);
        trvec[1].Apply(intpointtemp, intpoint);

        msp->ComputeShape(intpoint, datavec[1].x, datavec[1].jacobian, datavec[1].axes,
                          datavec[1].detjac, datavec[1].jacinv, datavec[1].phi, datavec[1].dphix);

        datavec[1].intPtIndex = int_ind;

        msp->ComputeRequiredData(datavec[1], intpoint);

        // =============================
        // Computation of global problem
        // =============================

		material->Contribute(datavec, weight, ek.fMat, ef.fMat);
	}
}

#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "tpzgraphelt2dmapped.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel.h"
#include "pzmeshid.h"

template <class TSHAPE, class TGeometry>
void TPZMultiphysicsCompElMHM<TSHAPE, TGeometry>::EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv), TPZVec<REAL> & errors,TPZBlock * flux)
{
	int NErrors = this->Material()->NEvalErrors();
	errors.Resize(NErrors);
	errors.Fill(0.);

	TPZAutoPointer<TPZMaterial> material = this->Material();
	if(!material)
    {
		PZError << "TPZMultiphysicsCompEl::EvaluateError : no material for this element\n";
		TPZMultiphysicsCompEl<TGeometry>::Print(PZError);
		return;
	}

    TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->fElementVec[0]);

	int problemdimension = this->Mesh()->Dimension();
	if( msp->Reference()->Dimension() < problemdimension )
            return;

	// Adjust the order of the integration rule
	//Cesar 2007-06-27 ==>> Begin
	//this->MaxOrder is usefull to evaluate polynomial function of the aproximation space.
	//fp can be any function and max order of the integration rule could produce best results
	int dim = this->Dimension();
    int nref = this->fElementVec.size();

    TPZVec<TPZMaterialData> datavec;
	datavec.resize(nref);
	InitMaterialData(datavec);

	TPZAutoPointer<TPZIntPoints> intrule;
    int maxIntOrder = intrule->GetMaxOrder();

	TPZManVector<REAL,3> intpoint(dim,0.), intpointtemp(dim,0.);
	REAL weight = 0.;

	TPZVec<int> ordervec;
	ordervec.resize(nref);
	for (int iref=0;  iref<nref; iref++) 
	{
		TPZInterpolationSpace *msp_temp  = dynamic_cast <TPZInterpolationSpace *>(this->fElementVec[iref]);

        datavec[iref].p = msp_temp->MaxOrder();
        ordervec[iref] = datavec[iref].p;
	}
	int order = material->IntegrationRuleOrder(ordervec);
    
	TPZGeoEl *ref = this->Reference();
	intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);

	TPZManVector<int,3> intorder(dim,order);
	intrule->SetOrder(intorder);	
	int intrulepoints = intrule->NPoints();

	int ndof = material->NStateVariables();
	int nflux = material->NFluxes();
	TPZManVector<REAL,10> u_exact(ndof);
	TPZFNMatrix<90> du_exact(dim,ndof);
	TPZManVector<REAL,10> values(NErrors);
	values.Fill(0.0);
	TPZManVector<REAL,9> flux_el(nflux,0.);

    TPZManVector<TPZTransform> trvec;
	this->AffineTransform(trvec);

	for(int nint = 0; nint < intrulepoints; nint++)
    {
        //std::cout << "IntPoint: " << nint << "\n";

		intrule->Point(nint, intpoint, weight);
        trvec[0].Apply(intpointtemp, intpoint);

        msp->ComputeShape(intpoint, datavec[0].x, datavec[0].jacobian, datavec[0].axes, datavec[0].detjac, datavec[0].jacinv, datavec[0].phi, datavec[0].dphix);

		weight *= fabs(datavec[0].detjac);

		//contribuicoes dos erros
		if(fp)
        {
			fp(datavec[0].x,u_exact,du_exact);
			//std::cout << " funcao exata calculada no pto X " << datavec[0].x << " valor " << u_exact << " dudx " << du_exact << std::endl;

			//tentando implementar um erro meu
			if(datavec[0].fVecShapeIndex.NElements())
			{
                ComputeSolution(intpoint, datavec[0].phi, datavec[0].dphix, datavec[0].axes, datavec[0].sol, datavec[0].dsol);
				material->ErrorsHdiv(datavec[0], u_exact, du_exact, values);
			}
			else
            {
				ComputeSolution(intpoint, datavec[0].phi, datavec[0].dphix, datavec[0].axes, datavec[0].sol, datavec[0].dsol);
				material->Errors(datavec[0].x, datavec[0].sol[0], datavec[0].dsol[0], datavec[0].axes, flux_el, u_exact, du_exact, values);
			}

			//	std::cout<<"erro depois de Hdiv 2 "<<values<<std::endl;
			for(int ier = 0; ier < NErrors; ier++)
				errors[ier] += values[ier] * weight;
		}
        else
            std::cout << "fp is not defined... \n";
	}

	//Norma sobre o elemento
	for(int ier = 0; ier < NErrors; ier++)
		errors[ier] = sqrt(errors[ier]);

//    std::cout << "error: ";
//    errors.Print();
}

//---------------------------------------------------------------

template class TPZMultiphysicsCompElMHM<pzshape::TPZShapePoint, pzgeom::TPZGeoPoint>;
template class TPZMultiphysicsCompElMHM<pzshape::TPZShapeLinear, pzgeom::TPZGeoLinear>;
template class TPZMultiphysicsCompElMHM<pzshape::TPZShapeTriang, pzgeom::TPZGeoTriangle>;
template class TPZMultiphysicsCompElMHM<pzshape::TPZShapeQuad, pzgeom::TPZGeoQuad>;
template class TPZMultiphysicsCompElMHM<pzshape::TPZShapeCube, pzgeom::TPZGeoCube>;
template class TPZMultiphysicsCompElMHM<pzshape::TPZShapePrism, pzgeom::TPZGeoPrism>;
template class TPZMultiphysicsCompElMHM<pzshape::TPZShapeTetra, pzgeom::TPZGeoTetrahedra>;
template class TPZMultiphysicsCompElMHM<pzshape::TPZShapePiram, pzgeom::TPZGeoPyramid>;

TPZCompEl * CreateMultiphysicsMHMPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZMultiphysicsCompElMHM<pzshape::TPZShapePoint, pzgeom::TPZGeoPoint>(mesh, gel, index); 
}


TPZCompEl * CreateMultiphysicsMHMLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZMultiphysicsCompElMHM<pzshape::TPZShapeLinear, pzgeom::TPZGeoLinear>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsMHMTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZMultiphysicsCompElMHM<pzshape::TPZShapeTriang, pzgeom::TPZGeoTriangle>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsMHMQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZMultiphysicsCompElMHM<pzshape::TPZShapeQuad, pzgeom::TPZGeoQuad>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsMHMCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZMultiphysicsCompElMHM<pzshape::TPZShapeCube, pzgeom::TPZGeoCube>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsMHMPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZMultiphysicsCompElMHM<pzshape::TPZShapePrism, pzgeom::TPZGeoPrism>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsMHMTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZMultiphysicsCompElMHM<pzshape::TPZShapeTetra, pzgeom::TPZGeoTetrahedra>(mesh,gel,index);
}

TPZCompEl * CreateMultiphysicsMHMPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index)
{
	return new TPZMultiphysicsCompElMHM<pzshape::TPZShapePiram, pzgeom::TPZGeoPyramid>(mesh,gel,index);
}