/*
 *  ElasticMatInterface2D.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/23/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZElasticMatInterface2D.h"
#include "pzdiscgal.h"
#include "pzelasmat.h"
#include "pzlog.h"
#include "pzaxestools.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.poroelastic2d"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic.data"));
#endif

TPZElasticMatInterface2D::TPZElasticMatInterface2D() : TPZRegisterClassId(&TPZElasticMatInterface2D::ClassId),TPZElasticityMaterial(){
	fkn = 1000000.0;
	fkt = 1000000.0;		
}

TPZElasticMatInterface2D::TPZElasticMatInterface2D(int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress) : TPZRegisterClassId(&TPZElasticMatInterface2D::ClassId),TPZElasticityMaterial(num, E, nu, fx, fy, plainstress){
	fkn = 1000000.0;
	fkt = 1000000.0;	
}

TPZElasticMatInterface2D::~TPZElasticMatInterface2D()
	{
	}

// Contribute Interface implementations

void TPZElasticMatInterface2D::SetPenalty(REAL kn, REAL kt)
{
	fkn = kn;
	fkt = kt;
}

void TPZElasticMatInterface2D::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
//	PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	
//	Definition of penalty constants note: this constans for nolinear analysis are funtion of normal and tangencial forces.
	REAL kn = this->fkn;
	REAL kt = this->fkt;
	
//	TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
//	TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
//	TPZManVector<REAL,3> &normal = data.normal;

//	TPZFNMatrix<660> dphiL, dphiR;
//	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
//	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
	
//	int &LeftPOrder=dataleft.p;
//	int &RightPOrder=dataright.p;
	
//	REAL &faceSize=data.HSize;
	
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr;	
	
	
	for (int in = 0; in < 2*(nrowl + nrowr) ; in++) {
		ef(in,0) = 0.0;
	}
	REAL t[2] = {-data.normal[1],data.normal[0]};
	REAL n[2] = {data.normal[0],data.normal[1]};
	TPZFNMatrix<4,REAL> nx(2,2),tx(2,2);
	for (int i=0; i<2; i++) {
		for (int j=0; j<2 ; j++) {
			nx(i,j) = n[i]*n[j];
			tx(i,j) = t[i]*t[j];
		}
	}
	
	// Left Left contribution
	for(il=0; il < nrowl; il++) {
		for(jl=0; jl < nrowl; jl++) {
			ek(2*il,2*jl) += weight * (+ kn *phiL(il)*phiL(jl)*nx(0,0) + kt*phiL(il)*phiL(jl)*tx(0,0));
			ek(2*il+1,2*jl+1) += weight * (+ kn *phiL(il)*phiL(jl)*nx(1,1) + kt*phiL(il)*phiL(jl)*tx(1,1));
			ek(2*il+1,2*jl) += weight * (+ kn *phiL(il)*phiL(jl)*nx(1,0) + kt*phiL(il)*phiL(jl)*tx(1,0));
			ek(2*il,2*jl+1) += weight * (+ kn *phiL(il)*phiL(jl)*nx(0,1) + kt*phiL(il)*phiL(jl)*tx(0,1));
		}
	}
	// Right Left contribution	
	for(ir=0; ir < nrowr; ir++) {
		for(jl=0; jl < nrowl; jl++) {
			ek(2*ir+2*nrowl,2*jl) += weight * (- kn*phiR(ir)*phiL(jl)*nx(0,0) - kt*phiR(ir)*phiL(jl)*tx(0,0));
			ek(2*ir+1+2*nrowl,2*jl+1) += weight * (- kn*phiR(ir)*phiL(jl)*nx(1,1) - kt*phiR(ir)*phiL(jl)*tx(1,1));
			ek(2*ir+1+2*nrowl,2*jl) += weight * (- kn*phiR(ir)*phiL(jl)*nx(1,0) - kt*phiR(ir)*phiL(jl)*tx(1,0));
			ek(2*ir+2*nrowl,2*jl+1) += weight * (- kn*phiR(ir)*phiL(jl)*nx(0,1) - kt*phiR(ir)*phiL(jl)*tx(0,1));			
		}
	}
	
	// Left Right contribution		
	for(il=0; il < nrowl; il++) {
		for(jr=0; jr < nrowr; jr++) {
			ek(2*il,2*jr+2*nrowl) += weight * (- kn*phiR(jr)*phiL(il)*nx(0,0) - kt*phiR(jr)*phiL(il)*tx(0,0));
			ek(2*il+1,2*jr+1+2*nrowl) += weight * (- kn*phiR(jr)*phiL(il)*nx(1,1) - kt*phiR(jr)*phiL(il)*tx(1,1));
			ek(2*il+1,2*jr+2*nrowl) += weight * (- kn*phiR(jr)*phiL(il)*nx(1,0) - kt*phiR(jr)*phiL(il)*tx(1,0));
			ek(2*il,2*jr+1+2*nrowl) += weight * (- kn*phiR(jr)*phiL(il)*nx(0,1) - kt*phiR(jr)*phiL(il)*tx(0,1));				
		}
	}
	
	// Right Right contribution		
	for(ir=0; ir < nrowr; ir++) {
		for(jr=0; jr < nrowr; jr++) {
			ek(2*ir+2*nrowl,2*jr+2*nrowl) += weight * (+ kn *phiR(ir)*phiR(jr)*nx(0,0) + kt*phiR(ir)*phiR(jr)*tx(0,0));
			ek(2*ir+1+2*nrowl,2*jr+1+2*nrowl) += weight * (+ kn *phiR(ir)*phiR(jr)*nx(1,1) + kt*phiR(ir)*phiR(jr)*tx(1,1));
			ek(2*ir+1+2*nrowl,2*jr+2*nrowl) += weight * (+ kn *phiR(ir)*phiR(jr)*nx(1,0) + kt*phiR(ir)*phiR(jr)*tx(1,0));
			ek(2*ir+2*nrowl,2*jr+1+2*nrowl) += weight * (+ kn *phiR(ir)*phiR(jr)*nx(0,1) + kt*phiR(ir)*phiR(jr)*tx(0,1));	
		}
	}
	
}

void TPZElasticMatInterface2D::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	DebugStop();
}

void TPZElasticMatInterface2D::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef){
	PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	DebugStop();	
}

void TPZElasticMatInterface2D::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &left, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	DebugStop();	
}

int TPZElasticMatInterface2D::ClassId() {
	//CLASSIDFRANreturn TPZElasticityMaterial::ClassId()^Hash("TPZElasticMatInterface2D");
	return 666;
}
