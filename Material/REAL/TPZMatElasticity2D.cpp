//
//  TPZMatElasticity2D.cpp
//  PZ
//
//  Created by Omar on 10/27/14.
//
//


#include <iostream>
#include <string>
#include "TPZMatElasticity2D.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzlog.h"
#include <cmath>


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif


TPZMatElasticity2D::TPZMatElasticity2D():TPZMaterial()
{
    fE = 0.;
    fnu = 0.;
    flambda = 0.;
    fmu = 0.;
    ff.resize(3);
    ff[0]=0.;
    ff[1]=0.;
    ff[2]=0.;
    fPlaneStress = 1.;
    fPreStressXX = 0.0;
    fPreStressXY = 0.0;
    fPreStressYY = 0.0;
    fPreStressZZ = 0.0;
    
}

TPZMatElasticity2D::TPZMatElasticity2D(int matid):TPZMaterial(matid)
{
    fE = 0.;
    fnu = 0.;
    flambda = 0.;
    fmu = 0.;
    ff.resize(3);
    ff[0]=0.;
    ff[1]=0.;
    ff[2]=0.;
    fPlaneStress = 1.;
    fPreStressXX = 0.0;
    fPreStressXY = 0.0;
    fPreStressYY = 0.0;
    fPreStressZZ = 0.0;
}

TPZMatElasticity2D::TPZMatElasticity2D(int matid, REAL E, REAL nu, REAL fx, REAL fy, int plainstress):TPZMaterial(matid)
{
    fE = E;
    fnu = nu;
    flambda = (E*nu)/((1+nu)*(1-2*nu));
    fmu = E/(2*(1+nu));
    ff.resize(3);
    ff[0]=fx;
    ff[1]=fy;
    ff[2]=0.0;
    fPlaneStress = plainstress;
    fPreStressXX = 0.0;
    fPreStressXY = 0.0;
    fPreStressYY = 0.0;
    fPreStressZZ = 0.0;
}

TPZMatElasticity2D::~TPZMatElasticity2D()
{
}


TPZMatElasticity2D::TPZMatElasticity2D(const TPZMatElasticity2D &copy) : TPZMaterial(copy)
{
    fE = copy.fE;
    fnu = copy.fnu;
    flambda = copy.flambda;
    fmu = copy.fmu;
    ff.resize(copy.ff.size());
    for (int i = 0; i < copy.ff.size(); i++) {
        ff[i] = copy.ff[i];
    }
    fPlaneStress = copy.fPlaneStress;
    fPreStressXX = copy.fPreStressXX;
    fPreStressXY = copy.fPreStressXY;
    fPreStressYY = copy.fPreStressYY;
    fPreStressZZ = copy.fPreStressZZ;
}

TPZMatElasticity2D & TPZMatElasticity2D::operator=(const TPZMatElasticity2D &copy)
{
	TPZMaterial::operator = (copy);
    fE = copy.fE;
    fnu = copy.fnu;
    flambda = copy.flambda;
    fmu = copy.fmu;
    fPreStressXX = copy.fPreStressXX;
    fPreStressXY = copy.fPreStressXY;
    fPreStressYY = copy.fPreStressYY;
    fPreStressZZ = copy.fPreStressZZ;
    ff.resize(copy.ff.size());
    for (int i = 0; i < copy.ff.size(); i++) {
        ff[i] = copy.ff[i];
    }
    fPlaneStress = copy.fPlaneStress;
    return *this;
}

int TPZMatElasticity2D::NStateVariables() {
    return 2;
}


void TPZMatElasticity2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE>  &ek, TPZFMatrix<STATE> &ef) {
     
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiU     =  data.phi;    // vector of shapefunctions (format is dependent on the value of shapetype)
    TPZFMatrix<REAL> &dphiU     =  data.dphix;  // values of the derivative of the shape functions
    int phrU = phiU.Rows();
    int FirstU  = 0;
    
    REAL LambdaL, MuL, E, nu;
    
    // Functions computed at point x_{k} for each integration point
//    LambdaL     = flambda;
//    MuL         = fmu;
    
    //data.XCenter
    
    // Forcing Function -> Stochastic
    //TPZVec<REAL> res; // vector for E and nu
    TPZVec<STATE> res(2,0.0);
    
        if(ForcingFunction())
        {
            //this->fForcingFunction->Execute(data.x,res, fCorrelationMatrix);
            this->fForcingFunction->Execute(data.x,res);
            E = res[0];
            nu = fnu;
            
            LambdaL = (E*nu)/((1+nu)*(1-2*nu));
            MuL = E/(2*(1+nu));
            
            flambda=LambdaL;
            fmu = MuL;
            
        }
    
    
        //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
        //  Contribution of domain integrals for Jacobian matrix
        //  Elasticity Block (Equation for elasticity )
        //	Elastic equation
        //	Linear strain operator
        //	Ke Matrix
        // axes indicating the directions of the derivatives of the shapefunctions
    
        TPZFMatrix<REAL>	du(2,2);
        for(int iu = 0; iu < phrU; iu++ )
        {
            //	Derivative for Ux
            du(0,0) = dphiU(0,iu)*data.axes(0,0)+dphiU(1,iu)*data.axes(1,0); // du/dx
            //	Derivative for Uy
            du(1,0) = dphiU(0,iu)*data.axes(0,1)+dphiU(1,iu)*data.axes(1,1); // du/dy
            
            for(int ju = 0; ju < phrU; ju++)
            {
                //	Derivative for Vx
                du(0,1) = dphiU(0,ju)*data.axes(0,0)+dphiU(1,ju)*data.axes(1,0); // dv/dx
                //	Derivative for Vy
                du(1,1) = dphiU(0,ju)*data.axes(0,1)+dphiU(1,ju)*data.axes(1,1); // dv/dy
                
                if (this->fPlaneStress == 1)
                {
                    /* Plain stress state */
                    ek(2*iu + FirstU, 2*ju + FirstU)	     += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(0,0)*du(0,1)		+ (2*MuL)*du(1,0)*du(1,1));
                    
                    ek(2*iu + FirstU, 2*ju+1 + FirstU)       += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(0,0)*du(1,1)			+ (2*MuL)*du(1,0)*du(0,1));
                    
                    ek(2*iu+1 + FirstU, 2*ju + FirstU)       += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(1,0)*du(0,1)			+ (2*MuL)*du(0,0)*du(1,1));
                    
                    ek(2*iu+1 + FirstU, 2*ju+1 + FirstU)     += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(1,0)*du(1,1)		+ (2*MuL)*du(0,0)*du(0,1));
                }
                else
                {
                    /* Plain Strain State */
                    ek(2*iu + FirstU,2*ju + FirstU)         += weight*	((LambdaL + 2*MuL)*du(0,0)*du(0,1)	+ (MuL)*du(1,0)*du(1,1));
                    
                    ek(2*iu + FirstU,2*ju+1 + FirstU)       += weight*	(LambdaL*du(0,0)*du(1,1)			+ (MuL)*du(1,0)*du(0,1));
                    
                    ek(2*iu+1 + FirstU,2*ju + FirstU)       += weight*	(LambdaL*du(1,0)*du(0,1)			+ (MuL)*du(0,0)*du(1,1));
                    
                    ek(2*iu+1 + FirstU,2*ju+1 + FirstU)     += weight*	((LambdaL + 2*MuL)*du(1,0)*du(1,1)	+ (MuL)*du(0,0)*du(0,1));
                    
                }
            }
        }
        //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    this->Contribute(data,weight,ef);
    
    
}

void TPZMatElasticity2D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) {
    
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiU =  data.phi;
    TPZFMatrix<REAL> &dphiU = data.dphix;
    int phrU = phiU.Rows();
    int FirstU  = 0;
    
    TPZManVector<STATE,3> sol_u =data.sol[0]; /// vector of the solutions at the integration point
    TPZFNMatrix<4,STATE> dsol_u = data.dsol[0]; /// vector of the derivatives of the solution at the integration point
    
    
    //REAL LambdaL, MuL, E, nu;
    
    REAL LambdaL, MuL;
    
    // Functions computed at point x_{k} for each integration point
    LambdaL     = flambda;
    MuL         = fmu;
    
    
    //data.XCenter
//    
//    // Forcing Function -> Stochastic
//    //TPZVec<REAL> res; // vector for E and nu
//    TPZVec<STATE> res(2,0.0);
//    
//    if(ForcingFunction())
//    {
//        this->fForcingFunction->Execute(data.x,res);
//        E = res[0];
//        nu = res[1];
//        
//        LambdaL = (E*nu)/((1+nu)*(1-2*nu));
//        MuL = E/(2*(1+nu));
//        
//    }
//    
//    TPZVec<STATE> P(1,0.0);
//    TPZFMatrix<STATE> GradP(2,1,0.0);
    
//    if(this->HasffBCForcingFunction())
//    {
//        fForcingFunction->Execute(data.x,P,GradP);
////        REAL Pressure = P[0];
//    }
    
    
    //*************** fPreStress by Analytical Solution ****************//
    
    REAL SigmaX = 0., SigmaY = 0., SigmaXY = 0., SigmaZ = 0., SigmaXZ = 0., SigmaYZ = 0.;
    
    //data.dphix.Print(" Socorro ");
    
    if (fAnalytics == 1) {
        REAL theta = 0.;
        REAL coordY = 0.;
        REAL coordX = 0.;
        REAL r = 0.;
        coordX = data.x[0];
        coordY = data.x[1];
        theta = atan2(coordY,coordX);
        r = sqrt((coordX*coordX)+(coordY*coordY));
        
        AnalyticalWellboreSolution(SigmaX, SigmaY, SigmaXY, SigmaZ, SigmaXZ ,SigmaYZ ,theta, r);
    }
    else{
        SigmaX = fPreStressXX;
        SigmaY = fPreStressYY;
        SigmaXY = fPreStressXY;
        SigmaZ = fPreStressZZ;
    }
    
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  Contribution of domain integrals for Residual Vector
    //  Elastic equation
    //  Linear strain operator
    //  Ke Matrix
    
    TPZFNMatrix<4,REAL>    du(2,2);
    TPZFNMatrix<4,REAL>    dux(2,2,0.);
    TPZFNMatrix<4,REAL>    duy(2,2,0.);
    
    
    // Required check out of this implementation
    //  Derivative for Ux
    dux(0,1) = dsol_u(0,0)*data.axes(0,0)+dsol_u(1,0)*data.axes(1,0); // dUx/dx
    dux(1,1) = dsol_u(0,0)*data.axes(0,1)+dsol_u(1,0)*data.axes(1,1); // dUx/dy
    
    //  Derivative for Uy
    duy(0,1) = dsol_u(0,1)*data.axes(0,0)+dsol_u(1,1)*data.axes(1,0); // dUy/dx
    duy(1,1) = dsol_u(0,1)*data.axes(0,1)+dsol_u(1,1)*data.axes(1,1); // dUy/dy
    
    for(int iu = 0; iu < phrU; iu++ )
    {
        //  Derivative for Vx
        du(0,0) = dphiU(0,iu)*data.axes(0,0)+dphiU(1,iu)*data.axes(1,0); // dv/dx
        //  Derivative for Vy
        du(1,0) = dphiU(0,iu)*data.axes(0,1)+dphiU(1,iu)*data.axes(1,1); // dv/dy
        
//          Vector Force right hand term
             ef(2*iu + FirstU)     +=    weight*(ff[0]*phiU(iu, 0)- (du(0,0)* SigmaX + du(1,0)* SigmaXY));    // direcao x
             ef(2*iu+1 + FirstU)   +=    weight*(ff[1]*phiU(iu, 0)- (du(0,0)* SigmaXY + du(1,0)* SigmaY));    // direcao y

        if (fPlaneStress == 1)
        {
            /* Plain stress state */
            ef(2*iu + FirstU)           += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(0,0)*dux(0,1)     + (2*MuL)*du(1,0)*dux(1,1));
            
            ef(2*iu + FirstU)           += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(0,0)*duy(1,1)         + (2*MuL)*du(1,0)*duy(0,1));
            
            ef(2*iu+1 + FirstU)         += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(1,0)*dux(0,1)         + (2*MuL)*du(0,0)*dux(1,1));
            
            ef(2*iu+1 + FirstU)         += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(1,0)*duy(1,1)     + (2*MuL)*du(0,0)*duy(0,1));
        }
        else
        {
            /* Plain Strain State */
            ef(2*iu + FirstU)           += weight*  ((LambdaL + 2*MuL)*du(0,0)*dux(0,1)  + (MuL)*du(1,0)*(dux(1,1)));
            
            ef(2*iu + FirstU)           += weight*  (LambdaL*du(0,0)*duy(1,1)            + (MuL)*du(1,0)*(duy(0,1)));
            
            ef(2*iu+1 + FirstU)         += weight*  (LambdaL*du(1,0)*dux(0,1)            + (MuL)*du(0,0)*(dux(1,1)));
            
            ef(2*iu+1 + FirstU)         += weight*  ((LambdaL + 2*MuL)*du(1,0)*duy(1,1)  + (MuL)*du(0,0)*(duy(0,1)));
        }
    }
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    
}



void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{

    
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<REAL,3> sol_u = data.sol[0];
    TPZFMatrix<REAL> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    REAL uy = sol_u[1];
    
    int phru = phiu.Rows();
    short in,jn;
    TPZManVector<STATE,3> v2(3);
    TPZFMatrix<STATE> &v1 = bc.Val1();
    v2[0] = bc.Val2()(0,0);	//	Ux displacement or Tnx
    v2[1] = bc.Val2()(1,0);	//	Uy displacement or Tny
    
    //	Here each digit represent an individual boundary condition corresponding to each state variable.
    //	0 means Dirichlet condition on x-y
    //	1 means Neumann condition
    //	7 means Dirichlet condition on x
    //	8 means Dirichlet condition on y
    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    switch (bc.Type())
    {
        case 0 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)      += (BIGNUMBER*(ux - v2[0])*phiu(in,0))*weight;	// X displacement Value
                ef(2*in+1,0)	+= (BIGNUMBER*(uy - v2[1])*phiu(in,0))*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)       += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
            
        case 1 :
        {
            //	Neumann condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;		//	Tnx
                ef(2*in+1,0)	+= -1.0*v2[1]*phiu(in,0)*weight;		//	Tny
            }
            break;
        }
            
        case 2 :
        {
            //	Mixed condition for each state variable no used here
            //	Elasticity Equation
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += bc.Val1()(i,j)*data.sol[0][j];
            }
            
            for(in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += weight * ((v2[0]-res(0,0)) * phiu(in,0));
                ef(2*in+1,0) += weight * ((v2[1]-res(1,0)) * phiu(in,0));
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf< this->Dimension(); jdf++)
                    {
                        ek(2*in+idf,2*jn+jdf) += v1(idf,jdf)*phiu(in,0)*phiu(jn,0)*weight;
                        //      Not Complete with val2? HERE! PHIL!!!!
                        //      DebugStop();
                    }
                }
            }
            
            break;
        }
            
        case 3 :
        {
            //	Null Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)      += (BIGNUMBER*(0.0 - v2[0])*phiu(in,0))*weight;	// X displacement Value
                ef(2*in+1,0)	+= (BIGNUMBER*(0.0 - v2[1])*phiu(in,0))*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)       += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                    ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
            
        case 4 :
        {
            // Stress field
            
            // BC as a function of the Analytic Solution
            if (fAnalytics == 1||2) {
                REAL theta = 0.;
                REAL coordY = 0.;
                REAL coordX = 0.;
                REAL r = 0.;
                coordX = data.x[0];
                coordY = data.x[1];
                theta = atan2(coordY,coordX);
                r = sqrt((coordX*coordX)+(coordY*coordY));
                
                REAL Sx=0., Sy=0., Sxy=0., Sz=0., Sxz=0., Syz=0.;
                
                AnalyticalWellboreSolution(Sx, Sy, Sxy, Sz, Sxz, Syz, theta, r);
                
                v1(0,0) = Sx;
                v1(0,1) = Sxy;
                v1(1,0) = Sxy;
                v1(1,1) = Sy;
                
                //v1.Print(" Valor de v1 ");
                
            }
            
            for(in = 0; in < this->Dimension(); in ++){
                v2[in] = ( v1(in,0) * data.normal[0] + v1(in,1) * data.normal[1]);
            }
            
            for(in = 0 ; in <phru; in++)
            {
                //	Normal Tension Components on neumann boundary
                ef(2*in,0)      += 1.0*v2[0]*phiu(in,0)*weight;      //	Tnx
                ef(2*in+1,0)	+= 1.0*v2[1]*phiu(in,0)*weight;		//	Tny
            }
            
            break;
        }
            
        case 5 :
        {
            DebugStop();
        }
        break;
            
            
        case 6 :
            //	Normal Pressure condition Pressure value Should be inserted in v2[0]
            //	Elasticity Equation
            {
                
                TPZManVector<REAL> n = data.normal;
//                TPZManVector<REAL> n_ab = data.normal;
//                
//                REAL lxx = 0., lxy = 0., lxz = 0., lyx =0., lyy = 0., lyz = 0., lzx = 0., lzy = 0., lzz = 0.;
//                
//                // x-diretion
//                lxx = cos(falpha)*cos(fbeta);
//                lxy = sin(falpha)*cos(fbeta);
//                lxz = -sin(fbeta);
//                // y-direction
//                lyx = -sin(falpha);
//                lyy = cos(falpha);
//                lyz = 0;
//                // z-direction
//                lzx = cos(falpha)*sin(fbeta);
//                lzy = sin(falpha)*sin(fbeta);
//                lzz = cos(fbeta);
//                
//                n_ab[0] = lxx*n[0] + lxy*n[1] + lxz*n[2];
//                n_ab[1] = lyx*n[0] + lyy*n[1] + lyz*n[2];
//                n_ab[2] = lzx*n[0] + lzy*n[1] + lzz*n[2];
                
                TPZFNMatrix<2,STATE> Tn(2,1,0.);
                for(int i=0; i<2; i++)
                {
                    for(int j=0; j<2; j++)
                    {
                        Tn(i,0) += bc.Val1()(i,j)*n[j];
                    }
                }
                
                for(int in = 0 ; in < phru; in++)
                {
                    ef(2*in+0,0) += weight * Tn(0,0)* phiu(in,0);
                    ef(2*in+1,0) += weight * Tn(1,0) * phiu(in,0);
                }
            }
            break;
            
            
        case 7 :
        {
            //	Dirichlet condition for each state variable
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in,0)		+= (BIGNUMBER*(ux - v2[0])*phiu(in,0))*weight;	// X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
                }
            }
            
            break;
        }
            
        case 8 :
        {
            //	Dirichlet condition for uy
            //	Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //	Contribution for load Vector
                ef(2*in+1,0)	+= (BIGNUMBER*(uy - v2[1])*phiu(in,0))*weight;	// y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
                }
            }
            
            break;
        }
            
        default:
        {
            PZError << "TPZMatElasticity2D::ContributeBC error - Wrong boundary condition type" << std::endl;
            DebugStop();
        }
            break;
    }

}




void TPZMatElasticity2D::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNeighborSol = true;
    data.fNeedsNeighborCenter = false;
    data.fNeedsNormal = true;
}

void TPZMatElasticity2D::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}


void TPZMatElasticity2D::Print(std::ostream &out)
{
    out << "Material Name : " << Name() << "\n";
    out << "Plane Problem (fPlaneStress = 0, for Plane Strain conditions) " << fPlaneStress << std::endl;
    out << "Properties for elasticity: \n";
    out << "\t Young modulus   = "											<< fE		<< std::endl;
    out << "\t Poisson Ratio   = "											<< fnu		<< std::endl;
    out << "\t First Lamé Parameter   = "									<< flambda	<< std::endl;
    out << "\t Second Lamé Parameter   = "									<< fmu		<< std::endl;
    out << "\t Body force vector B {X-direction, Y-direction}   = "			<< ff[0] << ' ' << ff[1]   << std::endl;
    out << "\t fPreStressXX   = "                                           << fPreStressXX << std::endl;
    out << "\t fPreStressXY   = "                                           << fPreStressXY << std::endl;
    out << "\t fPreStressYY   = "                                           << fPreStressYY << std::endl;
    out << "\t fPreStressZZ   = "                                           << fPreStressZZ << std::endl;
    out << "Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
    
}

/** Returns the variable index associated with the name */
int TPZMatElasticity2D::VariableIndex(const std::string &name)
{
    //	Elasticity Variables
    if(!strcmp("Displacement",name.c_str()))			  	    return	1;
    if(!strcmp("SolidPressure",name.c_str()))				    return	2;
    if(!strcmp("SigmaX",name.c_str()))						    return	3;
    if(!strcmp("SigmaY",name.c_str()))						    return	4;
    if(!strcmp("SigmaZ",name.c_str()))						    return	5;
    if(!strcmp("TauXY",name.c_str()))						    return	6;
    if(!strcmp("SigmaXAnalytic",name.c_str()))				    return	7;
    if(!strcmp("SigmaYAnalytic",name.c_str()))				    return	8;
    if(!strcmp("SigmaZAnalytic",name.c_str()))				    return	9;
    if(!strcmp("TauXYAnalytic",name.c_str()))				    return	10;
    if(!strcmp("SolidPressureAnalytic",name.c_str()))		    return	11;
    if(!strcmp("SigmaXProjected",name.c_str()))		            return	12;
    if(!strcmp("SigmaYProjected",name.c_str()))				    return	13;
    if(!strcmp("SigmaZProjected",name.c_str()))				    return	14;
    if(!strcmp("TauXYProjected",name.c_str()))				    return	15;
    if(!strcmp("SolidPressureProjected",name.c_str()))		    return	16;
    if(!strcmp("SigmaXAnalyticProjected",name.c_str()))		    return	17;
    if(!strcmp("SigmaYAnalyticProjected",name.c_str()))		    return	18;
    if(!strcmp("SigmaZAnalyticProjected",name.c_str()))		    return	19;
    if(!strcmp("TauXYAnalyticProjected",name.c_str()))		    return	20;
    if(!strcmp("SolidPressureAnalyticProjected",name.c_str()))	return	21;
    if(!strcmp("ExxAnalytic",name.c_str()))                     return	22;
    if(!strcmp("EyyAnalytic",name.c_str()))                     return	23;
    if(!strcmp("ExyAnalytic",name.c_str()))                     return	24;
    if(!strcmp("Exx",name.c_str()))                             return	25;
    if(!strcmp("Eyy",name.c_str()))                             return	26;
    if(!strcmp("Exy",name.c_str()))                             return	27;
    if(!strcmp("J2",name.c_str()))                              return	28;
    if(!strcmp("F1",name.c_str()))                              return	29;
    if(!strcmp("I1",name.c_str()))                   return	30;
    if(!strcmp("Sigma1",name.c_str()))                          return	31;
    if(!strcmp("Sigma2",name.c_str()))                          return	32;
    if(!strcmp("Sigma3",name.c_str()))                          return	33;
    if(!strcmp("CheckingVM1",name.c_str()))                     return	34;
    if(!strcmp("CheckingVM2",name.c_str()))                     return	35;
    if(!strcmp("CheckingVM3",name.c_str()))                     return	36;
    if(!strcmp("F_Mogi-Coulomb",name.c_str()))                  return	37;
    if(!strcmp("J2_Projected",name.c_str()))                     return	38;
    if(!strcmp("F1_Projected",name.c_str()))                     return	39;
    if(!strcmp("I1_Projected",name.c_str()))                  return	40;
    
    PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
    return -1;
    
    return TPZMaterial::VariableIndex(name);
}

/**
 * Save the element data to a stream
 */
void TPZMatElasticity2D::Write(TPZStream &buf, int withclassid)
{
    TPZMaterial::Write(buf,withclassid);
    buf.Write(&fE);
    buf.Write(&fnu);
    buf.Write(&flambda);
    buf.Write(&fmu);
    TPZSaveable::WriteObjects(buf, ff);
    buf.Write(&fPreStressXX);
    buf.Write(&fPreStressXY);
    buf.Write(&fPreStressYY);
    buf.Write(&fPreStressZZ);
    buf.Write(&fPlaneStress);
    
}

/**
 * Read the element data from a stream
 */
void TPZMatElasticity2D::Read(TPZStream &buf, void *context)
{
    TPZMaterial::Read(buf,context);
    buf.Read(&fE);
    buf.Read(&fnu);
    buf.Read(&flambda);
    buf.Read(&fmu);
    TPZSaveable::ReadObjects(buf, ff);
    buf.Read(&fPreStressXX);
    buf.Read(&fPreStressXY);
    buf.Read(&fPreStressYY);
    buf.Read(&fPreStressZZ);
    buf.Read(&fPlaneStress);
    
}

int TPZMatElasticity2D::NSolutionVariables(int var){
    if(var == 1)	return 3;
    if(var == 2)	return 1;
    if(var == 3)	return 1;
    if(var == 4)	return 1;
    if(var == 5)	return 1;
    if(var == 6)	return 1;
    if(var == 7)	return 1;
    if(var == 8)	return 1;
    if(var == 9)	return 1;
    if(var == 10)	return 1;
    if(var == 11)	return 1;
    if(var == 12)	return 1;
    if(var == 13)	return 1;
    if(var == 14)	return 1;
    if(var == 15)	return 1;
    if(var == 16)	return 1;
    if(var == 17)	return 1;
    if(var == 18)	return 1;
    if(var == 19)	return 1;
    if(var == 20)	return 1;
    if(var == 21)	return 1;
    if(var == 22)	return 1;
    if(var == 23)	return 1;
    if(var == 24)	return 1;
    if(var == 25)	return 1;
    if(var == 26)	return 1;
    if(var == 27)	return 1;
    if(var == 28)	return 1;
    if(var == 29)	return 1;
    if(var == 30)	return 1;
    if(var == 31)	return 1;
    if(var == 32)	return 1;
    if(var == 33)	return 1;
    if(var == 34)	return 1;
    if(var == 35)	return 1;
    if(var == 36)	return 1;
    if(var == 37)	return 1;
    if(var == 38)	return 1;
    if(var == 39)	return 1;
    if(var == 40)	return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

//	Calculate Secondary variables based on ux, uy, Pore pressure and their derivatives
void TPZMatElasticity2D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize(this->NSolutionVariables(var));
    
    TPZManVector<STATE,3> SolU, SolP;
    TPZFNMatrix <6,STATE> DSolU, DSolP;
    TPZFNMatrix <9> axesU, axesP;
    
    TPZVec<REAL> ptx(3);
    TPZVec<STATE> solExata(3);
    TPZFMatrix<STATE> flux(5,1);
    
    if (data.sol.size() != 1) {
        DebugStop();
    }
    
    STATE x = data.x[0];
    STATE y = data.x[1];
    
    SolU	=	data.sol[0];
    DSolU	=	data.dsol[0];
    axesU	=	data.axes;
    
    
    //	Displacements
    if(var == 1){
        Solout[0] = SolU[0];
        Solout[1] = SolU[1];
        Solout[2] = 0.0;
        return;
    }
    
    
    //***** Analytic PreStresses *******//
    REAL SigmaX = 0., SigmaY = 0., SigmaXY = 0., SigmaZ = 0., SigmaXZ = 0., SigmaYZ = 0.;

    if (fAnalytics == 1) { //E fAnalytics == 2???
        REAL theta = 0.;
        REAL coordY = 0.;
        REAL coordX = 0.;
        REAL r = 0.;
        coordX = data.x[0];
        coordY = data.x[1];
        theta = atan2(coordY,coordX);
        r = sqrt((coordX*coordX)+(coordY*coordY));
        
        AnalyticalWellboreSolution(SigmaX, SigmaY, SigmaXY, SigmaZ, SigmaXZ, SigmaYZ, theta, r);
    }
    
    else {
        
        SigmaX = fPreStressXX;
        SigmaY = fPreStressYY;
        SigmaXY = fPreStressXY;
        SigmaZ = fPreStressZZ;
    }
    
    
//    std::cout<< SigmaX << std::endl;
//    std::cout<< SigmaY << std::endl;
//    std::cout<< SigmaXY << std::endl;
//    std::cout<< SigmaZ << std::endl;

    
    
    //************ Calculates Variation ***********//
    
    REAL epsx;
    REAL epsy;
    REAL epsxy;
    REAL SigX;
    REAL SigY;
    REAL SigZ;
    REAL Tau, DSolxy[2][2];
    REAL divu;
    
    DSolxy[0][0] = DSolU(0,0)*axesU(0,0)+DSolU(1,0)*axesU(1,0); // dUx/dx
    DSolxy[1][0] = DSolU(0,0)*axesU(0,1)+DSolU(1,0)*axesU(1,1); // dUx/dy
    
    DSolxy[0][1] = DSolU(0,1)*axesU(0,0)+DSolU(1,1)*axesU(1,0); // dUy/dx
    DSolxy[1][1] = DSolU(0,1)*axesU(0,1)+DSolU(1,1)*axesU(1,1); // dUy/dy
    
    
    epsx = DSolxy[0][0];// du/dx
    epsy = DSolxy[1][1];// dv/dy
    epsxy = 0.5*(DSolxy[1][0]+DSolxy[0][1]);
    REAL C11 = 4.0*(fmu)*(flambda+fmu)/(flambda+2.0*fmu);
    REAL C22 = 2.0*(fmu)*(flambda)/(flambda+2.0*fmu);
    
    if (this->fPlaneStress)
    {
        SigX = C11*epsx+C22*epsy;
        SigY = C11*epsy+C22*epsx;
        SigZ = 0.0;
        Tau = 2.0*fmu*epsxy;
    }
    else
    {
        SigX = ((flambda + 2.0*fmu)*(epsx) + (flambda)*epsy);
        SigY = ((flambda + 2.0*fmu)*(epsy) + (flambda)*epsx);
        SigZ = fnu*(SigX+SigY);
        Tau = 2.0*fmu*epsxy;		
    }
    
    
//    std::cout<< SigX << std::endl;
//    std::cout<< SigY << std::endl;
//    std::cout<< Tau << std::endl;
//    std::cout<< SigZ << std::endl;

    
    
    //********** Stresses Solution ***********//
    
    //	Hydrostatic stress
    if(var == 2) 
    {
        Solout[0] = ((SigX + SigmaX)+(SigY + SigmaY)+(SigZ+SigmaZ))/3.0;
        return;
    }
    
    //	Effective Stress x-direction
    if(var == 3) {
        Solout[0] = SigX + SigmaX;
        return;
    }
    
    //	Effective Stress y-direction	
    if(var == 4) {
        Solout[0] = SigY + SigmaY;
        return;
    }
    
    //	Effective Stress y-direction
    if(var == 5) {
        Solout[0] = SigZ + SigmaZ;
        return;
    }
    
    //	Shear Stress	
    if(var == 6) {
        Solout[0] = Tau + SigmaXY;
        return;
    }
    
    
    
    /******** Analytical Solution ********/
    
    //***** Analytic Stresses *******//
    REAL SigmaAnX = 0., SigmaAnY = 0., SigmaAnXY = 0., SigmaAnZ = 0., SigmaAnXZ = 0., SigmaAnYZ = 0.;
    
    REAL theta = 0.;
    REAL coordY = 0.;
    REAL coordX = 0.;
    REAL r = 0.;
    coordX = data.x[0];
    coordY = data.x[1];
    theta = atan2(coordY,coordX);
    r = sqrt((coordX*coordX)+(coordY*coordY));
    
            AnalyticalWellboreSolution(SigmaAnX, SigmaAnY, SigmaAnXY, SigmaAnZ, SigmaAnXZ, SigmaAnYZ, theta, r);
    
    
    //	Analytical Stress x-direction
    if(var == 7) {
        Solout[0] = SigmaAnX;
        return;
    }
    
    //	Analytical Stress y-direction
    if(var == 8) {
        Solout[0] = SigmaAnY;
        return;
    }
    
    // Analytical Stress z-direction
    if(var == 9) {
        Solout[0] = SigmaAnZ;
        return;
    }
    
    
    //	Analytical Shear Stress
    if(var == 10) {
        Solout[0] = SigmaAnXY;
        return;
    }
    
    
    //	Hydrostatic Analytical stress
    if(var == 11)
    {
        Solout[0] = ((SigmaAnX)+(SigmaAnY)+(SigmaAnZ))/3.0;
        return;
    }
   
    
    
    
    /******** Projected Solution ********/
    //******* Cod feito com rotacao inversa ***********//
    
    REAL SigXP = 0., SigYP = 0., TauP = 0., SigZP = 0.;
    
    // *********** What????? *******************
    //******** SigXZ e SigYZ nao sao calculados!!! **********
    REAL SigXZ = SigmaAnXZ; //*** Modificado: considerando Sol da Eq Analitica ***//
    REAL SigYZ = SigmaAnYZ; //*** Modificado: considerando Sol da Eq Analitica ***//
    
    // Cria variavel para solucao FEM
    REAL SolX = 0., SolY = 0., SolZ = 0., SolTau = 0.;
    
    SolX =    SigX+SigmaX;
    SolY =    SigY+SigmaY;
    SolZ =    SigZ+SigmaZ;
    SolTau =  Tau+SigmaXY;
    
    
    //FEM Sol: Sig
    //PreStress: Sigma
    
    SigXP  = (SolY)*pow(sin(falpha),2) - sin(2.0*falpha)*((SolTau)*cos(fbeta) + SigYZ*sin(fbeta)) + pow(cos(falpha),2)*((SolX)*pow(cos(fbeta),2) + (SolZ)*pow(sin(fbeta),2) + SigXZ*sin(2.0*fbeta));
       //Sigy*Power(Sin(Degree*\[Alpha]),2) - Sin(2*Degree*\[Alpha])*(Sigxy*Cos(Degree*\[Beta]) + Sigyz*Sin(Degree*\[Beta])) + Power(Cos(Degree*\[Alpha]),2)*(Sigx*Power(Cos(Degree*\[Beta]),2) + Sigz*Power(Sin(Degree*\[Beta]),2) + Sigxz*Sin(2*Degree*\[Beta]))
    
    SigYP  = (SolY)*pow(cos(falpha),2) + (SolX)*pow(cos(fbeta),2)*pow(sin(falpha),2) + (SolTau)*cos(fbeta)*sin(2.0*falpha) + SigYZ*sin(2.0*falpha)*sin(fbeta) + (SolZ)*pow(sin(falpha),2)*pow(sin(fbeta),2) + SigXZ*pow(sin(falpha),2)*sin(2.0*fbeta);
       //Sigy*Power(Cos(Degree*\[Alpha]),2) + Sigx*Power(Cos(Degree*\[Beta]),2)*Power(Sin(Degree*\[Alpha]),2) + Sigxy*Cos(Degree*\[Beta])*Sin(2*Degree*\[Alpha]) + Sigyz*Sin(2*Degree*\[Alpha])*Sin(Degree*\[Beta]) + Sigz*Power(Sin(Degree*\[Alpha]),2)*Power(Sin(Degree*\[Beta]),2) + Sigxz*Power(Sin(Degree*\[Alpha]),2)*Sin(2*Degree*\[Beta])
    
    SigZP  = (SolZ)*pow(cos(fbeta),2) + (SolX)*pow(sin(fbeta),2) - SigXZ*sin(2.0*fbeta);
        //Sigz*Power(Cos(Degree*\[Beta]),2) + Sigx*Power(Sin(Degree*\[Beta]),2) - Sigxz*Sin(2*Degree*\[Beta])
    
    TauP   = pow(cos(falpha),2)*((SolTau)*cos(fbeta) + SigYZ*sin(fbeta)) - pow(sin(falpha),2)*((SolTau)*cos(fbeta) + SigYZ*sin(fbeta)) + (sin(2.0*falpha)*((SolX) - 2.0*(SolY) + (SolZ) + ((SolX) - (SolZ))*cos(2.0*fbeta) + 2.0*SigXZ*sin(2.0*fbeta)))/4.;
        //Power(Cos(Degree*\[Alpha]),2)*(Sigxy*Cos(Degree*\[Beta]) + Sigyz*Sin(Degree*\[Beta])) - Power(Sin(Degree*\[Alpha]),2)*(Sigxy*Cos(Degree*\[Beta]) + Sigyz*Sin(Degree*\[Beta])) +     (Sin(2*Degree*\[Alpha])*(Sigx - 2*Sigy + Sigz + (Sigx - Sigz)*Cos(2*Degree*\[Beta]) + 2*Sigxz*Sin(2*Degree*\[Beta])))/4.
    
    
    //	Projected Stress x-direction
    if(var == 12) {
        Solout[0] = SigXP;
        return;
    }
    
    //	Projected Stress y-direction
    if(var == 13) {
        Solout[0] = SigYP;
        return;
    }
    
    // Projected Stress z-direction
    if(var == 14) {
        Solout[0] = SigZP;
        return;
    }
    
    
    //	Projected Shear Stress
    if(var == 15) {
        Solout[0] = TauP;
        return;
    }
    
    
    //	Projected Hydrostatic Stress
    if(var == 16)
    {
        Solout[0] = ((SigXP)+(SigYP)+(SigZP))/3.0;
        return;
    }
    
    
    /******** Analytical Projected Solution ********/
    //*******  Cod feito com rotacao inversa ***********//
    
    REAL SigmaXP = 0., SigmaYP = 0., SigmaXYP = 0., SigmaZP = 0., SigmaXZP=0., SigmaYZP = 0.;
    
    // *********** What????? *******************
    //******** SigXZ e SigYZ nao sao calculados!!! **********
    
    //FEM Sol: Sig
    //Analytic Stress: SigmaAn
    
    
    SigmaXP  = SigmaAnY*pow(sin(falpha),2) - sin(2.0*falpha)*(SigmaAnXY*cos(fbeta) + SigmaAnYZ*sin(fbeta)) + pow(cos(falpha),2)*(SigmaAnX*pow(cos(fbeta),2) + SigmaAnZ*pow(sin(fbeta),2) + SigmaAnXZ*sin(2.0*fbeta));
    //Sigy*Power(Sin(Degree*\[Alpha]),2) - Sin(2*Degree*\[Alpha])*(Sigxy*Cos(Degree*\[Beta]) + Sigyz*Sin(Degree*\[Beta])) + Power(Cos(Degree*\[Alpha]),2)*(Sigx*Power(Cos(Degree*\[Beta]),2) + Sigz*Power(Sin(Degree*\[Beta]),2) + Sigxz*Sin(2*Degree*\[Beta]))
    
    SigmaYP  = SigmaAnY*pow(cos(falpha),2) + SigmaAnX*pow(cos(fbeta),2)*pow(sin(falpha),2) + SigmaAnXY*cos(fbeta)*sin(2.0*falpha) + SigmaAnYZ*sin(2.0*falpha)*sin(fbeta) + SigmaAnZ*pow(sin(falpha),2)*pow(sin(fbeta),2) + SigmaAnXZ*pow(sin(falpha),2)*sin(2.0*fbeta);
    //Sigy*Power(Cos(Degree*\[Alpha]),2) + Sigx*Power(Cos(Degree*\[Beta]),2)*Power(Sin(Degree*\[Alpha]),2) + Sigxy*Cos(Degree*\[Beta])*Sin(2*Degree*\[Alpha]) + Sigyz*Sin(2*Degree*\[Alpha])*Sin(Degree*\[Beta]) + Sigz*Power(Sin(Degree*\[Alpha]),2)*Power(Sin(Degree*\[Beta]),2) + Sigxz*Power(Sin(Degree*\[Alpha]),2)*Sin(2*Degree*\[Beta])
    
    SigmaZP  = SigmaAnZ*pow(cos(fbeta),2) + SigmaAnX*pow(sin(fbeta),2) - SigmaAnXZ*sin(2.0*fbeta);
    //Sigz*Power(Cos(Degree*\[Beta]),2) + Sigx*Power(Sin(Degree*\[Beta]),2) - Sigxz*Sin(2*Degree*\[Beta])
    
    SigmaXYP   = pow(cos(falpha),2)*(SigmaAnXY*cos(fbeta) + SigmaAnYZ*sin(fbeta)) - pow(sin(falpha),2)*(SigmaAnXY*cos(fbeta) + SigmaAnYZ*sin(fbeta)) + (sin(2.0*falpha)*(SigmaAnX - 2.0*SigmaAnY + SigmaAnZ + (SigmaAnX - SigmaAnZ)*cos(2.0*fbeta) + 2.0*SigmaAnXZ*sin(2.0*fbeta)))/4.;
    //Power(Cos(Degree*\[Alpha]),2)*(Sigxy*Cos(Degree*\[Beta]) + Sigyz*Sin(Degree*\[Beta])) - Power(Sin(Degree*\[Alpha]),2)*(Sigxy*Cos(Degree*\[Beta]) + Sigyz*Sin(Degree*\[Beta])) +     (Sin(2*Degree*\[Alpha])*(Sigx - 2*Sigy + Sigz + (Sigx - Sigz)*Cos(2*Degree*\[Beta]) + 2*Sigxz*Sin(2*Degree*\[Beta])))/4.
    
    SigmaXZP = sin(falpha)*(-(SigmaAnYZ*cos(fbeta)) + SigmaAnXY*sin(fbeta)) + (cos(falpha)*(2.0*SigmaAnXZ*cos(2.0*fbeta) + (-SigmaAnX + SigmaAnZ)*sin(2.0*fbeta)))/2.;
    //Sin(Degree*\[Alpha])*(-(Sigyz*Cos(Degree*\[Beta])) + Sigxy*Sin(Degree*\[Beta])) + (Cos(Degree*\[Alpha])*(2*Sigxz*Cos(2*Degree*\[Beta]) + (-Sigx + Sigz)*Sin(2*Degree*\[Beta])))/2.
    
    
    SigmaYZP = cos(falpha)*(SigmaAnYZ*cos(fbeta) - SigmaAnXY*sin(fbeta)) + (sin(falpha)*(2.0*SigmaAnXZ*cos(2.0*fbeta) + (-SigmaAnX + SigmaAnZ)*sin(2.0*fbeta)))/2.;
    //Cos(Degree*\[Alpha])*(Sigyz*Cos(Degree*\[Beta]) - Sigxy*Sin(Degree*\[Beta])) + (Sin(Degree*\[Alpha])*(2*Sigxz*Cos(2*Degree*\[Beta]) + (-Sigx + Sigz)*Sin(2*Degree*\[Beta])))/2.
    
    
    //	Projected Stress x-direction
    if(var == 17) {
        Solout[0] = SigmaXP;
        return;
    }
    
    //	Projected Stress y-direction
    if(var == 18) {
        Solout[0] = SigmaYP;
        return;
    }
    
    // Projected Stress z-direction
    if(var == 19) {
        Solout[0] = SigmaZP;
        return;
    }
    
    
    //	Projected Shear Stress
    if(var == 20) {
        Solout[0] = SigmaXYP;
        return;
    }
    
    
    //	Projected Hydrostatic Stress
    if(var == 21)
    {
        Solout[0] = ((SigmaXP)+(SigmaYP)+(SigmaZP))/3.0;
        return;
    }
    
    
    
    /******** Analytical Deformation ********/
    
    REAL exx = 0., eyy = 0., exy = 0.;
    
    exx= (SigmaAnX-fnu*(SigmaAnY+SigmaAnZ))/fE;
    //(\[Sigma]xx - \[Nu]*(\[Sigma]yy + \[Sigma]zz))/El
    
    eyy= (SigmaAnY-fnu*(SigmaAnX+SigmaAnZ))/fE;
    //(\[Sigma]yy - \[Nu]*(\[Sigma]xx + \[Sigma]zz))/El
    
    exy = ((1.0+fnu)*SigmaAnXY)/fE;
    //((1 + \[Nu])*\[Sigma]xy)/El
    
    
    //	exx analytical
    if(var == 22)
    {
        Solout[0] = exx;
        return;
    }
    
    //	eyy analytical
    if(var == 23)
    {
        Solout[0] = eyy;
        return;
    }
    
        //	exy analytical
    if(var == 24)
    {
        Solout[0] = exy;
        return;
    }
    
    
    
    /******** Numerical Deformation ********/
    
    //	exx Numerical
    if(var == 25)
    {
        Solout[0] = epsx;
        return;
    }
    
    //	eyy Numerical
    if(var == 26)
    {
        Solout[0] = epsy;
        return;
    }
    
    //	exy Numerical
    if(var == 27)
    {
        Solout[0] = epsxy;
        return;
    }

    
    
    
    ///////////////************************************************ ELASTOPLASTICITY ****************************************************///////////////
    
    
    //Stress Tensor
    TPZFNMatrix<9> T(3,3,0.);
  
    
    //T = STRESS TENSOR
    T.PutVal(0,0, (SigX+SigmaX));
    T.PutVal(0,1, (Tau+SigmaXY));
    T.PutVal(0,2, SigmaAnXZ);  /// ZERO MESMO????
    T.PutVal(1,0, (Tau+SigmaXY));
    T.PutVal(1,1, (SigY+SigmaY));
    T.PutVal(1,2, SigmaAnYZ); /// ZERO MESMO????
    T.PutVal(2,0, SigmaAnXZ); /// ZERO MESMO????
    T.PutVal(2,1, SigmaAnYZ); /// ZERO MESMO????
    T.PutVal(2,2, (SigZ+SigmaZ));
    
    //    std::cout << T << std::endl;
    
    
    long NumIt = 1000;
    REAL tol = 1.E-5;
    TPZVec<REAL> EigValues(3,0.);
    TPZFNMatrix<9,REAL> EigVectors(3,3,0.);
    bool EigenWorks;
    EigenWorks = T.SolveEigensystemJacobi(NumIt, tol, EigValues, EigVectors);
    
   // std::cout << EigValues << std::endl;
    
    REAL Sigma1 = 0., Sigma2 = 0., Sigma3 = 0.;
    
    //**********// Criar metodo para garantir que Sig1 > Sig2 > Sig3 (em modulo)
        REAL temp;
    for (int i = 0; i < EigValues.size() - 1; i++) {
        for (int j = 1; j < EigValues.size() - i; j++) {
            if (fabs(EigValues[j - 1]) < fabs(EigValues[j])) {
                temp = EigValues[j - 1];
                EigValues[j - 1] = EigValues[j];
                EigValues[j] = temp;
            }
        }
 }

        Sigma1 = EigValues[0];
        Sigma2 = EigValues[1];
        Sigma3 = EigValues[2];
    
    
 //   std::cout << EigValues << std::endl;
    
   // std::cout << (EigValues[0]+EigValues[1]+EigValues[2])/3 << std::endl;
    
   
    
    REAL i1, i2, i3, j1, j2, j3;
    
    i1 = T(0,0) + T(1,1) + T(2,2);
    TPZFMatrix<REAL> T2(3,3,0.);
    T.Multiply(T,T2); // Multiplica o Tensor?
    
    i2 = 0.5*( i1*i1 - (T2(0,0) + T2(1,1) + T2(2,2)) );
    
    i3 = T(0,0)*T(1,1)*T(2,2) + T(0,1)*T(1,2)*T(2,0) + T(0,2)*T(1,0)*T(2,1) - T(0,2)*T(1,1)*T(2,0) - T(1,2)*T(2,1)*T(0,0) - T(2,2)*T(0,1)*T(1,0);
    
    j1 = 0.;
    j2 = 1./3.*(i1*i1 - 3.*i2);
    j3 = 1./27.*(2.*i1*i1*i1 - 9.*i1*i2 + 27.*i3);
    
    
    
       // std::cout << sqrt(j2) << std::endl;
    
       //  std::cout << j2 << std::endl;
    
      //  std::cout << sqrt(3*j2) << std::endl;
    
    
    
    //********* J2  *********//
    
   
    //	Sig1 (<ou>0)
    if(var == 28)
    {
        Solout[0] = j2;
        return;
    }

    

    //********* Sandler-DiMaggio  *********// If F1 <=0 --> ok
    
    REAL F1 = 0., Ff = 0., OmBeta = 0.;
    
//    std::cout << fA << std::endl;
//    std::cout << fC << std::endl;
//    std::cout << fB << std::endl;
//    
    
//    Ff = fA - pow(fC,fB*i1);
    
    Ff = fA - (fC * exp(fB*i1));
    
    OmBeta = 1.0; // Modelo clássico de Sandler-DiMaggio (tese Diogo)
    
    F1 = sqrt(j2) - (Ff/OmBeta);
    
   // std::cout << F1 << std::endl;
    
    
    /******** Sandler-DiMaggio ********/
    
    //	F1 <= 0   -> Failure
    if(var == 29)
    {
        Solout[0] = F1;
        return;
    }

    
    
    //********* Tensões Princiais *********//
    
    
    //	Sig1 (<ou>0)
    if(var == 30)
    {
        Solout[0] = i1;
        return;
    }
    
    
    if(var == 31)
    {
        Solout[0] = Sigma1;
        return;
    }
    

    if(var == 32)
    {
        Solout[0] = Sigma2;
        return;
    }
    
    
    if(var == 33)
    {
        Solout[0] = Sigma3;
        return;
    }
    
    
    
    
    //********* Checking Von Misses *********//
    
    
    //	Sig1 (<ou>0)
    if(var == 34)
    {
        if (fabs(Sigma1) < fabs(sqrt(3*j2))){
            Solout[0] = -1;
        }
        else{
            Solout[0] = 1;
        }
        
        return;
    }
    
    
    if(var == 35)
    {
        if (fabs(Sigma2) < fabs(sqrt(3*j2))){
            Solout[0] = -1;
        }
        else{
            Solout[0] = 1;
        }
        
        return;
    }
    
    
    if(var == 36)
    {
        if (fabs(Sigma3) < fabs(sqrt(3*j2))){
            Solout[0] = -1;
        }
        else{
            Solout[0] = 1;
        }
        
        return;
        
    }
    
//    //********* Checking Von Misses *********//
//    
//    
//    //	Sig1 (<ou>0)
//    if(var == 34)
//    {
//        if (Sigma1 < sqrt(3*j2)){
//            Solout[0] = -1;
//        }
//        else{
//            Solout[0] = 1;
//        }
//        
//        return;
//    }
//    
//    
//    if(var == 35)
//    {
//        if (Sigma2 < sqrt(3*j2)){
//            Solout[0] = -1;
//        }
//        else{
//            Solout[0] = 1;
//        }
//        
//        return;
//    }
//    
//    
//    if(var == 36)
//    {
//        if (Sigma3 < sqrt(3*j2)){
//            Solout[0] = -1;
//        }
//        else{
//            Solout[0] = 1;
//        }
//        
//        return;
//    }
    
    
    
    //********* Mogi-Coulomb *********////	F (<ou=0) //Referencia: Aslannezhad, Manshad, Jalalifar
    
    
    REAL SigM2 = 0., TauOct = 0., a = 0., b = 0., FMC = 0.;
    
    
    SigM2 = (EigValues[0]+EigValues[2])/2;
    
    TauOct = 1/3*(sqrt(pow((EigValues[0]-EigValues[1]), 2)+pow((EigValues[1]-EigValues[2]), 2)+pow((EigValues[2]-EigValues[0]), 2)));
                  
    a = ((2*sqrt(2))/3)*fc*cos(ffriction);
    
    b = ((2*sqrt(2))/3)*sin(ffriction);
    
    FMC = (a + b * SigM2) - TauOct;
    
    
    if(var == 37)
    {
         Solout[0] = FMC;
             return;
    }
    
    
    
    
    
    ///////////////************************************************     ELASTOPLASTICITY PROJECTED!!!!!!!!!!     ****************************************************///////////////
    
    
    //Stress Tensor
    TPZFNMatrix<9> TP(3,3,0.);
    
    
    //T = STRESS TENSOR
    TP.PutVal(0,0, (SigXP));
    TP.PutVal(0,1, (TauP));
    TP.PutVal(0,2, SigmaXZP);  /// ZERO MESMO????
    TP.PutVal(1,0, (TauP));
    TP.PutVal(1,1, (SigYP));
    TP.PutVal(1,2, SigmaYZP); /// ZERO MESMO????
    TP.PutVal(2,0, SigmaXZP); /// ZERO MESMO????
    TP.PutVal(2,1, SigmaYZP); /// ZERO MESMO????
    TP.PutVal(2,2, (SigZP));
    
    
    //    std::cout << TP << std::endl;
    
    
//    long NumIt = 1000;
//    REAL tol = 1.E-5;
    TPZVec<REAL> EigValuesP(3,0.);
    TPZFNMatrix<9,REAL> EigVectorsP(3,3,0.);
    bool EigenWorksP;
    EigenWorksP = TP.SolveEigensystemJacobi(NumIt, tol, EigValuesP, EigVectorsP);
    
    // std::cout << EigValues << std::endl;
    
    REAL Sigma1P = 0., Sigma2P = 0., Sigma3P = 0.;
    
    //**********// Criar metodo para garantir que Sig1 > Sig2 > Sig3 (em modulo)
    REAL tempP;
    for (int i = 0; i < EigValuesP.size() - 1; i++) {
        for (int j = 1; j < EigValuesP.size() - i; j++) {
            if (fabs(EigValuesP[j - 1]) < fabs(EigValuesP[j])) {
                tempP = EigValuesP[j - 1];
                EigValuesP[j - 1] = EigValuesP[j];
                EigValuesP[j] = tempP;
            }
        }
    }
    
    Sigma1P = EigValuesP[0];
    Sigma2P = EigValuesP[1];
    Sigma3P = EigValuesP[2];
    
    
    //   std::cout << EigValues << std::endl;
    
    // std::cout << (EigValues[0]+EigValues[1]+EigValues[2])/3 << std::endl;
    
    
    
    REAL i1P, i2P, i3P, j1P, j2P, j3P;
    
    i1P = TP(0,0) + TP(1,1) + TP(2,2);
    TPZFMatrix<REAL> T2P(3,3,0.);
    TP.Multiply(TP,T2P); // Multiplica o Tensor?
    
    i2P = 0.5*( i1P*i1P - (T2P(0,0) + T2P(1,1) + T2P(2,2)) );
    
    i3P = TP(0,0)*TP(1,1)*TP(2,2) + TP(0,1)*TP(1,2)*TP(2,0) + TP(0,2)*TP(1,0)*TP(2,1) - TP(0,2)*TP(1,1)*TP(2,0) - TP(1,2)*TP(2,1)*TP(0,0) - TP(2,2)*TP(0,1)*TP(1,0);
    
    j1P = 0.;
    j2P = 1./3.*(i1P*i1P - 3.*i2P);
    j3P = 1./27.*(2.*i1P*i1P*i1P - 9.*i1P*i2P + 27.*i3P);
    
    
    
    // std::cout << sqrt(j2) << std::endl;
    
    //  std::cout << j2 << std::endl;
    
    //  std::cout << sqrt(3*j2) << std::endl;
    
    
    
    //********* J2  *********//
    
    
    //	Sig1 (<ou>0)
    if(var == 38)
    {
        Solout[0] = j2P;
        return;
    }
    
    
    
    //********* Sandler-DiMaggio  *********// If F1 <=0 --> ok
    
    REAL F1P = 0., FfP = 0., OmBetaP = 0.;
    
    //    std::cout << fA << std::endl;
    //    std::cout << fC << std::endl;
    //    std::cout << fB << std::endl;
    //
    
    //    Ff = fA - pow(fC,fB*i1);
    
    FfP = fA - (fC * exp(fB*i1P));
    
    OmBetaP = 1.0; // Modelo clássico de Sandler-DiMaggio (tese Diogo)
    
    F1P = sqrt(j2P) - (FfP/OmBetaP);
    
    // std::cout << F1 << std::endl;
    
    
    /******** Sandler-DiMaggio ********/
    
    //	F1 <= 0   -> Failure
    if(var == 39)
    {
        Solout[0] = F1P;
        return;
    }
    
    
    //	Sig1 (<ou>0)
    if(var == 40)
    {
        Solout[0] = i1P;
        return;
    }
    
   
    
    
    
    
}




/*********************************************************** NOT USING ********************************************************************/

/*
 void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc)
 {
 TPZFMatrix<REAL> &phi = data.phi;
 const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
 int dim = Dimension();
 int nstate = NStateVariables();
 
 const int phr = phi.Rows();
 int in,jn,idf,jdf;
 REAL v2[2];
 v2[0] = bc.Val2()(0,0);
 v2[1] = bc.Val2()(1,0);
 
 if (this->fForcingFunction) {
 
 }
 
 TPZFMatrix<REAL> &v1 = bc.Val1();
 switch (bc.Type()){
 case 0: // Dirichlet condition
 for(in = 0 ; in < phr; in++){
 ef(nstate*in+0,0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * phi(in,0) * weight;
 ef(nstate*in+1,0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * phi(in,0) * weight;
 
 for (jn = 0 ; jn < phr; jn++) {
 ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
 ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
 
 }//jn
 }//in
 break;
 
 case 1: // Neumann condition
 for(in = 0 ; in < phi.Rows(); in++) {
 ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
 ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
 }
 break;
 
 case 2: // Mixed condition
 {
 TPZFNMatrix<2,STATE> res(2,1,0.);
 for(int i=0; i<2; i++) for(int j=0; j<2; j++)
 {
 res(i,0) += bc.Val1()(i,j)*data.sol[0][j];
 }
 
 for(in = 0 ; in < phi.Rows(); in++) {
 ef(nstate*in+0,0) += (v2[0]-res(0,0)) * phi(in,0) * weight;
 ef(nstate*in+1,0) += (v2[1]-res(1,0)) * phi(in,0) * weight;
 for(jn=0; jn<phi.Rows(); jn++)
 {
 for(idf=0; idf<2; idf++) for(jdf=0; jdf<2; jdf++)
 {
 ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
 //BUG FALTA COLOCAR VAL2
 //DebugStop();
 }
 }
 }//in
 }
 break;
 
 case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
 for(in = 0 ; in < phr; in++) {
 ef(nstate*in+0,0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in,0) * weight;
 ef(nstate*in+1,0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in,0) * weight;
 for (jn = 0 ; jn < phr; jn++) {
 ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[0];
 ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[1];
 }//jn
 }//in
 break;
 
 case 4: // stressField Neumann condition
 for(in = 0; in < dim; in ++)
 v2[in] = ( v1(in,0) * data.normal[0] +
 v1(in,1) * data.normal[1]);
 // The normal vector points towards the neighbour. The negative sign is there to
 // reflect the outward normal vector.
 for(in = 0 ; in < phi.Rows(); in++) {
 ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
 ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
 //	cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << endl;
 //	cout << "val2:  " << v2[0]  << endl;
 }
 break;
 
 case 6://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
 {
 TPZFNMatrix<2,STATE> res(2,1,0.);
 for(int i=0; i<2; i++) for(int j=0; j<2; j++)
 {
 res(i,0) += bc.Val1()(i,j)*data.sol[0][j];
 }
 for(in = 0 ; in < phi.Rows(); in++)
 {
 ef(nstate*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phi(in,0) * weight ;
 ef(nstate*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phi(in,0) * weight ;
 for(jn=0; jn<phi.Rows(); jn++)
 {
 for(idf=0; idf<2; idf++) for(jdf=0; jdf<2; jdf++)
 {
 ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
 //BUG FALTA COLOCAR VAL2
 //                        DebugStop();
 }
 }
 
 }
 
 }
 break;
 case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
 {
 TPZFNMatrix<2,STATE> res(2,1,0.);
 for(int i=0; i<2; i++) for(int j=0; j<2; j++)
 {
 res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
 }
 for(in = 0 ; in < phi.Rows(); in++)
 {
 ef(nstate*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phi(in,0) * weight ;
 ef(nstate*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phi(in,0) * weight ;
 for(jn=0; jn<phi.Rows(); jn++)
 {
 for(idf=0; idf<2; idf++) for(jdf=0; jdf<2; jdf++)
 {
 ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*data.normal[idf]*data.normal[jdf]*phi(in,0)*phi(jn,0)*weight;
 //BUG FALTA COLOCAR VAL2
 //                        DebugStop();
 }
 }
 
 }
 }
 break;
 
 default:
 PZError << "TPZMatElastoPlastic2D::ContributeBC error - Wrong boundary condition type" << std::endl;
 }
 //cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
 //cout << "val2:  " << v2[0] << endl;
 }
 */


//void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
//{
//    
//    DebugStop();
//    
//    TPZFMatrix<REAL>  &phiu = data.phi;
//    TPZManVector<REAL,3> sol_u = data.sol[0];
//    TPZFMatrix<REAL> dsol_u = data.dsol[0];
//    
//    REAL ux = sol_u[0];
//    REAL uy = sol_u[1];
//    
//   
//    
//    int phru = phiu.Rows();
//    short in;
//    STATE v2[3]; TPZFMatrix<STATE> &v1 = bc.Val1();
//    v2[0] = bc.Val2()(0,0);	//	Ux displacement or Tnx
//    v2[1] = bc.Val2()(1,0);	//	Uy displacement or Tny
//    
//    //	Here each digit represent an individual boundary condition corresponding to each state variable.
//    //	0 means Dirichlet condition on x-y
//    //	1 means Neumann condition
//    //	7 means Dirichlet condition on x
//    //	8 means Dirichlet condition on y
//    
//    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
//    switch (bc.Type())
//    {
//        case 0 :
//        {
//            //	Dirichlet condition for each state variable
//            //	Elasticity Equation
//            for(in = 0 ; in < phru; in++)
//            {
//                //	Contribution for load Vector
//                ef(2*in,0)      += BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;	// X displacement Value
//                ef(2*in+1,0)	+= BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;	// y displacement Value
//                
//            }
//            
//            break;
//        }
//        case 1 :
//        {
//            //	Neumann condition for each state variable
//            //	Elasticity Equation
//            for(in = 0 ; in <phru; in++)
//            {
//                //	Normal Tension Components on neumann boundary
//                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;		//	Tnx
//                ef(2*in+1,0)	+= -1.0*v2[1]*phiu(in,0)*weight;		//	Tny
//            }
//            break;
//        }
//        case 2 :
//        {
//            //	Mixed condition for each state variable no used here
//            //	Elasticity Equation
//            TPZFNMatrix<2,STATE> res(2,1,0.);
//            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
//            {
//                res(i,0) += bc.Val1()(i,j)*data.sol[0][j];
//            }
//            
//            for(in = 0 ; in < phru; in++)
//            {
//                ef(2*in+0,0) += weight * (v2[0]-res(0,0)) * phiu(in,0);
//                ef(2*in+1,0) += weight * (v2[1]-res(1,0)) * phiu(in,0);
//                
//            }
//            
//            break;
//        }
//        case 3 :
//        {
//            //	Null Dirichlet condition for each state variable
//            //	Elasticity Equation
//            for(in = 0 ; in < phru; in++)
//            {
//                //	Contribution for load Vector
//                ef(2*in,0)      += BIGNUMBER*(0.0 - v2[0])*phiu(in,0)*weight;	// X displacement Value
//                ef(2*in+1,0)	+= BIGNUMBER*(0.0 - v2[1])*phiu(in,0)*weight;	// y displacement Value
//                
//            }
//            
//            break;
//        }
//        case 4 :
//        {
//            //	Stress Field as Neumann condition for each state variable
//            //	Elasticity Equation
//            
//            for(in = 0; in < this->Dimension(); in ++){ v2[in] = ( v1(in,0) * data.normal[0] + v1(in,1) * data.normal[1]);}
//            
//            for(in = 0 ; in <phru; in++)
//            {
//                //	Normal Tension Components on neumann boundary
//                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;        //	Tnx
//                ef(2*in+1,0)	+= -1.0*v2[1]*phiu(in,0)*weight;		//	Tny
//            }
//            
//            break;
//        }
//        case 5 :
//            //	Normal Pressure condition Pressure value Should be inserted in v2[0]
//            //	Elasticity Equation
//        {
//            TPZFNMatrix<2,STATE> res(2,1,0.);
//            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
//            {
//                res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
//            }
//            for(int in = 0 ; in < phru; in++)
//            {
//                ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
//                ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
//            }
//        }
//            break;
//        case 6 :
//            //	Normal Pressure condition Pressure value Should be inserted in v2[0]
//            //	Elasticity Equation
//        {
//            TPZFNMatrix<2,STATE> res(2,1,0.);
//            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
//            {
//                res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
//            }
//            for(int in = 0 ; in < phru; in++)
//            {
//                ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
//                ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
//            }
//        }
//            break;
//        case 7 :
//        {
//            //	Dirichlet condition for each state variable
//            //	Elasticity Equation
//            for(in = 0 ; in < phru; in++)
//            {
//                //	Contribution for load Vector
//                ef(2*in,0)		+= BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;	// X displacement Value
//            }
//            
//            break;
//        }
//        case 8 :
//        {
//            //	Dirichlet condition for each state variable
//            //	Elasticity Equation
//            for(in = 0 ; in < phru; in++)
//            {
//                //	Contribution for load Vector
//                ef(2*in+1,0)	+= BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;	// y displacement Value
//            }
//            
//            break;
//        }
//        default:
//        {
//            PZError << "TPZMatElasticity2D::ContributeBC error - Wrong boundary condition type" << std::endl;
//            DebugStop();
//        }
//            break;
//    }
//    
//}
//





//
//void TPZMatElasticity2D::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
//{
//
//TPZFMatrix<REAL>  &phiu = data.phi;
//TPZManVector<REAL,3> sol_u = data.sol[0];
//TPZFMatrix<REAL> dsol_u = data.dsol[0];
//
//REAL ux = sol_u[0];
//REAL uy = sol_u[1];
//
//int phru = phiu.Rows();
//short in,jn;
//TPZManVector<STATE,3> v2(3);
//TPZFMatrix<STATE> &v1 = bc.Val1();
//v2[0] = bc.Val2()(0,0);	//	Ux displacement or Tnx
//v2[1] = bc.Val2()(1,0);	//	Uy displacement or Tny
//
////	Here each digit represent an individual boundary condition corresponding to each state variable.
////	0 means Dirichlet condition on x-y
////	1 means Neumann condition
////	7 means Dirichlet condition on x
////	8 means Dirichlet condition on y
//
//const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
//switch (bc.Type())
//{
//    case 0 :
//    {
//        //	Dirichlet condition for each state variable
//        //	Elasticity Equation
//        for(in = 0 ; in < phru; in++)
//        {
//            //	Contribution for load Vector
//            ef(2*in,0)      += (BIGNUMBER*(ux - v2[0])*phiu(in,0))*weight;	// X displacement Value
//            ef(2*in+1,0)	+= (BIGNUMBER*(uy - v2[1])*phiu(in,0))*weight;	// y displacement Value
//            
////            for (jn = 0 ; jn < phru; jn++)
////            {
////                //	Contribution for Stiffness Matrix
////                ek(2*in,2*jn)       += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
////                ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
////            }
//        }
//        
//        break;
//    }
//        
//    case 1 :
//    {
//        //	Neumann condition for each state variable
//        //	Elasticity Equation
//        for(in = 0 ; in <phru; in++)
//        {
//            //	Normal Tension Components on neumann boundary
//            ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;		//	Tnx
//            ef(2*in+1,0)	+= -1.0*v2[1]*phiu(in,0)*weight;		//	Tny
//        }
//        break;
//    }
//        
//    case 2 :
//    {
//        //	Mixed condition for each state variable no used here
//        //	Elasticity Equation
//        TPZFNMatrix<2,STATE> res(2,1,0.);
//        for(int i=0; i<2; i++) for(int j=0; j<2; j++)
//        {
//            res(i,0) += bc.Val1()(i,j)*data.sol[0][j];
//        }
//        
//        for(in = 0 ; in < phru; in++)
//        {
//            ef(2*in+0,0) += weight * ((v2[0]-res(0,0)) * phiu(in,0));
//            ef(2*in+1,0) += weight * ((v2[1]-res(1,0)) * phiu(in,0));
//            
////            for (jn = 0 ; jn < phru; jn++)
////            {
////                for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf< this->Dimension(); jdf++)
////                {
////                    ek(2*in+idf,2*jn+jdf) += v1(idf,jdf)*phiu(in,0)*phiu(jn,0)*weight;
////                    //      Not Complete with val2? HERE! PHIL!!!!
////                    //      DebugStop();
////                }
////            }
//        }
//        
//        break;
//    }
//        
//    case 3 :
//    {
//        //	Null Dirichlet condition for each state variable
//        //	Elasticity Equation
//        for(in = 0 ; in < phru; in++)
//        {
//            //	Contribution for load Vector
//            ef(2*in,0)      += (BIGNUMBER*(0.0 - v2[0])*phiu(in,0))*weight;	// X displacement Value
//            ef(2*in+1,0)	+= (BIGNUMBER*(0.0 - v2[1])*phiu(in,0))*weight;	// y displacement Value
//            
////            for (jn = 0 ; jn < phru; jn++)
////            {
////                //	Contribution for Stiffness Matrix
////                ek(2*in,2*jn)       += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
////                ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
////            }
//        }
//        
//        break;
//    }
//        
//    case 4 :
//    {
//        // Stress field
//        
//        // BC as a function of the Analytic Solution
//        if (fAnalytics == 1||2) {
//            REAL theta = 0.;
//            REAL coordY = 0.;
//            REAL coordX = 0.;
//            REAL r = 0.;
//            coordX = data.x[0];
//            coordY = data.x[1];
//            theta = atan2(coordY,coordX);
//            r = sqrt((coordX*coordX)+(coordY*coordY));
//            
//            REAL Sx=0., Sy=0., Sxy=0., Sz=0., Sxz=0., Syz=0.;
//            
//            AnalyticalWellboreSolution(Sx, Sy, Sxy, Sz, Sxz, Syz, theta, r);
//            
//            v1(0,0) = Sx;
//            v1(0,1) = Sxy;
//            v1(1,0) = Sxy;
//            v1(1,1) = Sy;
//            
//            //v1.Print(" Valor de v1 ");
//            
//        }
//        
//        for(in = 0; in < this->Dimension(); in ++){
//            v2[in] = ( v1(in,0) * data.normal[0] + v1(in,1) * data.normal[1]);
//        }
//        
//        for(in = 0 ; in <phru; in++)
//        {
//            //	Normal Tension Components on neumann boundary
//            ef(2*in,0)      += 1.0*v2[0]*phiu(in,0)*weight;      //	Tnx
//            ef(2*in+1,0)	+= 1.0*v2[1]*phiu(in,0)*weight;		//	Tny
//        }
//        
//        break;
//    }
//        
//    case 5 :
//    {
//        DebugStop();
//    }
//        break;
//        
//        
//    case 6 :
//        //	Normal Pressure condition Pressure value Should be inserted in v2[0]
//        //	Elasticity Equation
//    {
//        
//        TPZManVector<REAL> n = data.normal;
//        //                TPZManVector<REAL> n_ab = data.normal;
//        //
//        //                REAL lxx = 0., lxy = 0., lxz = 0., lyx =0., lyy = 0., lyz = 0., lzx = 0., lzy = 0., lzz = 0.;
//        //
//        //                // x-diretion
//        //                lxx = cos(falpha)*cos(fbeta);
//        //                lxy = sin(falpha)*cos(fbeta);
//        //                lxz = -sin(fbeta);
//        //                // y-direction
//        //                lyx = -sin(falpha);
//        //                lyy = cos(falpha);
//        //                lyz = 0;
//        //                // z-direction
//        //                lzx = cos(falpha)*sin(fbeta);
//        //                lzy = sin(falpha)*sin(fbeta);
//        //                lzz = cos(fbeta);
//        //
//        //                n_ab[0] = lxx*n[0] + lxy*n[1] + lxz*n[2];
//        //                n_ab[1] = lyx*n[0] + lyy*n[1] + lyz*n[2];
//        //                n_ab[2] = lzx*n[0] + lzy*n[1] + lzz*n[2];
//        
//        TPZFNMatrix<2,STATE> Tn(2,1,0.);
//        for(int i=0; i<2; i++)
//        {
//            for(int j=0; j<2; j++)
//            {
//                Tn(i,0) += bc.Val1()(i,j)*n[j];
//            }
//        }
//        
//        for(int in = 0 ; in < phru; in++)
//        {
//            ef(2*in+0,0) += weight * Tn(0,0)* phiu(in,0);
//            ef(2*in+1,0) += weight * Tn(1,0) * phiu(in,0);
//        }
//    }
//        break;
//        
//        
//    case 7 :
//    {
//        //	Dirichlet condition for each state variable
//        //	Elasticity Equation
//        for(in = 0 ; in < phru; in++)
//        {
//            //	Contribution for load Vector
//            ef(2*in,0)		+= (BIGNUMBER*(ux - v2[0])*phiu(in,0))*weight;	// X displacement Value
//            
////            for (jn = 0 ; jn < phru; jn++)
////            {
////                //	Contribution for Stiffness Matrix
////                ek(2*in,2*jn)		+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// X displacement
////            }
//        }
//        
//        break;
//    }
//        
//    case 8 :
//    {
//        //	Dirichlet condition for uy
//        //	Elasticity Equation
//        for(in = 0 ; in < phru; in++)
//        {
//            //	Contribution for load Vector
//            ef(2*in+1,0)	+= (BIGNUMBER*(uy - v2[1])*phiu(in,0))*weight;	// y displacement Value
//            
////            for (jn = 0 ; jn < phru; jn++)
////            {
////                //	Contribution for Stiffness Matrix
////                ek(2*in+1,2*jn+1)	+= BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;	// Y displacement
////            }
//        }
//        
//        break;
//    }
//        
//    default:
//    {
//        PZError << "TPZMatElasticity2D::ContributeBC error - Wrong boundary condition type" << std::endl;
//        DebugStop();
//        }
//        break;
//        }
//}
//
//

