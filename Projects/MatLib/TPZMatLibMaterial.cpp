//
//  TPZMatLibMaterial.cpp
//  PZ
//
//  Created by Philippe Devloo on 7/11/16.
//
//

#include "TPZMatLibMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"

zorglib::ParameterSet TPZMatLibMaterial::fExternalState;



/** To be implemented only in the proper materials. */
int TPZMatLibMaterial::PushMemItem(int sourceIndex)
{
    int newindex = TPZMatWithMem<TPZMatLibMemory>::PushMemItem(sourceIndex);
    if (sourceIndex < 0) {
        fMatModel.initState(MemItem(newindex).fPrevState);
        fMatModel.initState(MemItem(newindex).fNextState);
        MemItem(newindex).fK.resize(fMatModel.model().nExtVar());
        MemItem(newindex).fK = 0.;
    }
    return newindex;
}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */

void TPZMatLibMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    long memindex = data.intGlobPtIndex;
    TPZMatLibMemory &mem = MemItem(memindex);
    
#ifdef PZDEBUG
    if(mem.fPrevState.grad.size() != 6)
    {
        DebugStop();
    }
#endif
    
    int nstate = fDimension;
    
    int matmodelnvar = fMatModel.model().nExtVar();
    int matmodeldim = 3;
    if (matmodelnvar == 6) {
        matmodeldim = 3;
    }
    else if(matmodelnvar == 4)
    {
        matmodeldim = 2;
    }
    else if(matmodelnvar == 3)
    {
        matmodeldim = 1;
    }
    
    TPZFNMatrix<9,STATE> dsolxy;
    TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0],dsolxy,data.axes);
    
    const int _XX_ = 0;
    const int _XY_ = 1;
    int _YY_ = 2;
    const int _XZ_ = 3;
    const int _YZ_ = 4;
    int _ZZ_ = 5;
    if (matmodeldim == 2) {
        _ZZ_ = 3;
    }
    if (matmodeldim == 1) {
        _YY_ = 1;
        _ZZ_ = 2;
    }

    mem.fNextState.grad = 0.;
    mem.fNextState.grad(_XX_) = dsolxy(0,0);
    if (fDimension > 1)
    {
        mem.fNextState.grad(_XY_) = dsolxy(1,0) + dsolxy(0,1);
        mem.fNextState.grad(_YY_) = dsolxy(1,1);
    }
    if (fDimension > 2)
    {
        mem.fNextState.grad(_XZ_) = dsolxy(0,2)+dsolxy(2,0);
        mem.fNextState.grad(_YZ_) = dsolxy(1,2)+dsolxy(2,1);
        mem.fNextState.grad(_ZZ_) = dsolxy(2,2);
    }
    bool shouldcomputetangent = true;
    REAL delt = 1.0;

    fMatModel.updateState(fExternalState,mem.fPrevState,mem.fNextState,delt,mem.fK,shouldcomputetangent);

//    std::cout << mem.fNextState.grad << std::endl;
//    std::cout << mem.fNextState.flux << std::endl;
    
    TPZFNMatrix<200,REAL> dphiXYZ;
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix,dphiXYZ,data.axes);
    
    
    int phr = data.phi.Rows();
    
    if (fDimension == 3)
    {
        for(int in = 0; in < phr; in++) { //in: test function index
            
            // fForce represents the gravity acceleration
            //First equation: fb and fk
            STATE val  = 0.; // fb
            val -= mem.fNextState.flux(_XX_) * dphiXYZ(0,in); // |
            val -= mem.fNextState.flux(_XY_) * dphiXYZ(1,in); // fk
            val -= mem.fNextState.flux(_XZ_) * dphiXYZ(2,in); // |
            ef(in*nstate+0,0) += weight * val;
            
            //Second equation: fb and fk
            val  = 0.; // fb
            val -= mem.fNextState.flux(_XY_) * dphiXYZ(0,in); // |
            val -= mem.fNextState.flux(_YY_) * dphiXYZ(1,in); // fk
            val -= mem.fNextState.flux(_YZ_) * dphiXYZ(2,in); // |
            ef(in*nstate+1,0) += weight * val;
            
            //third equation: fb and fk
            val  = 0.; // fb
            val -= mem.fNextState.flux(_XZ_) * dphiXYZ(0,in); // |
            val -= mem.fNextState.flux(_YZ_) * dphiXYZ(1,in); // fk
            val -= mem.fNextState.flux(_ZZ_) * dphiXYZ(2,in); // |
            ef(in*nstate+2,0) += weight * val;
            
            for( int jn = 0; jn < phr; jn++ ) {
                //jn: trial function index
                //this matrix will store
                //{{dvdx*dudx, dvdx*dudy, dvdx*dudz},
                //{dvdy*dudx, dvdy*dudy, dvdy*dudz},
                //{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
                //Compute Deriv matrix
                TPZFNMatrix<9,REAL> Deriv(3,3);
                for(int ud = 0; ud < 3; ud++){
                    for(int vd = 0; vd < 3; vd++){
                        Deriv(vd,ud) = dphiXYZ(vd,in)*dphiXYZ(ud,jn);
                    }//ud
                }//vd
                
                
                //#define _XX_ 0
                //#define _XY_ 1
                //#define _XZ_ 2
                //#define _YY_ 3
                //#define _YZ_ 4
                //#define _ZZ_ 5
                //First equation Dot[Sigma1, gradV1]
                STATE val2  = mem.fK(_XX_,_XX_) * Deriv(0,0);//dvdx*dudx
                val2 +=      mem.fK(_XX_,_XY_) * Deriv(0,1);//dvdx*dudy
                val2 +=	   mem.fK(_XX_,_XZ_) * Deriv(0,2);//dvdx*dudz
                val2 +=    mem.fK(_XY_,_XX_) * Deriv(1,0);//dvdy*dudx
                val2 +=      mem.fK(_XY_,_XY_) * Deriv(1,1);//dvdy*dudy
                val2 +=      mem.fK(_XY_,_XZ_) * Deriv(1,2);//dvdy*dudz
                val2 +=      mem.fK(_XZ_,_XX_) * Deriv(2,0);//dvdz*dudx
                val2 +=      mem.fK(_XZ_,_XY_) * Deriv(2,1);//dvdz*dudy
                val2 +=      mem.fK(_XZ_,_XZ_) * Deriv(2,2);//dvdz*dudz
                ek(in*nstate+0,jn*nstate+0) += weight * val2;
                
                STATE val3  =      mem.fK(_XX_,_XY_) * Deriv(0,0);
                val3 +=      mem.fK(_XX_,_YY_) * Deriv(0,1);
                val3 +=      mem.fK(_XX_,_YZ_) * Deriv(0,2);
                val3 +=      mem.fK(_XY_,_XY_) * Deriv(1,0);
                val3 +=      mem.fK(_XY_,_YY_) * Deriv(1,1);
                val3 +=      mem.fK(_XY_,_YZ_) * Deriv(1,2);
                val3 +=      mem.fK(_XZ_,_XY_) * Deriv(2,0);
                val3 +=      mem.fK(_XZ_,_YY_) * Deriv(2,1);
                val3 +=      mem.fK(_XZ_,_YZ_) * Deriv(2,2);
                ek(in*nstate+0,jn*nstate+1) += weight * val3;
                
                STATE val4  =      mem.fK(_XX_,_XZ_) * Deriv(0,0);
                val4 +=      mem.fK(_XX_,_YZ_) * Deriv(0,1);
                val4 +=      mem.fK(_XX_,_ZZ_) * Deriv(0,2);//
                val4 +=      mem.fK(_XY_,_XZ_) * Deriv(1,0);
                val4 +=      mem.fK(_XY_,_YZ_) * Deriv(1,1);
                val4 +=      mem.fK(_XY_,_ZZ_) * Deriv(1,2);//
                val4 +=      mem.fK(_XZ_,_XZ_) * Deriv(2,0);
                val4 +=      mem.fK(_XZ_,_YZ_) * Deriv(2,1);
                val4 +=      mem.fK(_XZ_,_ZZ_) * Deriv(2,2);
                ek(in*nstate+0,jn*nstate+2) += weight * val4;
                
                //Second equation Dot[Sigma2, gradV2]
                STATE val5  =  mem.fK(_XY_,_XX_) * Deriv(0,0);
                val5 +=      mem.fK(_XY_,_XY_) * Deriv(0,1);
                val5 +=      mem.fK(_XY_,_XZ_) * Deriv(0,2);
                val5 +=      mem.fK(_YY_,_XX_) * Deriv(1,0);
                val5 +=      mem.fK(_YY_,_XY_) * Deriv(1,1);
                val5 +=      mem.fK(_YY_,_XZ_) * Deriv(1,2);
                val5 +=      mem.fK(_YZ_,_XX_) * Deriv(2,0);
                val5 +=      mem.fK(_YZ_,_XY_) * Deriv(2,1);
                val5 +=      mem.fK(_YZ_,_XZ_) * Deriv(2,2);
                ek(in*nstate+1,jn*nstate+0) += weight * val5;
                
                STATE val6  =      mem.fK(_XY_,_XY_) * Deriv(0,0);
                val6 +=      mem.fK(_XY_,_YY_) * Deriv(0,1);
                val6 +=      mem.fK(_XY_,_YZ_) * Deriv(0,2);
                val6 +=      mem.fK(_YY_,_XY_) * Deriv(1,0);
                val6 +=      mem.fK(_YY_,_YY_) * Deriv(1,1);
                val6 +=      mem.fK(_YY_,_YZ_) * Deriv(1,2);
                val6 +=      mem.fK(_YZ_,_XY_) * Deriv(2,0);
                val6 +=      mem.fK(_YZ_,_YY_) * Deriv(2,1);
                val6 +=      mem.fK(_YZ_,_YZ_) * Deriv(2,2);
                ek(in*nstate+1,jn*nstate+1) += weight * val6;
                
                STATE val7  =      mem.fK(_XY_,_XZ_) * Deriv(0,0);
                val7 +=      mem.fK(_XY_,_YZ_) * Deriv(0,1);
                val7 +=      mem.fK(_XY_,_ZZ_) * Deriv(0,2);//
                val7 +=      mem.fK(_YY_,_XZ_) * Deriv(1,0);
                val7 +=      mem.fK(_YY_,_YZ_) * Deriv(1,1);
                val7 +=      mem.fK(_YY_,_ZZ_) * Deriv(1,2);//
                val7 +=      mem.fK(_YZ_,_XZ_) * Deriv(2,0);
                val7 +=      mem.fK(_YZ_,_YZ_) * Deriv(2,1);
                val7 +=      mem.fK(_YZ_,_ZZ_) * Deriv(2,2);
                ek(in*nstate+1,jn*nstate+2) += weight * val7;
                
                //Third equation Dot[Sigma3, gradV3]
                STATE val8  =  mem.fK(_XZ_,_XX_) * Deriv(0,0);
                val8 +=      mem.fK(_XZ_,_XY_) * Deriv(0,1);
                val8 +=      mem.fK(_XZ_,_XZ_) * Deriv(0,2);
                val8 +=      mem.fK(_YZ_,_XX_) * Deriv(1,0);
                val8 +=      mem.fK(_YZ_,_XY_) * Deriv(1,1);
                val8 +=      mem.fK(_YZ_,_XZ_) * Deriv(1,2);
                val8 +=      mem.fK(_ZZ_,_XX_) * Deriv(2,0);//
                val8 +=      mem.fK(_ZZ_,_XY_) * Deriv(2,1);//
                val8 +=      mem.fK(_ZZ_,_XZ_) * Deriv(2,2);
                ek(in*nstate+2,jn*nstate+0) += weight * val8;
                
                STATE val9  =      mem.fK(_XZ_,_XY_) * Deriv(0,0);
                val9 +=      mem.fK(_XZ_,_YY_) * Deriv(0,1);
                val9 +=      mem.fK(_XZ_,_YZ_) * Deriv(0,2);
                val9 +=      mem.fK(_YZ_,_XY_) * Deriv(1,0);
                val9 +=      mem.fK(_YZ_,_YY_) * Deriv(1,1);
                val9 +=      mem.fK(_YZ_,_YZ_) * Deriv(1,2);
                val9 +=      mem.fK(_ZZ_,_XY_) * Deriv(2,0);//
                val9 +=      mem.fK(_ZZ_,_YY_) * Deriv(2,1);//
                val9 +=      mem.fK(_ZZ_,_YZ_) * Deriv(2,2);
                ek(in*nstate+2,jn*nstate+1) += weight * val9;
                
                STATE val10  =      mem.fK(_XZ_,_XZ_) * Deriv(0,0);
                val10 +=      mem.fK(_XZ_,_YZ_) * Deriv(0,1);
                val10 +=      mem.fK(_XZ_,_ZZ_) * Deriv(0,2);
                val10 +=      mem.fK(_YZ_,_XZ_) * Deriv(1,0);
                val10 +=      mem.fK(_YZ_,_YZ_) * Deriv(1,1);
                val10 +=      mem.fK(_YZ_,_ZZ_) * Deriv(1,2);
                val10 +=      mem.fK(_ZZ_,_XZ_) * Deriv(2,0);
                val10 +=      mem.fK(_ZZ_,_YZ_) * Deriv(2,1);
                val10 +=      mem.fK(_ZZ_,_ZZ_) * Deriv(2,2);//
                ek(in*nstate+2,jn*nstate+2) += weight * val10;
                
            }//jn
        }//in
    }
    else if(fDimension == 2)
    {
        for(int in = 0; in < phr; in++) { //in: test function index
            
            // fForce represents the gravity acceleration
            //First equation: fb and fk
            STATE val  = 0.; // fb
            val -= mem.fNextState.flux(_XX_) * dphiXYZ(0,in); // |
            val -= mem.fNextState.flux(_XY_) * dphiXYZ(1,in); // fk
            ef(in*nstate+0,0) += weight * val;
            
            //Second equation: fb and fk
            val  = 0.; // fb
            val -= mem.fNextState.flux(_XY_) * dphiXYZ(0,in); // |
            val -= mem.fNextState.flux(_YY_) * dphiXYZ(1,in); // fk
            ef(in*nstate+1,0) += weight * val;
            
            for( int jn = 0; jn < phr; jn++ ) {
                //jn: trial function index
                //this matrix will store
                //{{dvdx*dudx, dvdx*dudy, dvdx*dudz},
                //{dvdy*dudx, dvdy*dudy, dvdy*dudz},
                //{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
                //Compute Deriv matrix
                TPZFNMatrix<9,REAL> Deriv(fDimension,fDimension);
                for(int ud = 0; ud < fDimension; ud++){
                    for(int vd = 0; vd < fDimension; vd++){
                        Deriv(vd,ud) = dphiXYZ(vd,in)*dphiXYZ(ud,jn);
                    }//ud
                }//vd
                
                
                //#define _XX_ 0
                //#define _XY_ 1
                //#define _XZ_ 2
                //#define _YY_ 3
                //#define _YZ_ 4
                //#define _ZZ_ 5
                //First equation Dot[Sigma1, gradV1]
                STATE val2  =  mem.fK(_XX_,_XX_) * Deriv(0,0);//dvdx*dudx
                val2 +=  mem.fK(_XX_,_XY_) * Deriv(0,1);//dvdx*dudy
                val2 +=  mem.fK(_XY_,_XX_) * Deriv(1,0);//dvdy*dudx
                val2 +=  mem.fK(_XY_,_XY_) * Deriv(1,1);//dvdy*dudy
                ek(in*nstate+0,jn*nstate+0) += weight * val2;
                
                STATE val3  =  mem.fK(_XX_,_XY_) * Deriv(0,0);
                val3 +=  mem.fK(_XX_,_YY_) * Deriv(0,1);
                val3 +=  mem.fK(_XY_,_XY_) * Deriv(1,0);
                val3 +=  mem.fK(_XY_,_YY_) * Deriv(1,1);
                ek(in*nstate+0,jn*nstate+1) += weight * val3;
                
                
                //Second equation Dot[Sigma2, gradV2]
                STATE val5  =  mem.fK(_XY_,_XX_) * Deriv(0,0);
                val5 +=  mem.fK(_XY_,_XY_) * Deriv(0,1);
                val5 +=  mem.fK(_YY_,_XX_) * Deriv(1,0);
                val5 +=  mem.fK(_YY_,_XY_) * Deriv(1,1);
                ek(in*nstate+1,jn*nstate+0) += weight * val5;
                
                STATE val6  =  mem.fK(_XY_,_XY_) * Deriv(0,0);
                val6 +=  mem.fK(_XY_,_YY_) * Deriv(0,1);
                val6 +=  mem.fK(_YY_,_XY_) * Deriv(1,0);
                val6 +=  mem.fK(_YY_,_YY_) * Deriv(1,1);
                ek(in*nstate+1,jn*nstate+1) += weight * val6;
                
            }//jn
        }//in
        
    }
    else if(fDimension == 1)
    {
        for(int in = 0; in < phr; in++) { //in: test function index
            
            // fForce represents the gravity acceleration
            //First equation: fb and fk
            STATE val  = 0.; // fb
            val -= mem.fNextState.flux(_XX_) * dphiXYZ(0,in); // |
            ef(in*nstate+0,0) += weight * val;
            
            
            for( int jn = 0; jn < phr; jn++ ) {
                //jn: trial function index
                //this matrix will store
                //{{dvdx*dudx, dvdx*dudy, dvdx*dudz},
                //{dvdy*dudx, dvdy*dudy, dvdy*dudz},
                //{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
                //Compute Deriv matrix
                TPZFNMatrix<9,REAL> Deriv(fDimension,fDimension);
                for(int ud = 0; ud < fDimension; ud++){
                    for(int vd = 0; vd < fDimension; vd++){
                        Deriv(vd,ud) = dphiXYZ(vd,in)*dphiXYZ(ud,jn);
                    }//ud
                }//vd
                
                
                //#define _XX_ 0
                //#define _XY_ 1
                //#define _XZ_ 2
                //#define _YY_ 3
                //#define _YZ_ 4
                //#define _ZZ_ 5
                //First equation Dot[Sigma1, gradV1]
                STATE val2  = 2. * mem.fK(_XX_,_XX_) * Deriv(0,0);//dvdx*dudx
                val2 *= 0.5;
                ek(in*nstate+0,jn*nstate+0) += weight * val2;
                
                
            }//jn
        }//in

    }

}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition material
 * @since October 07, 2011
 */
void TPZMatLibMaterial::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype==data.EVecShape){
        DebugStop();
        return;
    }
    
    TPZFMatrix<REAL> &phi = data.phi;
    
    const STATE BIGNUMBER  = 1.e12;
    
    const int phr = phi.Rows();
    int in,jn,idf,jdf;
    TPZManVector<STATE,3> v2(fDimension);
    for (int i=0; i<fDimension; i++) {
        v2[i] = bc.Val2()(i,0);
    }
    TPZFMatrix<STATE> &v1 = bc.Val1();
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
    
    TPZVec<STATE> &sol = data.sol[0];
    int nstate = fDimension;
    
    switch (bc.Type()) {
        case 0: // Dirichlet condition
            for (int i=0; i<fDimension; i++)
            {
                for(in = 0 ; in < phr; in++) {
                    ef(nstate*in+i,0) += BIGNUMBER * (v2[i]-sol[i]) * phi(in,0) * weight;
                    for (jn = 0 ; jn < phr; jn++) {
                        ek(nstate*in+i,nstate*jn+i) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                    }//jn
                }//in
            }
            break;
            
        case 1: // Neumann condition
            for(in = 0 ; in < phi.Rows(); in++) {
                for (int i=0; i<fDimension; i++)
                {
                    ef(nstate*in+i,0) += v2[i] * phi(in,0) * weight;
                }
            }//in
            break;
        case 2: // Mixed condition
        {
            TPZManVector<STATE,3> force(fDimension,0);
            for (int i=0; i<fDimension; i++) {
                for (int j=0; j<fDimension; j++) {
                    force[i] += bc.Val1()(i,j)*(v2[j]-sol[j]);
                }
            }
            for(in = 0 ; in < phi.Rows(); in++) {
                for (int i=0; i<fDimension; i++)
                {
                    ef(nstate*in+i,0) += force[i] * phi(in,0) * weight;
                }
                for(jn=0; jn<phi.Rows(); jn++)
                {
                    for(idf=0; idf<fDimension; idf++) for(jdf=0; jdf<fDimension; jdf++)
                    {
                        ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*weight*phi(in,0)*phi(jn,0);
                    }
                }
            }//in
        }
            break;
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
        {
            for(in = 0 ; in < phr; in++) {
                for (int i=0; i<fDimension; i++)
                {
                    ef(nstate*in+i) += -BIGNUMBER * weight * sol[i] * v2[i] * phi(in,0);
                }
                for (jn = 0 ; jn < phr; jn++) {
                    for (int i=0; i<fDimension; i++)
                    {
                        ek(nstate*in+i,nstate*jn+i) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[i];
                    }
                }//jn
            }//in
        }
            break;
            
        case 4: // stressField Neumann condition
            for(in = 0; in < fDimension; in ++)
            {
                v2[in] = 0;
                for (int j=0; j<fDimension; j++)
                {
                    v2[in] -=  v1(in,j) * data.normal[j];
                }
            }
            // The normal vector points towards the neighbour. The negative sign is there to
            // reflect the outward normal vector.
            for(in = 0 ; in < phi.Rows(); in++) {
                for (int i=0; i<fDimension; i++)
                {
                    ef(nstate*in+0,0) += v2[i] * phi(in,0) * weight;
                }
                //cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
                //cout << "val2:  " << v2[0]          << ' ' << v2[1]          << ' ' << v2[2]          << endl;
            }
            break;
        default:
            PZError << "TPZElastitity3D::ContributeBC error - Wrong boundary condition type" << std::endl;
    }//switch
}//method


/// Consider the iterative process converged : the nextstate becomes previous state
// this method will loop over all the integration points and update the plastic state
void TPZMatLibMaterial::LoadState()
{
    TPZAdmChunkVector<TPZMatLibMemory> &mem = GetMemory();
    long nel = mem.NElements();
    for (long el=0; el<nel; el++) {
        mem[el].fPrevState = mem[el].fNextState;
    }
}

/** @brief Returns the variable index associated with the name */
int TPZMatLibMaterial::VariableIndex(const std::string &name)
{
    for (int i=0; i<fZorgLibVars.size(); i++) {
        if (name == fZorgLibVars[i].fName) {
            return fZorgLibVars[i].fPZIndex;
        }
    }
    return TPZMaterial::VariableIndex(name);
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatLibMaterial::NSolutionVariables(int var)
{
    if (var < 10 || var >= 99) {
        return TPZMaterial::NSolutionVariables(var);
    }
    if(var-10 >= fZorgLibVars.size())
    {
        DebugStop();
    }
    return fZorgLibVars[var-10].fNumVars;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatLibMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    
    if (var < 10 || var >= 99) {
        TPZMaterial::Solution(data, var, Solout);
        return;
    }
    long intpointindex = data.intGlobPtIndex;
    if (intpointindex == -1) {
        DebugStop();
    }
    TPZMatLibMemory &mem = MemItem(intpointindex);
    
    ZorgLibVar a = fZorgLibVars[var-10];
    int numvar = a.fNumVars;
    Solout.resize(numvar);

    TransferData(a, mem.fNextState, Solout);
}

void TPZMatLibMaterial::TransferData(ZorgLibVar &a, zorglib::MaterialState &state, TPZVec<STATE> &solution)
{
    const int _XX_ = 0;
    const int _XY_ = 1;
    int _YY_ = 2;
    const int _XZ_ = 3;
    const int _YZ_ = 4;
    int _ZZ_ = 5;
    if (fZorgLibDimension == 2) {
        _ZZ_ = 3;
    }
    if (fZorgLibDimension == 1) {
        _YY_ = 1;
        _ZZ_ = 2;
    }

    int firstindex = a.fZorgIndex;
    if(a.fNumVars < 6)
    {
        if (a.fIsExternal) {
            // help
            if(a.fIsGrad)
            {
                for (int i=0; i<a.fNumVars; i++)
                {
                    solution[i] = state.grad(i+firstindex);
                }
            }
            else
            {
                for (int i=0; i<a.fNumVars; i++)
                {
                    solution[i] = state.flux(i+firstindex);
                }
                
            }
        }
        else
        {
            
            for (int i=0; i<a.fNumVars; i++) {
                solution[i] = state.internal(i+firstindex);
            }
        }

    }
    else
    {
        solution.Fill(0.);
        zorglib::MatLibArray *orig = 0;
        if (a.fIsExternal) {
            // help
            if(a.fIsGrad)
            {
                orig = &state.grad;
            }
            else
            {
                orig = &state.flux;
            }
        }
        else
        {
            orig = &state.internal;
        }
        solution[0] = (*orig)(firstindex+_XX_);
        solution[4] = (*orig)(firstindex+_YY_);
        solution[8] = (*orig)(firstindex+_ZZ_);
        if (fZorgLibDimension > 1) {
            solution[1] = (*orig)(firstindex+_XY_);
            solution[3] = solution[1];
            if (fZorgLibDimension > 2) {
                solution[2] = (*orig)(_XZ_);
                solution[6] = solution[2];
                solution[5] = (*orig)(_YZ_);
                solution[7] = solution[5];
            }
        }
        
    }
}
