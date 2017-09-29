/*
 *  TPZStokesMaterial.cpp
 *  PZ
 *
 *  Created by Thiago Dias dos Santos on 12/01/2015.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZStokesMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzmatwithmem.h"
#include "pzfmatrix.h"

//SetSpace
//#define IsHDivQ
//#define IsH1
#define IsDGM


TPZStokesMaterial::TPZStokesMaterial() : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(){
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1;
    
}

////////////////////////////////////////////////////////////////////

TPZStokesMaterial::TPZStokesMaterial(int matid, int dimension, REAL viscosity, REAL theta, STATE Sigma) : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(matid),fViscosity(viscosity),fTheta(theta),fDimension(dimension),fSigma(Sigma)
{
    // symmetric version
    //fTheta = -1;
    
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1;
    
}

////////////////////////////////////////////////////////////////////

TPZStokesMaterial::TPZStokesMaterial(const TPZStokesMaterial &mat) : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(mat), fViscosity(mat.fViscosity), fTheta(mat.fTheta),fDimension(mat.fDimension),fSigma(mat.fSigma)
{
    fk= mat.fk;
    
}

////////////////////////////////////////////////////////////////////

TPZStokesMaterial::~TPZStokesMaterial(){
    
    
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

////////////////////////////////////////////////////////////////////

int TPZStokesMaterial::VariableIndex(const std::string &name) {
    
    if (!strcmp("P", name.c_str()))  return 0;
    if (!strcmp("V", name.c_str()))  return 1;
    if (!strcmp("f", name.c_str()))         return 2;
    if (!strcmp("V_exact", name.c_str()))   return 3;
    if (!strcmp("P_exact", name.c_str()))   return 4;
    //    if (!strcmp("V_exactBC", name.c_str()))   return 5;
    
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

////////////////////////////////////////////////////////////////////

int TPZStokesMaterial::NSolutionVariables(int var) {
    
    switch(var) {
            
        case 0:
            return 1; // Pressure, Scalar
        case 1:
            return this->Dimension(); // Velocity, Vector
        case 2:
            return this->Dimension(); // f, Vector
        case 3:
            return this->Dimension(); // V_exact, Vector
        case 4:
            return this->Dimension(); // P_exact, Vector
            //        case 5:
            //            return this->Dimension(); // V_exactBC, Vector
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    
    //itapopo conferir esse metodo
    
    int vindex = this->VIndex();
    int pindex = this->PIndex();
    
    TPZManVector<REAL,3> v_h = datavec[vindex].sol[0];
    REAL p_h = datavec[pindex].sol[0][0];
    
    // TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    // TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
            
        case 0: //Pressure
        {
            Solout[0] = p_h;
        }
            break;
            
        case 1: //Velocity
        {
            Solout[0] = v_h[0]; // Vx
            Solout[1] = v_h[1]; // Vy
        }
            break;
        case 2: //f
        {
            TPZVec<double> f;
            if(this->HasForcingFunction()){
                this->ForcingFunction()->Execute(datavec[vindex].x, f);
            }
            Solout[0] = f[0]; // fx
            Solout[1] = f[1]; // fy
        }
            break;
            
        case 3: //v_exact
        {
            TPZVec<double> v;
            if(this->HasfForcingFunctionExact()){
                this->ForcingFunctionExact()->Execute(datavec[vindex].x, v);
            }
            Solout[0] = v[0]; // vx
            Solout[1] = v[1]; // vy
        }
            break;
            
        case 4: //p_exact
        {
            TPZVec<double> p;
            if(this->HasfForcingFunctionExact()){
                this->ForcingFunctionExactPressure()->Execute(datavec[pindex].x, p);
            }
            Solout[0] = p[0]; // px
            
        }
            break;
            
            //        case 5: //v_exact
            //        {
            //            TPZVec<double> vbc;
            //            if(this->HasffBCForcingFunction()){
            //                this->ForcingFunctionBC()->Execute(datavec[vindex].x, vbc);
            //            }
            //            Solout[0] = vbc[0]; // vbcx
            //            Solout[1] = vbc[1]; // vbcy
            //        }
            //            break;
            
            
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

////////////////////////////////////////////////////////////////////

// Divergence on deformed element
void TPZStokesMaterial::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)
{
    
    //itapopo conferir esse método. Foi copiado do TPZDarcyFlow3D
    
    int ublock = 0;
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFMatrix<STATE> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFMatrix<STATE> Qaxes = datavec[ublock].axes;
    TPZFMatrix<STATE> QaxesT;
    TPZFMatrix<STATE> Jacobian = datavec[ublock].jacobian;
    TPZFMatrix<STATE> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFMatrix<STATE> GradOfX;
    TPZFMatrix<STATE> GradOfXInverse;
    TPZFMatrix<STATE> VectorOnMaster;
    TPZFMatrix<STATE> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    if (HDivPiola == 1)
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            VectorOnXYZ(0,0) = datavec[ublock].fNormalVec(0,ivectorindex);
            VectorOnXYZ(1,0) = datavec[ublock].fNormalVec(1,ivectorindex);
            VectorOnXYZ(2,0) = datavec[ublock].fNormalVec(2,ivectorindex);
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            DivergenceofPhi(iq,0) =  (1.0/JacobianDet) * ( dphiuH1(0,ishapeindex)*VectorOnMaster(0,0) +
                                                          dphiuH1(1,ishapeindex)*VectorOnMaster(1,0) +
                                                          dphiuH1(2,ishapeindex)*VectorOnMaster(2,0) );
        }
    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            DivergenceofPhi(iq,0) =  datavec[ublock].fNormalVec(0,ivectorindex)*GradphiuH1(0,ishapeindex) +
            datavec[ublock].fNormalVec(1,ivectorindex)*GradphiuH1(1,ishapeindex) +
            datavec[ublock].fNormalVec(2,ivectorindex)*GradphiuH1(2,ishapeindex) ;
        }
    }
    
    return;
    
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
    
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::Read(TPZStream &buf, void *context) {
    
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<STATE> > &GradPhi){
    
    
    TPZFMatrix<STATE> &dphiV = dataV.dphix;
    
    const int dim = this->Dimension();
    
    GradPhi.clear();
    GradPhi.resize(dim);
    
    //for each shape
    for(int shape = 0; shape < dphiV.Rows(); shape++){
        
        TPZFMatrix<STATE> GPhi(dim,dim,0.);
        
        for(int i = 0; i < dim; i++){
            
            for(int j = 0; j < dim; j++){
                
                GPhi(i,j) = dphiV(j,shape);// itapopo H1 ??
                
            }//j
        }//i
        
        GradPhi[shape] = GPhi;
        
    }//shape
    
}

// Contricucao dos elementos internos

void TPZStokesMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavec[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[vindex]);
    }
    // Setting forcing function
    /*STATE force = 0.;
     if(this->fForcingFunction) {
     TPZManVector<STATE> res(1);
     fForcingFunction->Execute(datavec[pindex].x,res);
     force = res[0];
     }*/
    
    //Gravity
    STATE rhoi = 900.; //itapopo
    STATE g = 9.81; //itapopo
    STATE force = rhoi*g;
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    
    TPZVec<double> f;
    TPZFMatrix<STATE> phiVi(fDimension,1,0.0);
    
    
    for(int i = 0; i < nshapeV; i++ )
    {
        int iphi = datavec[vindex].fVecShapeIndex[i].second;
        int ivec = datavec[vindex].fVecShapeIndex[i].first;
        TPZFNMatrix<4> GradVi(fDimension,fDimension),GradVit(fDimension,fDimension),Dui(fDimension,fDimension);
        for (int e=0; e<fDimension; e++) {
            phiVi(e,0) = phiV(iphi,0)*datavec[vindex].fNormalVec(e,ivec);
            for (int f=0; f<fDimension; f++) {
                GradVi(e,f) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                //termo transposto:
                GradVit(f,e) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
            }
        }
        
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Dui(e,f)= 0.5 * (GradVi(e,f) + GradVit(e,f));
            }
        }
        
        
        if(this->HasForcingFunction()){
            this->ForcingFunction()->Execute(datavec[vindex].x, f);
        }
        
        
        
        STATE phi_dot_f = 0.0;
        for (int e=0; e<fDimension; e++) {
            phi_dot_f += phiVi(e)*f[e];
        }
        
        ef(i) += weight * phi_dot_f;
        
        
        
        // matrix A - gradV
        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            
            TPZFNMatrix<4> GradVj(fDimension,fDimension),GradVjt(fDimension,fDimension),Duj(fDimension,fDimension);
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradVj(e,f) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                    //termo transposto:
                    GradVjt(f,e) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Duj(e,f)= 0.5 * (GradVj(e,f) + GradVjt(e,f));
                }
            }
            
            
            STATE val = Inner(Dui, Duj);
            ek(i,j) += 2. * weight * fViscosity * val ; ///Visc*(GradU+GradU^T):GradPhi
            
        }
        
        
        
        // matrix B - pressure and velocity
        for (int j = 0; j < nshapeP; j++) {
            
            TPZManVector<REAL,3> GradPj(fDimension);
            for (int e=0; e<fDimension; e++) {
                GradPj[e] = dphiPx(e,j);
            }
            
            STATE fact = (-1.) * weight * phiP(j,0) * Tr( GradVi ); ///p*div(U)
            
            
            // colocar vectoriais vezes pressao
            // Matrix B
            ek(i, nshapeV+j) += fact;
            
            // colocar pressao vezes vectoriais
            // Matrix B^T
            ek(nshapeV+j,i) += fact;
            
        }//j
        
    }
    
    

    
}


void TPZStokesMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavec[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[vindex]);
    }
    // Setting forcing function
    /*STATE force = 0.;
     if(this->fForcingFunction) {
     TPZManVector<STATE> res(1);
     fForcingFunction->Execute(datavec[pindex].x,res);
     force = res[0];
     }*/
    
    //Gravity
    //    STATE rhoi = 900.; //itapopo
    //    STATE g = 9.81; //itapopo
    //    STATE force = rhoi*g;
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    //    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    //    //Normal
    //    TPZManVector<REAL,3> &normal = datavec[vindex].normal;
    
    // Getting the linear combination or finite element approximations
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    //Adaptação para Hdiv
    int ekr= ek.Rows();
    
    //Vefifica se HDiv
    if(ekr!=nshapeP+nshapeV){
        nshapeV=nshapeV/2;
    }
    
    
    int gy=v_h.size();
    
    
    TPZFNMatrix<9> phiVi(fDimension,1,0.),phiVni(1,1,0.), phiVj(fDimension,1,0.),phiVnj(1,1,0.), phiPi(fDimension,1),phiPj(fDimension,1);
    
    TPZFMatrix<STATE> v_2=bc.Val2();
    TPZFMatrix<STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    switch (bc.Type()) {
        case 0: //Dirichlet for continuous formulation
        {
            
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                bc.ForcingFunction()->Execute(datavec[vindex].x,vbc);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D = vbc[2];
            }
            
#ifdef IsHDivQ
            
            for(int i = 0; i < nshapeV; i++ )
            {
                
                //Adaptação para Hdiv
                
                TPZManVector<REAL> n = datavec[0].normal;
                
                REAL vh_n = v_h[0];
                REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                
                ef(i,0) += -weight * gBigNumber * (vh_n - v_n) * phiV(i,0);
                
                for(int j = 0; j < nshapeV; j++){
                    
                    ek(i,j) += weight * gBigNumber * phiV(j,0) * phiV(i,0);
                    
                }
                
            }
            
            
            
#else
            
            for(int i = 0; i < nshapeV; i++ )
            {
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*phiV(iphi,0);
                }
                
                
                //Adaptação para Hdiv
                
                STATE factef=0.0;
                for(int is=0; is<gy ; is++){
                    factef += v_2(is,0) * phiVi(is,0);
                }
                
                ef(i,0) += weight * fSigma * factef;
                
                //                    ef(i,0) += weight * gBigNumber * ( (v_h[0] - vx_D) * phiVi(0,0) + (v_h[1] - vy_D) * phiVi(1,0) );
                
                for(int j = 0; j < nshapeV; j++){
                    int jphi = datavec[vindex].fVecShapeIndex[j].second;
                    int jvec = datavec[vindex].fVecShapeIndex[j].first;
                    
                    
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*phiV(jphi,0);
                    }
                    
                    //Adaptação para Hdiv
                    
                    STATE factek = 0.0;
                    for(int is=0; is<gy ; is++){
                        factek += phiVj(is,0) * phiVi(is,0);
                    }
                    
                    ek(i,j) += weight * fSigma * factek;
                    
                    //                    ek(i,j) += weight * gBigNumber * (phiVj(0,0) * phiVi(0,0) + phiVj(1,0) * phiVi(1,0) );
                    
                    //std::cout<<ek<<std::endl;
                    
                }
                
            }
#endif
            
            //pressao
            
            //            for(int i = 0; i < nshapeP; i++ )
            //            {
            //
            //                for (int e=0; e<fDimension; e++) {
            //                    phiPi(e,0)=phiP(i,0);
            //                }
            //
            //                ef(i+nshapeV) += -weight * gBigNumber * ( (p_h[0] - p_D) * phiPi(0,0));
            //
            //                for(int j = 0; j < nshapeP; j++){
            //
            //                    for (int e=0; e<fDimension; e++) {
            //                        phiPj(e,0)=phiP(j,0);
            //                    }
            //
            //                    ek(i+nshapeV,j+nshapeV) += weight * gBigNumber * (phiPj(0,0) * phiPi(0,0));
            //
            //                }
            //
            //            }
            
            
            
        }
            break;
            
        case 1: //Neumann for continuous formulation
        {
            
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> nbc(4);
                bc.ForcingFunction()->Execute(datavec[vindex].x,nbc);
                v_1(0,0) = nbc[0];
                v_1(0,1) = nbc[1];
                v_1(1,0) = nbc[2];
                v_1(1,1) = nbc[3];
                
                //                p_D = nbc[2];
            }
            
            
            for(int i = 0; i < nshapeV; i++ )
            {
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*phiV(iphi,0);
                }
                
                TPZManVector<REAL> n = datavec[0].normal;
                
                TPZFNMatrix<9> pn(fDimension,1);
                
                
                for (int e=0; e<fDimension; e++) {
                    for (int f=0; f<fDimension; f++) {
                        pn(e,0)=n[f]*v_1(e,f);
                    }
                }
                
                //Adaptação para Hdiv
                
                STATE factef=0.0;
                for(int is=0; is<gy ; is++){
                    factef += (pn(is,0))* phiVi(is,0);
                }
                
                ef(i,0) += weight * factef;
                
            }
            
            
        }
            
            
            
            break;
            
        case 2: //Condição Mista
        {
            
            TPZFMatrix<STATE> v_2=bc.Val2();
            TPZFMatrix<STATE> v_1=bc.Val1();
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> vbc(2),pbc(1);
                bc.ForcingFunction()->Execute(datavec[vindex].x,vbc);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                v_1(0,0) = pbc[0];
                
            }
            
            
            
            //            STATE vx_D = bc.Val2()(0,0);
            //            STATE vy_D = bc.Val2()(1,0);
            
            for(int i = 0; i < ekr; i++ )
            {
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*phiV(iphi,0);
                }
                
                STATE factef=0.0;
                for(int is=0; is<gy ; is++){
                    factef += (v_2(is,0)) * phiVi(is,0);
                }
                
                ef(i,0) += weight * factef;
                
                //                    ef(i,0) += weight * gBigNumber * ( (v_h[0] - vx_D) * phiVi(0,0) + (v_h[1] - vy_D) * phiVi(1,0) );
                
                for(int j = 0; j < ekr; j++){
                    int jphi = datavec[vindex].fVecShapeIndex[j].second;
                    int jvec = datavec[vindex].fVecShapeIndex[j].first;
                    
                    
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*phiV(jphi,0);
                    }
                    
                    //Adaptação para Hdiv
                    
                    STATE factek=0.0;
                    for(int is=0; is<gy ; is++){
                        factek += (v_1(is,0))* phiVj(is,0) * phiVi(is,0);
                    }
                    
                    ek(i,j) += weight * factek;
                    
                    //                    ek(i,j) += weight * gBigNumber * (phiVj(0,0) * phiVi(0,0) + phiVj(1,0) * phiVi(1,0) );
                    
                }
                
            }
            
        }
            break;
            
            
        case 3: //Contribuicao ponto no x
        {
            
            REAL p_D = v_2(0,0);
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> pbc(1);
                bc.ForcingFunction()->Execute(datavec[vindex].x,pbc);
                p_D = pbc[0];
                
            }
            
            TPZManVector<REAL> n = datavec[0].normal;
            
            REAL phiVi_n;
            
            
            for(int i = 0; i < 3; i++ )
            {
                phiVi_n = 0.0;
                
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*phiV(iphi,0);
                    phiVi_n += phiVi(e,0)*n[e];
                }
                
                ef(i*2,0) += -1.0*weight * p_D * phiVi_n;
                
                
            }
            
            
        }
            break;
            
        case 4: //Contribuicao ponto no y
        {
            
            REAL p_D = v_2(0,0);
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> pbc(1);
                bc.ForcingFunction()->Execute(datavec[vindex].x,pbc);
                p_D = pbc[0];
                
            }
            
            TPZManVector<REAL> n = datavec[0].normal;
            
            REAL phiVi_n;
            
            
            for(int i = 0; i < 3; i++ )
            {
                phiVi_n = 0.0;
                
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*phiV(iphi,0);
                    phiVi_n += phiVi(e,0)*n[e];
                }
                
                ef(i*2+1,0) += -1.0*weight * p_D * phiVi_n;
                
                
            }
            
            
        }
            break;
            
        case 5: //Ponto pressao
        {
            p_D = bc.Val2()(0,0);
            
            
            for(int i = 0; i < nshapeP; i++ )
            {
                
                //                for (int e=0; e<fDimension; e++) {
                //                    phiPi(e,0)=phiP(i,0);
                //                }
                
                ef(i) += 1.0 * p_D * phiP(i,0);
                
                for(int j = 0; j < nshapeP; j++){
                    
                    //                    for (int e=0; e<fDimension; e++) {
                    //                        phiPj(e,0)=phiP(j,0);
                    //                    }
                    
                    ek(i,j) += 1.0 * (phiP(i,0) * phiP(j,0));
                    
                }
                
            }
            
        }
            break;
            

            
        default:
        {
            std::cout << "Boundary not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    

}




////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    
    // Verificar que
    // os termos mistos devem estar sem viscosidade!
    
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space for datavecleft
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    // Setting forcing function
    /*STATE force = 0.;
     if(this->fForcingFunction) {
     TPZManVector<STATE> res(1);
     fForcingFunction->Execute(datavec[pindex].x,res);
     force = res[0];
     }*/
    
    //Gravity
    STATE rhoi = 900.; //itapopo
    STATE g = 9.81; //itapopo
    STATE force = rhoi*g;
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    // V - right
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    // P - right
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    //Detjac
    REAL Detjac=fabs(data.detjac);
    
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);
    
    //TPZManVector<REAL,3> normalx(fDimension,phiP2.Cols());
    //TPZAxesTools<REAL>::Axes2XYZ(normal, normalx, data.axes);
    
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    
    for(int i1 = 0; i1 < nshapeV1; i1++ )
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
        
        
        TPZFNMatrix<9> GradV1ni(fDimension,1,0.),phiV1i(fDimension,1),phiV1ni(1,1,0.);
        TPZFNMatrix<4> GradV1i(fDimension,fDimension,0.),GradV1it(fDimension,fDimension,0.),Du1i(fDimension,fDimension,0.),Du1ni(fDimension,1,0.);
        
        for (int e=0; e<fDimension; e++) {
            
            phiV1i(e,0)=datavecleft[vindex].fNormalVec(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
            phiV1ni(0,0)+=phiV1i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV1i(e,f) = datavecleft[vindex].fNormalVec(e,ivec1)*dphiVx1(f,iphi1);
                //termo transposto:
                GradV1it(f,e) = datavecleft[vindex].fNormalVec(e,ivec1)*dphiVx1(f,iphi1);
                //GradV1ni(e,0)+=datavecleft[vindex].fNormalVec(e,ivec1)*dphiVx1(f,iphi1)*normal[f];
            }
        }
        
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du1i(e,f)= (1./2.) * (GradV1i(e,f) + GradV1it(e,f));
            }
        }
        
        //Du1ni
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du1ni(e,0) += Du1i(e,f)*normal[f] ;
            }
        }
        
        
        
        
        TPZFNMatrix<9> GradV1nj(fDimension,1,0.),phiV1j(fDimension,1),phiV1nj(1,1,0.);
        
        // K11 - (trial V left) * (test V left)
        for(int j1 = 0; j1 < nshapeV1; j1++){
            int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
            
            TPZFNMatrix<4> GradV1j(fDimension,fDimension,0.),GradV1jt(fDimension,fDimension,0.),Du1j(fDimension,fDimension,0.),Du1nj(fDimension,1,0.);
            
            for (int e=0; e<fDimension; e++) {
                
                phiV1j(e,0)=datavecleft[vindex].fNormalVec(e,jvec1)*datavecleft[vindex].phi(jphi1,0);
                phiV1nj(0,0)+=phiV1j(e,0)*normal[e];
                
                for (int f=0; f<fDimension; f++) {
                    GradV1j(e,f) = datavecleft[vindex].fNormalVec(e,jvec1)*dphiVx1(f,jphi1);
                    //termo transposto:
                    GradV1jt(f,e) = datavecleft[vindex].fNormalVec(e,jvec1)*dphiVx1(f,jphi1);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1j(e,f)= (1./2.) * (GradV1j(e,f) + GradV1jt(e,f));
                }
            }
            
            //Du1nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1nj(e,0) += Du1j(e,f)*normal[f] ;
                }
            }
            
            STATE fact = (-1./2.) * weight * 2.* fViscosity * InnerVec(phiV1i, Du1nj);
            
            ek(i1,j1) +=fact;
            ek(j1,i1) +=-fact*fTheta;
            
            
            //Penalidade:
            
            STATE penalty = fSigma * weight * fViscosity * InnerVec(phiV1i, phiV1j);
            ek(i1,j1) +=penalty;
            
            
        }
        
        // K12 e K21 - (trial V left) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            
            STATE fact = (1./2.) * weight * Inner(phiV1ni,phiP1j);
            
            ek(i1,j1+nshapeV1) += fact;
            ek(j1+nshapeV1,i1) += fact;
            
        }
        
        
        // K13 - (trial V left) * (test V right)
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            TPZFNMatrix<9> GradV2nj(fDimension,1),phiV2j(fDimension,1),phiV2nj(1,1,0.);
            //TPZManVector<REAL,3> phiP1j(fDimension);
            
            TPZFNMatrix<4> GradV2j(fDimension,fDimension,0.),GradV2jt(fDimension,fDimension,0.),Du2j(fDimension,fDimension,0.),Du2nj(fDimension,1,0.);
            
            for (int e=0; e<fDimension; e++) {

                phiV2j(e,0)=datavecright[vindex].fNormalVec(e,jvec2)*datavecright[vindex].phi(jphi2,0);
                phiV2nj(0,0)+=phiV2j(e,0)*normal[e];
                
                for (int f=0; f<fDimension; f++) {
                    GradV2j(e,f) = datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2);
                    //termo transposto:
                    GradV2jt(f,e) = datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2j(e,f)= (1./2.) * (GradV2j(e,f) + GradV2jt(e,f));
                }
            }
            
            //Du2nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2nj(e,0) += Du2j(e,f)*normal[f] ;
                }
            }
            
            
            STATE fact = (-1./2.) * weight * 2. * fViscosity * InnerVec(phiV1i,Du2nj);
            
            ek(i1,j2+nshapeV1+nshapeP1) += fact;
            ek(j2+nshapeV1+nshapeP1,i1) += -fact*fTheta;
            
            //Penalidade:
            
            STATE penalty = fSigma * weight * fViscosity * InnerVec(phiV1i, phiV2j);
            ek(i1,j2+nshapeV1+nshapeP1) += -penalty;
            
            
        }
        
        // K14 e K41 - (trial V left) * (test P right)
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            STATE fact = (1./2.) * weight * InnerVec(phiV1ni,phiP2j);
            
            ek(i1,j2+2*nshapeV1+nshapeP1) += fact;
            ek(j2+2*nshapeV1+nshapeP1,i1) += fact;
            
        }
        
    }
    
    
    for(int i2 = 0; i2 < nshapeV2; i2++ ){
        
        TPZFNMatrix<9> GradV2ni(fDimension,1),phiV2i(fDimension,1),phiV2ni(1,1,0.);
        TPZFNMatrix<4> GradV2i(fDimension,fDimension,0.),GradV2it(fDimension,fDimension,0.),Du2i(fDimension,fDimension,0.),Du2ni(fDimension,1,0.);
        
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        
        for (int e=0; e<fDimension; e++) {
            
            phiV2i(e,0)=datavecright[vindex].fNormalVec(e,ivec2)*datavecright[vindex].phi(iphi2,0);
            phiV2ni(0,0)+=phiV2i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV2i(e,f) = datavecright[vindex].fNormalVec(e,ivec2)*dphiVx2(f,iphi2);
                //termo transposto:
                GradV2it(f,e) = datavecright[vindex].fNormalVec(e,ivec2)*dphiVx2(f,iphi2);
            }
        }
        
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du2i(e,f)= (1./2.) * (GradV2i(e,f) + GradV2it(e,f));
            }
        }
        
        //Du2ni
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du2ni(e,0) += Du2i(e,f)*normal[f] ;
            }
        }
        
        
        
        // K31 - (trial V right) * (test V left)
        for(int j1 = 0; j1 < nshapeV1; j1++){
            int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
            
            TPZFNMatrix<4> GradV1j(fDimension,fDimension,0.),GradV1jt(fDimension,fDimension,0.),Du1j(fDimension,fDimension,0.),Du1nj(fDimension,1,0.);

            TPZFNMatrix<9> phiV1j(fDimension,1),phiV1nj(1,1,0.);
            
            for (int e=0; e<fDimension; e++) {
            
                phiV1j(e,0)=datavecleft[vindex].fNormalVec(e,jvec1)*datavecleft[vindex].phi(jphi1,0);
                phiV1nj(0,0)+=phiV1j(e,0)*normal[e];
                
                for (int f=0; f<fDimension; f++) {
                    GradV1j(e,f) = datavecleft[vindex].fNormalVec(e,jvec1)*dphiVx1(f,jphi1);
                    //termo transposto:
                    GradV1jt(f,e) = datavecleft[vindex].fNormalVec(e,jvec1)*dphiVx1(f,jphi1);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1j(e,f)= (1./2.) * (GradV1j(e,f) + GradV1jt(e,f));
                }
            }
            
            //Du1nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1nj(e,0) += Du1j(e,f)*normal[f] ;
                }
            }
            
            
            
            STATE fact = (1./2.) * weight * 2. * fViscosity * InnerVec(phiV2i, Du1nj);
            
            ek(i2+nshapeV1+nshapeP1,j1) += fact;
            ek(j1,i2+nshapeV1+nshapeP1) += -fact*fTheta;
            
            //Penalidade:
            
            STATE penalty = fSigma * weight * fViscosity * InnerVec(phiV2i, phiV1j);
            ek(i2+nshapeV1+nshapeP1,j1) += -penalty;


        }
        
        // K32 e K23 - (trial V right) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            STATE fact = (-1./2.) * weight * InnerVec(phiV2ni,phiP1j);
            
            ek(i2+nshapeV1+nshapeP1,j1+nshapeV1) += fact;
            ek(j1+nshapeV1,i2+nshapeV1+nshapeP1) += fact;
            
        }
        
        
        // K33 - (trial V right) * (test V right)
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            TPZFNMatrix<9> GradV2nj(fDimension,1);
            //TPZManVector<REAL,3> phiP1j(fDimension);
            
            TPZFNMatrix<4> GradV2j(fDimension,fDimension,0.),GradV2jt(fDimension,fDimension,0.),Du2j(fDimension,fDimension,0.),Du2nj(fDimension,1,0.);
            
            TPZFNMatrix<9> phiV2j(fDimension,1),phiV2nj(1,1,0.);
            
            
            for (int e=0; e<fDimension; e++) {

                phiV2j(e,0)=datavecright[vindex].fNormalVec(e,jvec2)*datavecright[vindex].phi(jphi2,0);
                phiV2nj(0,0)+=phiV2j(e,0)*normal[e];
                
                for (int f=0; f<fDimension; f++) {
                    GradV2j(e,f) = datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2);
                    //termo transposto:
                    GradV2jt(f,e) = datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2j(e,f)= (1./2.) * (GradV2j(e,f) + GradV2jt(e,f));
                }
            }
            
            //Du2nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2nj(e,0) += Du2j(e,f)*normal[f] ;
                }
            }
            
            
            STATE fact = (1./2.) * weight *  2. * fViscosity * InnerVec(phiV2i,Du2nj);
            
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += fact;
            ek(j2+nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += -fact*fTheta;
            
            //Penalidade:
            
            STATE penalty = fSigma * weight * fViscosity * InnerVec(phiV2i, phiV2j);
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) +=penalty;
            
            
        }
        
        // K34 e K43- (trial V right) * (test P right)
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            STATE fact = (-1./2.) * weight * InnerVec(phiV2ni,phiP2j);
            
            ek(i2+nshapeV1+nshapeP1,j2+2*nshapeV1+nshapeP1) += fact;
            ek(j2+2*nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact;
        }
        
    }
    
    
    
}


void TPZStokesMaterial::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    

    
#ifdef IsH1
    //Caso H1 -> return
    return;
#endif
    
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavec[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[vindex]);
    }
    
    
    //Gravity
    STATE rhoi = 900.; //itapopo
    STATE g = 9.81; //itapopo
    STATE force = rhoi*g;
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    
    int sizek=ek.Rows();
    
    //Dirichlet
    
    TPZFMatrix<STATE> v_2=bc.Val2();
    TPZFMatrix<STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> vbc(3);
        bc.ForcingFunction()->Execute(datavec[vindex].x,vbc);
        v_2(0,0) = vbc[0];
        v_2(1,0) = vbc[1];
        p_D=vbc[2];
        
    }
    
    
    
    for(int i = 0; i < nshapeV; i++ )
    {
        int iphi = datavec[vindex].fVecShapeIndex[i].second;
        int ivec = datavec[vindex].fVecShapeIndex[i].first;
        TPZFNMatrix<9> GradVni(fDimension,1,0.),phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
        GradVni.Zero();
        
        TPZFNMatrix<4> GradVi(fDimension,fDimension,0.),GradVit(fDimension,fDimension,0.),Dui(fDimension,fDimension,0.),Duni(fDimension,1,0.);
        
        for (int e=0; e<fDimension; e++) {
            
            for (int f=0; f<fDimension; f++) {
                GradVi(e,f) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                //termo transposto:
                GradVit(f,e) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);
                
            }
        }
        
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Dui(e,f)= (1./2.) * (GradVi(e,f) + GradVit(e,f));
            }
        }
        
        //Duni
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Duni(e,0) += Dui(e,f)*normal[f] ;
            }
        }
        
        //GradVni
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                GradVni(e,0) += GradVi(e,f)*normal[f] ;
            }
        }
        
        
        for (int e=0; e<fDimension; e++) {
            phiVi(e,0)=datavec[vindex].fNormalVec(e,ivec)*datavec[vindex].phi(iphi,0);
            phiVni(0,0)+=phiVi(e,0)*normal[e];
            
        }
        
        TPZManVector<REAL> n = data.normal;
        TPZManVector<REAL> t(2);
        t[0]=-n[1];
        t[1]=n[0];
        
        
        
        phiVti(0,0)= t[0] * phiVi(0,0) + t[1] * phiVi(1,0);
        TPZFNMatrix<9> phiVtit(fDimension,1,0.);
        phiVtit(0,0)=phiVti(0,0)*t[0];
        phiVtit(1,0)=phiVti(0,0)*t[1];
        
        TPZFNMatrix<9> phiVnin(fDimension,1,0.);
        phiVnin(0,0)=phiVni(0,0)*n[0];
        phiVnin(1,0)=phiVni(0,0)*n[1];
        
        
        
#ifdef IsHDivQ
        
        REAL vh_t = v_h[1];
        
        REAL v_t = t[0] * v_2[0] + t[1] * v_2[1];
        
        TPZManVector<REAL> v_tt(2);
        v_tt[0]=v_t*t[0];
        v_tt[1]=v_t*t[1];
        
        TPZManVector<REAL> vh_tt(2);
        vh_tt[0]=vh_t*t[0];
        vh_tt[1]=vh_t*t[1];
        
        TPZFNMatrix<9> diffvt(fDimension,1,0.);
        diffvt(0,0)=v_tt[0];
        diffvt(1,0)=v_tt[1];
        
        STATE factf=(-1.) * weight * fViscosity * InnerVec(diffvt,GradVni) ;
        
        ef(i,0) += fTheta*factf;
        
        
        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            
            TPZFNMatrix<9> GradVnj(fDimension,1),phiVtj(1,1,0.),phiVj(fDimension,1);
            
            for (int e=0; e<fDimension; e++) {
                phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*datavec[vindex].phi(jphi,0);
            }
            
            //std::cout<<phiVj<<std::endl;
            
            phiVtj(0,0)= t[0] * phiVj(0,0) + t[1] * phiVj(1,0);
            
            
            
            TPZFNMatrix<4> GradVj(fDimension,fDimension,0.),GradVjt(fDimension,fDimension,0.),Duj(fDimension,fDimension,0.),Dunj(fDimension,1,0.);
            
            for (int e=0; e<fDimension; e++) {
                
                for (int f=0; f<fDimension; f++) {
                    GradVj(e,f) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                    //termo transposto:
                    GradVjt(f,e) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Duj(e,f)= (1./2.) * (GradVj(e,f) + GradVjt(e,f));
                }
            }
            
            //Du2nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Dunj(e,0) += Duj(e,f)*normal[f] ;
                }
            }
            
            
            
            STATE fact =(-1.) * weight * 2. * fViscosity * InnerVec(phiVtit, Dunj) ;
            ek(i,j) += fact ;
            ek(j,i) += fTheta*fact;
            
        }
        
        
        
#endif
        
        
        
        
#ifdef IsDGM
        
        
        
        if(bc.HasForcingFunction())
        {
            TPZManVector<STATE> vbc(3);
            bc.ForcingFunction()->Execute(datavec[vindex].x,vbc);
            v_2(0,0) = vbc[0];
            v_2(1,0) = vbc[1];
            p_D=vbc[2];
            
        }
        
//        REAL vh_t = v_h[1];
//        
//        REAL v_t = t[0] * v_2[0] + t[1] * v_2[1];
        
//        TPZManVector<REAL> v_tt(2);
//        v_tt[0]=v_t*t[0];
//        v_tt[1]=v_t*t[1];
//        
//        TPZManVector<REAL> vh_tt(2);
//        vh_tt[0]=vh_t*t[0];
//        vh_tt[1]=vh_t*t[1];
//        
//        TPZFNMatrix<9> diffvt(fDimension,1,0.);
//        diffvt(0,0)=v_tt[0];
//        diffvt(1,0)=v_tt[1];
//        REAL vh_n = v_h[0];
//        
//        REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
//        
//        TPZManVector<REAL> v_nn(2);
//        v_nn[0]=v_n*n[0];
//        v_nn[1]=v_n*n[1];
//        
//        TPZManVector<REAL> vh_nn(2);
//        vh_nn[0]=vh_n*n[0];
//        vh_nn[1]=vh_n*n[1];
//        
//        TPZFNMatrix<9> diffvn(fDimension,1,0.);
//        diffvn(0,0)=v_nn[0];
//        diffvn(1,0)=v_nn[1];
//        
//        
//        STATE factfn= 2. * weight * fViscosity * InnerVec(diffvn, Duni) ;
//        
//        ef(i,0) += fTheta*factfn ;

        STATE fact= 2. * weight * fViscosity * InnerVec(v_2, Duni) ;
        
        ef(i,0) += fTheta*fact;

        
        
        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            
            TPZFNMatrix<9> GradVnj(fDimension,1,0.),phiVj(fDimension,1);
            TPZFNMatrix<4> GradVj(fDimension,fDimension,0.),GradVjt(fDimension,fDimension,0.),Duj(fDimension,fDimension,0.),Dunj(fDimension,1,0.);

            
            for (int e=0; e<fDimension; e++) {
                phiVj(e,0)=datavec[vindex].fNormalVec(e,jvec)*datavec[vindex].phi(jphi,0);
                
                for (int f=0; f<fDimension; f++) {
                    GradVj(e,f) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                    //termo transposto:
                    GradVjt(f,e) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Duj(e,f)= (1./2.) * (GradVj(e,f) + GradVjt(e,f));
                }
            }
            
            //Dunj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Dunj(e,0) += Duj(e,f)*normal[f] ;
                }
            }
            
            
            //GradVni
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradVnj(e,0) += GradVj(e,f)*normal[f] ;
                }
            }
            
            STATE fact1=(-2.) * weight * fViscosity * InnerVec(phiVi, Dunj) ;
            
            
            ek(i,j) += fact1 ;
            ek(j,i) += -fTheta*fact1;
            
     
            
            
//            STATE fact2=(-2.) * weight * fViscosity * InnerVec(phiVnin, Dunj) ;
//            
//            
//            ek(i,j) += fact2 ;
//            ek(j,i) += fTheta*fact2;
//            
//            STATE fact1=(-2.) * weight * fViscosity * InnerVec(phiVtit, Dunj) ;
//            
//            
//            ek(i,j) += fact1 ;
//            ek(j,i) += fTheta*fact1;
            
            
//            STATE fact2=(-2.) * weight * fViscosity * InnerVec(phiVnin, Dunj) ;
//            
//            
//            ek(i,j) += fact2 ;
//            ek(j,i) += fTheta*fact2;
            
            
            //std::cout<<ek<<std::endl;
            
        }
        
#endif
        
        
        //pressao fracamente
        
        //papapa
        
        
                // K12 e K21 - (trial V left) * (test P left)
                for(int j = 0; j < nshapeP; j++){
        
        
                    TPZFNMatrix<9> phiPj(1,1,0.),v_2n(1,1,0.);;
                    phiPj(0,0)=phiP(j,0);
        
        
                    for (int e=0; e<fDimension; e++) {
                        v_2n(0,0)+=v_2(e,0)*normal[e];
                    }
        
                   // std::cout<<v_2<<std::endl;
                   // std::cout<<normal<<std::endl;
        
        
                    //STATE factfp= (-1.) * weight * fViscosity * Inner(v_2n,phiPj);
        
                    //ef(j+nshapeV,0) += factfp ;
        
                    //Condição em f2:
                    
                    STATE factfp= (1.) * weight * fViscosity * v_2n(0,0)* phiPj(0,0);
        
                    ef(j+nshapeV,0) += factfp ;
        
        
                    STATE fact = (1.) * weight * Inner(phiVni,phiPj);
                    ek(i,j+nshapeV) += fact;
                    ek(j+nshapeV,i) += fact;
        //
        
        
                }
        
        
        
        
        
        
    }
    
    
    
    
}



////////////////////////////////////////////////////////////////////

STATE TPZStokesMaterial::Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
    
    //inner product of two tensors
    
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
    
    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0;
    
    for(int i = 0; i < S.Cols(); i++){
        for(int j = 0; j < S.Cols(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}


////////////////////////////////////////////////////////////////////

STATE TPZStokesMaterial::InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
    
    //inner product of two vectors
    
    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0;
    
    for(int j = 0; j < S.Cols(); j++){
        for(int i = 0; i < S.Rows(); i++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}



////////////////////////////////////////////////////////////////////

STATE TPZStokesMaterial::Tr( TPZFMatrix<REAL> &GradU ){
    
#ifdef DEBUG
    if( GradU.Rows() != GradU.Cols() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0.;
    
    for(int i = 0; i < GradU.Rows(); i++){
        Val += GradU(i,i);
    }
    
    return Val;
}


/// transform a H1 data structure to a vector data structure
void TPZStokesMaterial::FillVecShapeIndex(TPZMaterialData &data)
{
    data.fNormalVec.Resize(fDimension,fDimension);
    data.fNormalVec.Identity();
    data.fVecShapeIndex.Resize(fDimension*data.phi.Rows());
    for (int d=0; d<fDimension; d++) {
        for (int i=0; i<data.phi.Rows(); i++) {
            data.fVecShapeIndex[i*fDimension+d].first = d;
            data.fVecShapeIndex[i*fDimension+d].second = i;
        }
    }
}



void TPZStokesMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    
    // @omar:: ate this point just velocity norms
    
    //                             TPZVec<REAL> &x,TPZVec<STATE> &u,
    //                             TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &/*flux*/,
    //                             TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    TPZManVector<STATE> Velocity, Pressure;
    Velocity.Fill(0.0);
    Pressure.Fill(0.0);
    
    this->Solution(data,VariableIndex("V"), Velocity);
    this->Solution(data,VariableIndex("P"), Pressure);
    
    int vindex = this->VIndex();
    int pindex = this->PIndex();
    
    TPZFMatrix<REAL> dudx(Dimension(),Dimension());
    TPZFMatrix<REAL> &dsol = data[vindex].dsol[0];
    TPZFMatrix<REAL> &dsolp = data[pindex].dsol[0];
    //std::cout<<dsol<<std::endl;
    
    //Adaptação feita para Hdiv
    dsol.Resize(Dimension(),Dimension());
    
    TPZFNMatrix<2,STATE> dsolxy(2,2), dsolxyp(2,1);
    TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, data[vindex].axes);
    TPZAxesTools<STATE>::Axes2XYZ(dsolp, dsolxyp, data[pindex].axes);
    
    
    
    int shift = 3;
    // velocity
    
    //values[2] : erro norma L2
    REAL diff, diffp;
    errors[1] = 0.;
    for(int i=0; i<Dimension(); i++) {
        diff = Velocity[i] - u_exact[i];
        errors[1]  += diff*diff;
    }
    
    ////////////////////////////////////////////////// H1 / GD
    
#ifdef IsH1
    
    //values[2] : erro em semi norma H1
    errors[2] = 0.;
    TPZFMatrix<STATE> S(Dimension(),Dimension(),0.0);
    for(int i=0; i<Dimension(); i++) {
        for(int j=0; j<Dimension(); j++) {
            S(i,j) = dsolxy(i,j) - du_exact(i,j);
        }
    }
    
    diff = Inner(S, S);
    errors[2]  += diff;
    
    //values[0] : erro em norma H1 <=> norma Energia
    errors[0]  = errors[1]+errors[2];
    
#endif
    
#ifdef IsDGM
    
    //values[2] : erro em semi norma H1
    errors[2] = 0.;
    TPZFMatrix<STATE> S(Dimension(),Dimension(),0.0);
    for(int i=0; i<Dimension(); i++) {
        for(int j=0; j<Dimension(); j++) {
            S(i,j) = dsolxy(i,j) - du_exact(i,j);
        }
    }
    
    diff = Inner(S, S);
    errors[2]  += diff;
    
    //values[0] : erro em norma H1 <=> norma Energia
    errors[0]  = errors[1]+errors[2];
    
#endif
    
    
    ////////////////////////////////////////////////// H1 / GD
    
    // pressure
    
    /// values[1] : eror em norma L2
    diffp = Pressure[0] - u_exact[2];
    errors[shift+1]  = diffp*diffp;
    
    // pressure gradient error ....
    
    errors[shift+2] = 0.;
    TPZFMatrix<STATE> Sp(Dimension(),1,0.0);
    for(int i=0; i<Dimension(); i++) {
        Sp(i,0) = dsolxyp(i,0) - du_exact(2,i);
    }
    
    diffp = InnerVec(Sp, Sp);
    errors[shift+2]  += diffp;
    
    //values[0] : erro em norma H1 <=> norma Energia
    errors[shift]  = errors[1+shift]+errors[2+shift];
    
    ////////////////////////////////////////////////// HDIV
    
#ifdef IsHDivQ
    /// erro norma HDiv
    
    STATE Div_exact=0., Div=0.;
    for(int i=0; i<Dimension(); i++) {
        Div_exact+=du_exact(i,i);
        Div+=dsolxy(i,i);
    }
    
    diff = Div-Div_exact;
    
    errors[2]  = diff*diff;
    
    errors[0]  = errors[1]+errors[2];
    
#endif
    
    ////////////////////////////////////////////////// HDIV
    
}
