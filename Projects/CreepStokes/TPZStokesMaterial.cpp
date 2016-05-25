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


TPZStokesMaterial::TPZStokesMaterial() : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(){
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    
}

////////////////////////////////////////////////////////////////////

TPZStokesMaterial::TPZStokesMaterial(int matid, int dimension, REAL viscosity) : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(matid),fViscosity(viscosity), fDimension(dimension)
{
    // symmetric version
    fTheta = -1;
    
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);

}

////////////////////////////////////////////////////////////////////

TPZStokesMaterial::TPZStokesMaterial(const TPZStokesMaterial &mat) : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(mat), fViscosity(mat.fViscosity),fDimension(mat.fDimension), fTheta(mat.fTheta)
{

    
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
    
    if (!strcmp("Pressure", name.c_str())) return this->PIndex();
    if (!strcmp("Velocity", name.c_str())) return this->VIndex();
   
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
    
    int Vblock = this->VIndex();
    int Pblock = this->PIndex();
    
    TPZManVector<REAL,3> V = datavec[Vblock].sol[0];
    
    REAL P = datavec[Pblock].sol[0][0];
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        
        case 0: //Pressure
        {
            Solout[0] = P;
        }
            break;
        
        case 1: //Velocity
        {
            Solout[0] = V[0]; // Vx
            Solout[1] = V[1]; // Vy
            Solout[2] = V[2]; // Vy
        }
            break;
        
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
    
    

    
    for(int i = 0; i < nshapeV; i++ )
    {
        int iphi = datavec[vindex].fVecShapeIndex[i].second;
        int ivec = datavec[vindex].fVecShapeIndex[i].first;
        TPZFNMatrix<9> GradVi(fDimension,fDimension);
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                GradVi(e,f) = datavec[vindex].fNormalVec(e,ivec)*dphiVx(f,iphi);

            }
        }
        
        // matrix A - gradV
        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            TPZFNMatrix<9> GradVj(fDimension,fDimension);
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradVj(e,f) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi);
                }
            }
            
            ek(i,j) += weight * fViscosity * Inner(GradVi, GradVj) ; ///Visc*(GradU+GradU^T):GradPhi

            
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
    
    
    for (int ipressure = 0; ipressure < nshapeP; ipressure++) {
        TPZManVector<REAL,3> GradPi(fDimension);
        for (int e=0; e<fDimension; e++) {
            GradPi[e] = dphiPx(e,ipressure);
        }

        for (int jpressure = 0; jpressure < nshapeP; jpressure++) {
            // colocar as contribuicoes pressao - pressao aqui
            TPZManVector<REAL,3> GradPj(fDimension);
            for (int e=0; e<fDimension; e++) {
                GradPj[e] = dphiPx(e,jpressure);
            }
            // colocar os termos pressao pressao
             ek(nshapeV+ipressure, nshapeV+jpressure) += 0.;
            // talvez aqui nao tem nada???

        }
    }
    
    //teste: Zerar contribuições dos elementos
//    
//    for(int i=0;i<nshapeV+nshapeP;i++){
//        for(int j=0;j<nshapeV+nshapeP;j++){
//            ek(i,j)*=0.0;
//        }
//    
//    }
    
    
    
    std::cout<<ek<<std::endl;

    

}


void TPZStokesMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    

    DebugStop();
    
//#ifdef PZDEBUG
//    //2 = 1 Vel space + 1 Press space
//    int nrefl =  datavecleft.size();
//    int nrefr =  datavecright.size();
//    if (nrefl != 2 || nrefr != 2) {
//        std::cout << " Erro. The size of the datavec is different from 2 \n";
//        DebugStop();
//    }
//#endif
//    
//    
//    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
//        FillVecShapeIndex(datavecleft[vindex]);
//    }
//    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
//        FillVecShapeIndex(datavecright[vindex]);
//    }
    

    
    
}




////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){

    
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
    REAL Detjac=fabs(datavecleft[0].detjac);
    
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);


    
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
        
        
        for (int e=0; e<fDimension; e++) {
            
            phiV1i(e,0)=datavecleft[vindex].fNormalVec(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
            phiV1ni(0,0)+=phiV1i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV1ni(e,0)+=datavecleft[vindex].fNormalVec(e,ivec1)*dphiVx1(f,iphi1)*normal[f];
            }
        }
        
        std::cout<<phiV1i<<std::endl;
        std::cout<<phiV1ni<<std::endl;
        std::cout<<GradV1ni<<std::endl;
        
        TPZFNMatrix<9> GradV1nj(fDimension,1,0.);
        
        // K11 - (trial V left) * (test V left)
        for(int j1 = 0; j1 < nshapeV1; j1++){
            int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
            
            GradV1nj.Zero();
            
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradV1nj(e,0) += datavecleft[vindex].fNormalVec(e,jvec1)*dphiVx1(f,jphi1)*normal[f];
                }
            }
            
            ek(i1,j1) += (-1./2.) * weight * fViscosity * InnerVec(phiV1i, GradV1nj)*Detjac*0 ;
            
        
            
        }
        
        // K12 e K21 - (trial V left) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            

            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            STATE fact = (1./2.) * weight * fViscosity * Inner(phiV1ni,phiP1j);
            
            ek(i1,j1+nshapeV1) += fact;
            ek(j1+nshapeV1,i1) += fact;
         
        }
     
        
        // K13 - (trial V left) * (test V right)
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            TPZFNMatrix<9> GradV2nj(fDimension,1);
            //TPZManVector<REAL,3> phiP1j(fDimension);
            
            GradV2nj.Zero();
            
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradV2nj(e,0) += datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2)*normal[f];
                }
            }

            ek(i1,j2+nshapeV1+nshapeP1) += (-1./2.) * weight * fViscosity * InnerVec(phiV1i,GradV2nj)*Detjac*0;
            
        }

        // K14 e K41 - (trial V left) * (test P right)
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            STATE fact = (1./2.) * weight * fViscosity * InnerVec(phiV1ni,phiP2j)*Detjac*0.;
            
            ek(i1,j2+2*nshapeV1+nshapeP1) += fact;
            ek(j2+2*nshapeV1+nshapeP1,i1) += fact;
            
        }
       
    }

    
    for(int i2 = 0; i2 < nshapeV2; i2++ ){
        
        TPZFNMatrix<9> GradV2ni(fDimension,1),phiV2i(fDimension,1),phiV2ni(1,1,0.);
        
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        
        for (int e=0; e<fDimension; e++) {
            
            phiV2i(e,0)=datavecright[vindex].fNormalVec(e,ivec2)*datavecright[vindex].phi(iphi2,0);
            phiV2ni(0,0)+=phiV2i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV2ni(e,0) += datavecright[vindex].fNormalVec(e,ivec2)*dphiVx2(f,iphi2)*normal[f];
            }
        }
        
        // K31 - (trial V right) * (test V left)
        for(int j1 = 0; j1 < nshapeV1; j1++){
            int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
            TPZFNMatrix<9> GradV1nj(fDimension,1);
            GradV1nj.Zero();
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradV1nj(e,0) += datavecleft[vindex].fNormalVec(e,jvec1)*dphiVx1(f,jphi1)*normal[f];
                }
            }
            
            ek(i2+nshapeV1+nshapeP1,j1) += (1./2.) * weight * fViscosity * InnerVec(phiV2i, GradV1nj)*Detjac*0 ;
            
        }
        
        // K32 - (trial V right) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            STATE fact = (-1./2.) * weight * fViscosity * InnerVec(phiV2ni,phiP1j)*0;
            
            ek(i2+nshapeV1+nshapeP1,j1+nshapeV1) += fact;
            ek(j1+nshapeV1,i2+nshapeV1+nshapeP1) += fact;
            
        }
        
        
        // K33 - (trial V right) * (test V right)
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            TPZFNMatrix<9> GradV2nj(fDimension,1);
            //TPZManVector<REAL,3> phiP1j(fDimension);
            
            for (int e=0; e<fDimension; e++) {
                GradV2nj.Zero();
                for (int e=0; e<fDimension; e++) {
                    for (int f=0; f<fDimension; f++) {
                        GradV2nj(e,0) = datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2)*normal[f];
                    }
            }
            
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += (1./2.) * weight * fViscosity * InnerVec(phiV2i,GradV2nj);
            
        }
        
        // K34 - (trial V right) * (test P right)
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            STATE fact = (-1./2.) * weight * fViscosity * InnerVec(phiV2ni,phiP2j)*0;

            ek(i2+nshapeV1+nshapeP1,j2+2*nshapeV1+nshapeP1) += fact;
            ek(j2+2*nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact;
        }
        
    }

    
    std::cout<<ek<<std::endl;
    
    
    
//    const int pindex = this->PIndex();
//    const int vindex = this->VIndex();
//
//
//    TPZFMatrix<REAL> &dphiLdAxes = datavecleft[fb].dphix;
//    TPZFMatrix<REAL> &dphiRdAxes = datavecright[fb].dphix;
//    TPZFMatrix<REAL> &phiL = datavecleft[fb].phi;
//    TPZFMatrix<REAL> &phiR = datavecright[fb].phi;
//    TPZManVector<REAL,3> &normal = data.normal;
//
//    TPZFNMatrix<660> dphiL, dphiR;
//    TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
//    TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
//    
//    int &LeftPOrder=dataleft.p;
//    int &RightPOrder=dataright.p;
//    
//    REAL &faceSize=data.HSize;
//    
//    
//    int nrowl = phiL.Rows();
//    int nrowr = phiR.Rows();
//    int il,jl,ir,jr,id;
//    
//    //Convection term
//    REAL ConvNormal = 0.;
//    for(id=0; id<fDim; id++) ConvNormal += fC * fConvDir[id] * normal[id];
//    if(ConvNormal > 0.) {
//        for(il=0; il<nrowl; il++) {
//            for(jl=0; jl<nrowl; jl++) {
//                ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
//            }
//        }
//        for(ir=0; ir<nrowr; ir++) {
//            for(jl=0; jl<nrowl; jl++) {
//                ek(ir+nrowl,jl) -= weight * ConvNormal * phiR(ir) * phiL(jl);
//            }
//        }
//    } else {
//        for(ir=0; ir<nrowr; ir++) {
//            for(jr=0; jr<nrowr; jr++) {
//                ek(ir+nrowl,jr+nrowl) -= weight * ConvNormal * phiR(ir) * phiR(jr);
//            }
//        }
//        for(il=0; il<nrowl; il++) {
//            for(jr=0; jr<nrowr; jr++) {
//                ek(il,jr+nrowl) += weight * ConvNormal * phiL(il) * phiR(jr);
//            }
//        }
//    }
//    
//    if(IsZero(fK)) return;
//    //diffusion term
//    STATE leftK, rightK;
//    leftK  = this->fK;
//    rightK = this->fK;
//    
//    // 1) phi_I_left, phi_J_left
//    for(il=0; il<nrowl; il++) {
//        REAL dphiLinormal = 0.;
//        for(id=0; id<fDim; id++) {
//            dphiLinormal += dphiL(id,il)*normal[id];
//        }
//        for(jl=0; jl<nrowl; jl++) {
//            REAL dphiLjnormal = 0.;
//            for(id=0; id<fDim; id++) {
//                dphiLjnormal += dphiL(id,jl)*normal[id];
//            }
//            ek(il,jl) += (STATE)(weight * ( this->fSymmetry * (0.5)*dphiLinormal*phiL(jl,0)-(0.5)*dphiLjnormal*phiL(il,0))) * leftK;
//        }
//    }
//    
//    // 2) phi_I_right, phi_J_right
//    for(ir=0; ir<nrowr; ir++) {
//        REAL dphiRinormal = 0.;
//        for(id=0; id<fDim; id++) {
//            dphiRinormal += dphiR(id,ir)*normal[id];
//        }
//        for(jr=0; jr<nrowr; jr++) {
//            REAL dphiRjnormal = 0.;
//            for(id=0; id<fDim; id++) {
//                dphiRjnormal += dphiR(id,jr)*normal[id];
//            }
//            ek(ir+nrowl,jr+nrowl) += (STATE)(weight * (this->fSymmetry * ((-0.5) * dphiRinormal * phiR(jr) ) + (0.5) * dphiRjnormal * phiR(ir))) * rightK;
//        }
//    }
//    
//    // 3) phi_I_left, phi_J_right
//    for(il=0; il<nrowl; il++) {
//        REAL dphiLinormal = 0.;
//        for(id=0; id<fDim; id++) {
//            dphiLinormal += dphiL(id,il)*normal[id];
//        }
//        for(jr=0; jr<nrowr; jr++) {
//            REAL dphiRjnormal = 0.;
//            for(id=0; id<fDim; id++) {
//                dphiRjnormal += dphiR(id,jr)*normal[id];
//            }
//            ek(il,jr+nrowl) += (STATE)weight * ((STATE)fSymmetry * ((STATE)((-0.5) * dphiLinormal * phiR(jr)) * leftK ) - (STATE)((0.5) * dphiRjnormal * phiL(il))* rightK );
//        }
//    }
//    
//    // 4) phi_I_right, phi_J_left
//    for(ir=0; ir<nrowr; ir++) {
//        REAL dphiRinormal = 0.;
//        for(id=0; id<fDim; id++) {
//            dphiRinormal += dphiR(id,ir)*normal[id];
//        }
//        for(jl=0; jl<nrowl; jl++) {
//            REAL dphiLjnormal = 0.;
//            for(id=0; id<fDim; id++) {
//                dphiLjnormal += dphiL(id,jl)*normal[id];
//            }
//            ek(ir+nrowl,jl) += (STATE)weight * (
//                                                (STATE)(fSymmetry * (0.5) * dphiRinormal * phiL(jl)) * rightK + (STATE)((0.5) * dphiLjnormal * phiR(ir)) * leftK
//                                                );
//        }
//    }
//    
//    if (this->IsSymetric()){
//        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
//    }
//    
//    if (this->fPenaltyConstant == 0.) return;
//    
//    leftK  = this->fK;
//    rightK = this->fK;
//    
//    
//    
//    //penalty = <A p^2>/h
//    REAL penalty = fPenaltyConstant * (0.5 * (abs(leftK)*LeftPOrder*LeftPOrder + abs(rightK)*RightPOrder*RightPOrder)) / faceSize;
//    
//    if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){
//        
//        // 1) left i / left j
//        for(il=0; il<nrowl; il++) {
//            for(jl=0; jl<nrowl; jl++) {
//                ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
//            }
//        }
//        
//        // 2) right i / right j
//        for(ir=0; ir<nrowr; ir++) {
//            for(jr=0; jr<nrowr; jr++) {
//                ek(ir+nrowl,jr+nrowl) += weight * penalty * phiR(ir,0) * phiR(jr,0);
//            }
//        }
//        
//        // 3) left i / right j
//        for(il=0; il<nrowl; il++) {
//            for(jr=0; jr<nrowr; jr++) {
//                ek(il,jr+nrowl) += -1.0 * weight * penalty * phiR(jr,0) * phiL(il,0);
//            }
//        }
//        
//        // 4) right i / left j
//        for(ir=0; ir<nrowr; ir++) {
//            for(jl=0; jl<nrowl; jl++) {
//                ek(ir+nrowl,jl) += -1.0 * weight *  penalty * phiL(jl,0) * phiR(ir,0);
//            }
//        }
//        
//    }
//    
//    if (this->fPenaltyType == EFluxPenalty || this->fPenaltyType == EBoth){
//        
//        REAL NormalFlux_i = 0.;
//        REAL NormalFlux_j = 0.;
//        
//        // 1) left i / left j
//        for(il=0; il<nrowl; il++) {
//            NormalFlux_i = 0.;
//            for(id=0; id<fDim; id++) {
//                NormalFlux_i += dphiL(id,il)*normal[id];
//            }
//            for(jl=0; jl<nrowl; jl++) {
//                NormalFlux_j = 0.;
//                for(id=0; id<fDim; id++) {
//                    NormalFlux_j += dphiL(id,jl)*normal[id];
//                }
//                ek(il,jl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
//            }
//        }
//        
//        // 2) right i / right j
//        for(ir=0; ir<nrowr; ir++) {
//            NormalFlux_i = 0.;
//            for(id=0; id<fDim; id++) {
//                NormalFlux_i += dphiR(id,ir)*normal[id];
//            }
//            for(jr=0; jr<nrowr; jr++) {
//                NormalFlux_j = 0.;
//                for(id=0; id<fDim; id++) {
//                    NormalFlux_j += dphiR(id,jr)*normal[id];
//                }      
//                ek(ir+nrowl,jr+nrowl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
//            }
//        }
//        
//        // 3) left i / right j
//        for(il=0; il<nrowl; il++) {
//            NormalFlux_i = 0.;
//            for(id=0; id<fDim; id++) {
//                NormalFlux_i += dphiL(id,il)*normal[id];
//            }
//            for(jr=0; jr<nrowr; jr++) {
//                NormalFlux_j = 0.;
//                for(id=0; id<fDim; id++) {
//                    NormalFlux_j += dphiR(id,jr)*normal[id];
//                }      
//                ek(il,jr+nrowl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
//            }
//        }
//        
//        // 4) right i / left j
//        for(ir=0; ir<nrowr; ir++) {
//            NormalFlux_i = 0.;
//            for(id=0; id<fDim; id++) {
//                NormalFlux_i += dphiR(id,ir)*normal[id];
//            }
//            for(jl=0; jl<nrowl; jl++) {
//                NormalFlux_j = 0.;
//                for(id=0; id<fDim; id++) {
//                    NormalFlux_j += dphiL(id,jl)*normal[id];
//                }
//                ek(ir+nrowl,jl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
//            }
//        }
//        
//    }
//    
//    

   
    
}


void TPZStokesMaterial::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
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
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> phiVx(fDimension,phiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(phiV, phiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    
    
    
    
    for(int i = 0; i < nshapeV; i++ )
    {
        int iphi = datavec[vindex].fVecShapeIndex[i].second;
        int ivec = datavec[vindex].fVecShapeIndex[i].first;
        TPZFNMatrix<9> GradVnj(fDimension,fDimension),phiVi(fDimension,fDimension);
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                phiVi(e,f)=datavec[vindex].fNormalVec(e,ivec)*phiVx(f,iphi);
            }
        }
        
        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            TPZFNMatrix<9> GradVnj(fDimension,fDimension);
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradVnj(e,f) = datavec[vindex].fNormalVec(e,jvec)*dphiVx(f,jphi)*normal[f];
                    
                }
            }
            
            ek(i,j) += (-1.) * weight * fViscosity * Inner(phiVi, GradVnj) ;
            //ek(i,j) += ftheta*Transpose(ek);
        }
        
    }
    
    
    

    
    //std::cout<<ek<<std::endl;

    
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

