#include "TPZMatHCurlProjection.h"

#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif



TPZMatHCurlProjection::TPZMatHCurlProjection(int id) : TPZVecL2(id)
{

}

/** @brief Default constructor */
TPZMatHCurlProjection::TPZMatHCurlProjection() : TPZVecL2()
{

}


TPZMatHCurlProjection::TPZMatHCurlProjection(const TPZMatHCurlProjection &mat) : TPZVecL2(mat)
{
    
}

TPZMatHCurlProjection::~TPZMatHCurlProjection()
{
    
}
void TPZMatHCurlProjection::RotateForHCurl(TPZVec<REAL> normal , TPZFMatrix<REAL> vHdiv , TPZFMatrix<REAL> &vHcurl ){
    int nFunctions = vHdiv.Rows();
    vHcurl.Resize( vHdiv.Rows(), vHdiv.Cols());
    vHcurl.Zero();
    
    for (int i = 0 ; i < nFunctions; i++) {
        
        vHcurl(i,0) = normal[1]*vHdiv(i,2) - vHdiv(i,1)*normal[2];
        vHcurl(i,1) = normal[2]*vHdiv(i,0) - vHdiv(i,2)*normal[0];
        vHcurl(i,2) = normal[0]*vHdiv(i,1) - vHdiv(i,0)*normal[1];
    }
}
void TPZMatHCurlProjection::ComputeCurl(TPZFMatrix<REAL> gradScalarPhi , TPZFMatrix<REAL> ivecHCurl , TPZFMatrix<REAL> &curlPhi ){
    int nFunctions = gradScalarPhi.Rows();
    curlPhi.Resize( nFunctions , 3);
    curlPhi.Zero();
    TPZManVector<REAL,3> result(3,0.);
    
    for (int i = 0 ; i < nFunctions; i++) {
        curlPhi(i,0) = gradScalarPhi(i,1)*ivecHCurl(i,2) - ivecHCurl(i,1)*gradScalarPhi(i,2);
        curlPhi(i,1) = gradScalarPhi(i,2)*ivecHCurl(i,0) - ivecHCurl(i,2)*gradScalarPhi(i,0);
        curlPhi(i,2) = gradScalarPhi(i,0)*ivecHCurl(i,1) - ivecHCurl(i,0)*gradScalarPhi(i,1);
    }
}



void TPZMatHCurlProjection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    /*********************CREATE HDIV FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiScaHCurl = data.phi;    
    int phrq = data.fVecShapeIndex.NElements();
    //  std::cout<<"x"<<std::endl<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
    
    TPZFNMatrix< 36 , REAL > phiVecHDiv(phrq , 3 , 0.);
    for (int iq = 0 ; iq < phrq ; iq++) {
        int ivecind = data.fVecShapeIndex[iq].first;
        int ishapeind = data.fVecShapeIndex[iq].second;
        
        phiVecHDiv(iq , 0) = phiScaHCurl(ishapeind , 0) * data.fNormalVec(0,ivecind);
        phiVecHDiv(iq , 1) = phiScaHCurl(ishapeind , 0) * data.fNormalVec(1,ivecind);
        phiVecHDiv(iq , 2) = phiScaHCurl(ishapeind , 0) * data.fNormalVec(2,ivecind);
    }
    
    /*********************CALCULATE NORMAL VECTOR****************************/
    TPZManVector<REAL,3> ax1(3),ax2(3), elNormal(3);
    for (int i=0; i<3; i++) {
        ax1[i] = data.axes(0,i);
        ax2[i] = data.axes(1,i);
    }
    Cross(ax1, ax2, elNormal);
    
    /*********************CREATE HCURL FUNCTIONS*****************************/
    TPZFNMatrix< 36 , REAL > phiVecHCurl(phrq , 3 , 0.);
    RotateForHCurl(elNormal , phiVecHDiv , phiVecHCurl);
    /*********************COMPUTE CURL****************************/
    TPZFMatrix<REAL> &dphiQdaxes = data.dphix;
    TPZFNMatrix<3,REAL> dphiQ;
    TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);
    TPZFNMatrix<3,REAL> gradPhiForHCurl(phrq , 3 , 0.);
    TPZFNMatrix<3,REAL> ivecHCurl(phrq , 3 , 0.);
    TPZManVector<REAL,3> iVecHDiv(3,0.), ivecForCurl(3,0.);
    for (int iPhi = 0; iPhi < phrq; iPhi++) {
        int ivecind = data.fVecShapeIndex[iPhi].first;
        int ishapeind = data.fVecShapeIndex[iPhi].second;
        iVecHDiv[0] = data.fNormalVec(0,ivecind);
        iVecHDiv[1] = data.fNormalVec(1,ivecind);
        iVecHDiv[2] = data.fNormalVec(2,ivecind);
        Cross(elNormal, iVecHDiv, ivecForCurl);
        for (int i = 0; i<dphiQ.Rows(); i++) {
            gradPhiForHCurl(iPhi,i) = dphiQ(i,ishapeind);
            ivecHCurl(iPhi,i) = ivecForCurl[i];
        }
    }
    TPZFNMatrix<40,REAL> curlPhi;
    ComputeCurl(gradPhiForHCurl, ivecHCurl, curlPhi);
   	//*****************GET FORCING FUNCTION****************//
   	TPZManVector<STATE,3> force(3);
    force.Fill(0.);
    if(fForcingFunction) {
        fForcingFunction->Execute(data.x,force);
    }
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
    
    const int nHCurlFunctions  = phrq;

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        STATE load = 0.;
        load += phiVecHCurl(iVec , 0 ) * force[0];
        load += phiVecHCurl(iVec , 1 ) * force[1];
        load += phiVecHCurl(iVec , 2 ) * force[2];
        ef( iVec ) += load * weight;
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE stiff = 0.;
            
            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += phiVecHCurl(iVec , 0) * phiVecHCurl(jVec , 0);
            phiIdotPhiJ += phiVecHCurl(iVec , 1) * phiVecHCurl(jVec , 1);
            phiIdotPhiJ += phiVecHCurl(iVec , 2) * phiVecHCurl(jVec , 2);
            
            stiff = phiIdotPhiJ;
            
            
//            if( iVec == jVec ){
//                std::cout<<"stiffBtt "<<iVec<<" "<<jVec<<":"<<stiffBtt<<std::endl;
//            }
            //ek( firstHCurl + iVec , firstHCurl + jVec ) += curlIdotCurlJ * weight ;
            ek( iVec ,jVec ) += stiff * weight ;

        }
    }
}

void TPZMatHCurlProjection::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatHCurlProjection::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatHCurlProjection::ContributeForcingRTBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    DebugStop();
}

void TPZMatHCurlProjection::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
}

void TPZMatHCurlProjection::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();

}

void TPZMatHCurlProjection::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatHCurlProjection::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

int TPZMatHCurlProjection::IntegrationRuleOrder(int elPMaxOrder) const
{
    return 2+elPMaxOrder*2;
}


int TPZMatHCurlProjection::VariableIndex(const std::string &name)
{
    if( strcmp(name.c_str(), "Et") == 0) return 0;
    DebugStop();
    return 1;
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatHCurlProjection::NSolutionVariables(int var)
{
    switch (var) {
        case 0: //Et
            return 2;
            break;
        default:
            DebugStop();
            break;
    }
    return 1;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatHCurlProjection::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    
    TPZVec<STATE> et(3,0.);
    TPZManVector<STATE,3> ax1(3),ax2(3), normal(3);
    for (int i=0; i<3; i++) {
        ax1[i] = data.axes(0,i);
        ax2[i] = data.axes(1,i);
    }
    //ROTATE FOR HCURL
    Cross(ax1, ax2, normal);
    
    Cross(normal,data.sol[0], et);
    
    switch (var) {
        case 0: //Et
        {
            Solout = et;
        }
            break;
        default:
            DebugStop();
            break;
    }
}


/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatHCurlProjection::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    
    DebugStop();
}






