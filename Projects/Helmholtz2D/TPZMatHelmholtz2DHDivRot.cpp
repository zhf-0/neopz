              
                           
            
           
        
      
     
    
 
#include "TPZMatHelmholtz2DHDivRot.h"
#include "TPZVecL2.h"

#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif

TPZMatHelmholtz2DHDivRot::TPZMatHelmholtz2DHDivRot(int id,
                                     STATE (&cFunc)(const TPZVec<REAL> &))
    : TPZMatHelmholtz2D(id,cFunc){

}

TPZMatHelmholtz2DHDivRot::TPZMatHelmholtz2DHDivRot(int id) : TPZMatHelmholtz2D(id) {

}

/** @brief Default constructor */
TPZMatHelmholtz2DHDivRot::TPZMatHelmholtz2DHDivRot() : TPZMatHelmholtz2D(){
}

TPZMatHelmholtz2DHDivRot::TPZMatHelmholtz2DHDivRot(const TPZMatHelmholtz2DHDivRot &mat)
    : TPZMatHelmholtz2D(mat) {

}

TPZMatHelmholtz2DHDivRot::~TPZMatHelmholtz2DHDivRot() {}

void TPZMatHelmholtz2DHDivRot::RotateForHCurl(TPZVec<REAL> normal , TPZFMatrix<REAL> vHdiv , TPZFMatrix<REAL> &vHcurl ){
    int nFunctions = vHdiv.Rows();
    vHcurl.Resize( vHdiv.Rows(), vHdiv.Cols());
    vHcurl.Zero();
    
    for (int i = 0 ; i < nFunctions; i++) {
        
        vHcurl(i,0) = normal[1]*vHdiv(i,2) - vHdiv(i,1)*normal[2];
        vHcurl(i,1) = normal[2]*vHdiv(i,0) - vHdiv(i,2)*normal[0];
        vHcurl(i,2) = normal[0]*vHdiv(i,1) - vHdiv(i,0)*normal[1];
    }
}
void TPZMatHelmholtz2DHDivRot::ComputeCurl(TPZFMatrix<REAL> gradScalarPhi , TPZFMatrix<REAL> ivecHCurl , TPZFMatrix<REAL> &curlPhi ){
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

void TPZMatHelmholtz2DHDivRot::Contribute(TPZMaterialData &data, REAL weight,
                                   TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) {

    /*********************CREATE HDIV FUNCTIONS****************************/
    TPZFNMatrix<12,REAL> phiScaHCurl = data.phi;    
    int phrq = data.fVecShapeIndex.NElements();
    //  std::cout<<"x"<<std::endl<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
    
    TPZFNMatrix< 36 , REAL > phiHDiv(phrq , 3 , 0.);
    for (int iq = 0 ; iq < phrq ; iq++) {
        int ivecind = data.fVecShapeIndex[iq].first;
        int ishapeind = data.fVecShapeIndex[iq].second;
        
        phiHDiv(iq , 0) = phiScaHCurl(ishapeind , 0) * data.fNormalVec(0,ivecind);
        phiHDiv(iq , 1) = phiScaHCurl(ishapeind , 0) * data.fNormalVec(1,ivecind);
        phiHDiv(iq , 2) = phiScaHCurl(ishapeind , 0) * data.fNormalVec(2,ivecind);
    }
    
    /*********************CALCULATE NORMAL VECTOR****************************/
    TPZManVector<REAL,3> ax1(3),ax2(3), elNormal(3);
    for (int i=0; i<3; i++) {
        ax1[i] = data.axes(0,i);
        ax2[i] = data.axes(1,i);
    }
    Cross(ax1, ax2, elNormal);
    
    /*********************CREATE HCURL FUNCTIONS*****************************/
    TPZFNMatrix< 36 , REAL > phiHCurl(phrq , 3 , 0.);
    RotateForHCurl(elNormal , phiHDiv , phiHCurl);
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

    TPZManVector<REAL, 3> x = data.x;
    const STATE cVal = fC(x);

    //*****************GET FORCING FUNCTION****************//
    TPZManVector<STATE, 3> force(3);
    force.Fill(0.);
    if (fForcingFunction) {
        fForcingFunction->Execute(data.x, force);
    } else {
        DebugStop(); // RHS not set!
    }
	
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//

    const int nHCurlFunctions = phiHCurl.Rows();

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        STATE load = 0.;
        load += phiHCurl(iVec, 0) * force[0];
        load += phiHCurl(iVec, 1) * force[1];
        load += phiHCurl(iVec, 2) * force[2];
        ef(iVec) += load * weight;
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE stiff = 0.;

            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += phiHCurl(iVec, 0) * phiHCurl(jVec, 0);
            phiIdotPhiJ += phiHCurl(iVec, 1) * phiHCurl(jVec, 1);
            phiIdotPhiJ += phiHCurl(iVec, 2) * phiHCurl(jVec, 2);

            STATE curlPhiIvecCurlPhiJ = 0.;
            curlPhiIvecCurlPhiJ += curlPhi(iVec , 0) * curlPhi(jVec , 0);
            curlPhiIvecCurlPhiJ += curlPhi(iVec , 1) * curlPhi(jVec , 1);
            curlPhiIvecCurlPhiJ += curlPhi(iVec , 2) * curlPhi(jVec , 2);
//            curlPhiIvecCurlPhiJ += curlPhi(0, iVec) * curlPhi(0, jVec);
//            curlPhiIvecCurlPhiJ += curlPhi(1, iVec) * curlPhi(1, jVec);
//            curlPhiIvecCurlPhiJ += curlPhi(2, iVec) * curlPhi(2, jVec);

            stiff = curlPhiIvecCurlPhiJ + cVal * phiIdotPhiJ;
            ek(iVec, jVec) += stiff * weight;
        }
    }
}

void TPZMatHelmholtz2DHDivRot::ContributeBC(TPZMaterialData &data, REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
	TPZFMatrix<REAL> &phiHCurl = data.phi;

    int nHCurlFunctions = phiHCurl.Rows();
    REAL BIG = TPZMaterial::gBigNumber;

    // const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de
    // condicao mista
    const STATE v2 = bc.Val2()(0, 0); // sera posto no vetor F

    switch (bc.Type()) {
    case 0:
        for (int i = 0; i < nHCurlFunctions; i++) {
            const STATE rhs = phiHCurl(i, 0) * BIG * v2;
            ef(i, 0) += rhs * weight;
            for (int j = 0; j < nHCurlFunctions; j++) {
                const STATE stiff = phiHCurl(i, 0) * phiHCurl(j, 0) * BIG;
                ek(i, j) += stiff * weight;
            }
        }
        break;
    case 1:
        DebugStop();
        break;
    case 2:
        DebugStop();
        break;
    }
}

/** @brief Returns the solution associated with the var index based on the
 * finite element approximation */
void TPZMatHelmholtz2DHDivRot::Solution(TPZMaterialData &data, int var,
                                 TPZVec<STATE> &Solout) {
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
    case 0: // E
    {
        Solout = et;
    } break;
	case 1: // curlE
	{
        Solout[0] = data.dsol[0](0,0)+data.dsol[0](1,1);
    } break;
    default:
        DebugStop();
        break;
    }
}
