//
//  TPZLagrangeMultiplier.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#include "TPZLagrangeMultiplier.h"
#include "pzaxestools.h"



/** @brief Unique identifier for serialization purposes */
int TPZLagrangeMultiplier::ClassId() const
{
    return TPZLagrangeMultiplierID;
}

/** @brief Saves the element data to a stream */
void TPZLagrangeMultiplier::Write(TPZStream &buf, int withclassid)
{
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
}

/** @brief Reads the element data from a stream */
void TPZLagrangeMultiplier::Read(TPZStream &buf, void *context)
{
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fNStateVariables);
    
}

//Contribution of skeletal elements.
void TPZLagrangeMultiplier::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    int nmesh = datavec.size();
    if (nmesh!=2) DebugStop();

    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    int phrq = phiQ.Rows();
    int phrp = phiP.Rows();
    
//------- Block of matrix B ------
    int iq, jp;
	for(iq = 0; iq<phrq; iq++) {
		for(jp=0; jp<phrp; jp++) {
            ek(iq, phrq+jp) += fMultiplier*weight*phiQ(iq,0)*phiP(jp,0);
		}
	}
    
    
//------- Block of matrix B^T ------
    int ip, jq;
	for(ip=0; ip<phrp; ip++) {
		for(jq=0; jq<phrq; jq++) {
			ek(ip + phrq,jq) += fMultiplier*weight*phiP(ip,0)*phiQ(jq,0);
		}
	}
}

/**
 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZLagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if(dataleft.fVecShapeIndex.size())
    {
        ContributeInterfaceVecLeft(data, dataleft, dataright, weight, ek, ef);
        return;
    }
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
    int secondblock = ek.Rows()-phiR.Rows();
	int il,jl,ir,jr;
    
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl; il++) {
		for(jr=0; jr<nrowr; jr++) {
			ek(il,jr+secondblock) += weight * fMultiplier * (phiL(il) * phiR(jr));
		}
	}
	
    //	// 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr; ir++) {
		for(jl=0; jl<nrowl; jl++) {
			ek(ir+secondblock,jl) += weight * fMultiplier * (phiR(ir) * phiL(jl));
		}
	}
    
}

void TPZLagrangeMultiplier::ContributeInterfaceVecLeft(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> &phiL = dataleft.phi;
    TPZFMatrix<REAL> &phiR = dataright.phi;
    
    int nrowl = dataleft.fVecShapeIndex.size();
    int nrowr = phiR.Rows()*fNStateVariables;
    
#ifdef DEBUG
    if (ek.Rows() != nrowl+nrowr) {
        DebugStop();
    }
#endif
    int secondblock = nrowl;
    int il,jl,ir,jr;
    // 3) phi_I_left, phi_J_right
    for(il=0; il<nrowl; il++) {
        int ivec = dataleft.fVecShapeIndex[il].first;
        int iphindex = dataleft.fVecShapeIndex[il].second;
        TPZManVector<REAL,3> vec(fDimension,0);
        for(int i=0; i< fDimension; i++) vec[i] = dataleft.fNormalVec(i,ivec);
        for(jr=0; jr<nrowr; jr++) {
            REAL val = weight * fMultiplier * (phiL(iphindex) *vec[jr%fDimension] * phiR(jr/fDimension));
            ek(il,jr+secondblock) += val;
            ek(jr+secondblock,il) += val;
        }
    }
}

/**
 * @brief It computes a contribution to residual vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZLagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}


