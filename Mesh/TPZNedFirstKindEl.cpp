#include "TPZNedFirstKindEl.h"

void TPZNedFirstKindEl::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                         TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                         REAL &detjac, TPZFMatrix<REAL> &jacinv,
                                         TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx){
	TPZGeoEl * ref = this->Reference();
	if (!ref){
		DebugStop();
	} else if ( !(Reference()->Type() == ETriangle && Reference()->Dimension() == 2) ){
	//only triangular elements for now
		DebugStop();
	}

	ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
	this->Shape(intpoint,phi,dphi);
    this->Convert2Axes(dphi, jacinv, dphidx);   
}

void TPZNedFirstKindEl::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
    DebugStop();
    
    TPZGeoEl *ref = this->Reference();    
}

void TPZNedFirstKindEl::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
	DebugStop();
	data.intGlobPtIndex = -1;
    this->ComputeShape(qsi, data);
    
    if (data.fNeedsSol){
        if (data.phi.Rows()){//if shape functions are available
            this->ComputeSolution(qsi, data);
        }
        else{//if shape functions are not available
            this->ComputeSolution(qsi, data.sol, data.dsol, data.axes);
        }
    }//fNeedsSol
	
    data.x.Resize(3., 0.0);
    Reference()->X(qsi, data.x);

    TPZManVector<REAL,3> x_center(3,0.0);
    TPZVec<REAL> center_qsi(3,0.0);
    
    if (Reference()->Type() == EQuadrilateral && Reference()->Dimension() == 2)
    {
        center_qsi[0] = 0.0;
        center_qsi[1] = 0.0;
    }
    
    if (Reference()->Type() == ETriangle && Reference()->Dimension() == 2)
    {
        center_qsi[0] = 0.25;
        center_qsi[1] = 0.25;
    }

    Reference()->X(center_qsi, x_center);
    data.XCenter = x_center;
    
    if (data.fNeedsHSize){
        data.HSize = 2.*this->InnerRadius();
    }//fNeedHSize
    
    if (data.fNeedsNormal){
        this->ComputeNormal(data);
    }//fNeedsNormal
    

}//void TPZNedFirstKindEl<TSHAPE>::ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi)

template<class TSHAPE>
void TPZNedFirstKindEl<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data){
    DebugStop();
}//void TPZNedFirstKindEl<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data)

template class TPZCompElHDiv<TPZShapeTriang>;
template class TPZCompElHDiv<TPZShapeLinear>;
template class TPZCompElHDiv<TPZShapePoint>;

TPZCompEl * CreateNedFirstKindTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    return new TPZNedFirstKindEl< TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateNedFirstKindLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    return new TPZCompElHDiv< TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateNedFirstKindPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    return new TPZCompElHDiv< TPZShapePoint>(mesh,gel,index);
}