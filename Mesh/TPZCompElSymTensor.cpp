/**
 * @file
 * @brief Contains the implementation of the TPZCompElSymTensor methods.
 */

#include "pzcmesh.h"
#include "TPZCompElSymTensor.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"
#include "pzgeoquad.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzhdivpressure.h"
#include "pzshapepiram.h"



#include "pzshtmat.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElSymTensor"));
#endif

using namespace std;


template<class TSHAPE>
TPZCompElSymTensor<TSHAPE>::TPZCompElSymTensor(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) :
TPZIntelGen<TSHAPE>(mesh,gel,index,1), fSideOrient(TSHAPE::NFaces,1) {
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	gel->SetReference(this);
	
    for(int i=0;i< TSHAPE::NSides; i++)
	{
        int sideaux = i;
		this->fConnectIndexes[i] = this->CreateMidSideConnect(sideaux);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) 
        {
            std::stringstream sout;
            sout << "After creating last flux connect " << i << std::endl;
            //	this->Print(sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        if (i< TSHAPE::NCornerNodes)
        {
            TPZConnect &c = mesh.ConnectVec()[this->fConnectIndexes[i]];
            c.SetNState(3);
            long seqnum = c.SequenceNumber();
            mesh.Block().Set(seqnum, 3);
        }
        else
        {
        }
		mesh.ConnectVec()[this->fConnectIndexes[i]].IncrementElConnected();
    }
    
	
    /// initializing the integration rule
    int sideorder = SideOrder(TSHAPE::NSides-1);
    sideorder++;
	sideorder = 2*sideorder;
	if (sideorder > this->fIntRule.GetMaxOrder()) sideorder = this->fIntRule.GetMaxOrder();
	TPZManVector<int,3> order(3,sideorder);
	this->fIntRule.SetOrder(order);
    int firstside = TSHAPE::NSides-TSHAPE::NFaces-1;
    for(int side = firstside ; side < TSHAPE::NSides-1; side++ )
    {
        fSideOrient[side-firstside] = this->Reference()->NormalOrientation(side);
    }
    TPZMaterial *mat = this->Material();
    if (mat)
    {
        int order = mat->IntegrationRuleOrder(MaxOrder());
        TPZManVector<int,3> ord(gel->Dimension(),order);
        this->fIntRule.SetOrder(ord);
    }

    
}

template<class TSHAPE>
TPZCompElSymTensor<TSHAPE>::TPZCompElSymTensor(TPZCompMesh &mesh, const TPZCompElSymTensor<TSHAPE> &copy) :
TPZIntelGen<TSHAPE>(mesh,copy), fSideOrient(copy.fSideOrient)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		this-> fConnectIndexes[i] = copy.fConnectIndexes[i];
	}
    
}

template<class TSHAPE>
TPZCompElSymTensor<TSHAPE>::TPZCompElSymTensor(TPZCompMesh &mesh,
									 const TPZCompElSymTensor<TSHAPE> &copy,
									 std::map<long,long> & gl2lcConMap,
									 std::map<long,long> & gl2lcElMap) :
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fSideOrient(copy.fSideOrient)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		int lcIdx = -1;
		int glIdx = copy.fConnectIndexes[i];
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			this-> fConnectIndexes[i] = -1;
			return;
		}
		this-> fConnectIndexes[i] = lcIdx;
	}
}

template<class TSHAPE>
TPZCompElSymTensor<TSHAPE>::TPZCompElSymTensor() :
TPZIntelGen<TSHAPE>()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
		this-> fConnectIndexes[i] = -1;
	}
    
}

template<class TSHAPE>
TPZCompElSymTensor<TSHAPE>::~TPZCompElSymTensor(){
    TPZGeoEl *gel = this->Reference();
    if (gel->Reference() != this) {
        DebugStop();
    }
    gel->ResetReference();
}

template<class TSHAPE>
MElementType TPZCompElSymTensor<TSHAPE>::Type() {
	return TSHAPE::Type();
}



template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::SetConnectIndex(int i, long connectindex){
#ifndef NODEBUG
	if(i<0 || i>= this->NConnects()) {
		std::cout << " TPZCompElSymTensor<TSHAPE>::SetConnectIndex index " << i <<
		" out of range\n";
		DebugStop();
		return;
	}
#endif
	this-> fConnectIndexes[i] = connectindex;
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template<class TSHAPE>
int TPZCompElSymTensor<TSHAPE>::NConnectShapeF(int connect)const
{
    int ncon = this->NConnects();
    if (connect < TSHAPE::NCornerNodes) {
        return 1;
    }
    else if (connect < this->NConnects())
    {
        int order = ConnectOrder(connect);
        return TSHAPE::NConnectShapeF(connect,order);
    }

    DebugStop();
    return -1;
 }
 



////
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::SetIntegrationRule(int ord) {
	TPZManVector<int,3> order(TSHAPE::Dimension,ord);
	this->fIntRule.SetOrder(order);
}




template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
    ord.Resize(TPZIntelGen<TSHAPE>::NConnects());
	int i;
    for(i=0; i<TPZIntelGen<TSHAPE>::NConnects(); i++) {
		ord[i] = ConnectOrder(i);
	}
}


template<class TSHAPE>
int TPZCompElSymTensor<TSHAPE>::PreferredSideOrder(int side) {
	if(TSHAPE::SideDimension(side) < Dimension()-1)
	{
		PZError << __PRETTY_FUNCTION__ << " side " << side << std::endl;
	}
	int connect= side;
    if(connect<0 || connect > TPZIntelGen<TSHAPE>::NConnects()) {
		PZError << "TPZCompElSymTensor<TSHAPE>::PreferredSideOrder no polynomial associate " <<  endl;
		return -1;
	}
    if(connect<TPZIntelGen<TSHAPE>::NConnects()) {
			int order =this->fPreferredOrder;
			return order;//this->AdjustPreferredSideOrder(side,order);
	}
	PZError << "TPZCompElSymTensor<TSHAPE>::PreferredSideOrder called for connect = " << connect << "\n";
	return 0;
	
}

template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::SetPreferredOrder(int order)
{
    TPZIntelGen<TSHAPE>:: SetPreferredOrder(order);
}

template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::SetSideOrder(int side, int order){
	int connectaux= side;
	if(connectaux<0 || connectaux > this-> NConnects()) {
		PZError << "TPZCompElSymTensor::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef LOG4CXX
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
	
    }
    if (side < TSHAPE::NCornerNodes && order != 0) {
        DebugStop();
    }
	TPZConnect &c = this->Connect(connectaux);
    long cindex = this->ConnectIndex(connectaux);
    c.SetOrder(order,cindex);
    long seqnum = c.SequenceNumber();
    int nvar = 1;
    int nshape =this-> NConnectShapeF(connectaux);
    if (side < TSHAPE::NCornerNodes) {
        nshape = 1;
        nvar = 3;
    }
    if (side < TSHAPE::NSides-1) {
        nshape = TSHAPE::NConnectShapeF(side,order);
        nvar = 2;
    }
    if (side == TSHAPE::NSides - 1) {
        nvar = 3;
        int nshape_interior = TSHAPE::NConnectShapeF(side,order);
        int nside_shape = TSHAPE::NConnectShapeF(side-1,order);
        int totalshape = nshape_interior*3+3*nside_shape;
        nshape = nshape_interior+nside_shape;
    }
    c.SetNState(nvar);
    c.SetNShape(nshape);
	this-> Mesh()->Block().Set(seqnum,nshape*nvar);
}


template<class TSHAPE>
int TPZCompElSymTensor<TSHAPE>::ConnectOrder(int connect) const{
	if (connect < 0 || connect >= this->NConnects()){
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Connect index out of range connect " << connect <<
			" nconnects " << this->NConnects();
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return -1;
	}
	
	if (this->fConnectIndexes[connect] == -1) {
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " connect " << connect
		<< " is not initialized" << std::endl;
#ifdef LOG4CXX
		LOGPZ_ERROR(logger,sout.str());
#else
		std::cout << sout.str() << std::endl;
#endif
		return 0;
	}
	
    TPZConnect &c = this-> Connect(connect);
    return c.Order();
}

template<class TSHAPE>
int TPZCompElSymTensor<TSHAPE>::SideOrder(int side) const
{
    if(!TPZIntelGen<TSHAPE>::NSideConnects(side)) return -1;
	int corder = side;
	int maxorder = 0;
	int conectaux;
    if(corder>=0 || corder <= TPZIntelGen<TSHAPE>::NConnects()) return ConnectOrder(corder);
    
    DebugStop();
	TPZStack< int > high;
	TSHAPE::HigherDimensionSides(side, high);
	int highside= high.NElements();
	
	
	for(int j=0;j<highside;j++)
	{
		conectaux = high[j];
		maxorder = (ConnectOrder(conectaux) > maxorder) ? ConnectOrder(conectaux) : maxorder;
	}
	
	return maxorder;	
}

/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::FirstShapeIndex(TPZVec<long> &Index) const {
	
    int highorder = SideOrder(TSHAPE::NSides-1);
    TPZManVector<int> orders(TSHAPE::NSides-TSHAPE::NCornerNodes,highorder);
    Index[0] = 0;
	for(int iside=0;iside<TSHAPE::NSides;iside++)
	{
        int sideorder = 1;
        if (iside >= TSHAPE::NCornerNodes) {
            sideorder = orders[iside-TSHAPE::NCornerNodes];
        }
        int temp = Index[iside] + TSHAPE::NConnectShapeF(iside,sideorder);
        Index[iside+1] = temp;
	}
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "First  Index " << Index;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
}


template<class TSHAPE>
int TPZCompElSymTensor<TSHAPE>::NFluxShapeF() const{
    int in,result=0;
    int nn=TPZCompElSymTensor::NConnects();
    for(in=0;in<nn;in++){
//#ifdef LOG4CXX
//				std::stringstream sout;
//				sout << "conect " << in<< " seq number "<<seqnum<<" num func "<<TPZCompElSymTensor::NConnectShapeF(in);
//				LOGPZ_DEBUG(logger,sout.str())
//#endif
        result += TPZCompElSymTensor::NConnectShapeF(in);
    }
		
		
    DebugStop();
//#ifdef LOG4CXX
//    std::stringstream sout;
//    sout << "Num funcoes associada ao fluxo " << result;
//    LOGPZ_DEBUG(logger,sout.str())
//#endif
    return result;
		

}

/**
 * @brief Returns a matrix index of the shape and vector  associate to element
 * @param[in] VectorSide Indicates the side associated with each vector
 * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
 * @param[in] pressureorder Order of the pressure (to select shape functions?)
 */
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide, TPZVec<int> &bilinear, TPZVec<int> &direction, TPZVec<std::pair<int,long> > & IndexVecShape, int pressureorder)
{
    
        DebugStop();
        
}


template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,long> > & ShapeAndVec, int pressureorder)
{
    DebugStop();
}


/**
 * @brief It returns the normal orientation of the reference element by the side.
 * Only side that has dimension larger than zero and smaller than me.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
int TPZCompElSymTensor<TSHAPE>::GetSideOrient(int side){
    
    int firstside = TSHAPE::NSides-TSHAPE::NFaces-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    return fSideOrient[side-firstside];
}

/**
 * @brief It set the normal orientation of the element by the side.
 * Only side that has dimension equal to my dimension minus one.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::SetSideOrient(int side, int sideorient){
    
    int firstside = TSHAPE::NSides-TSHAPE::NFaces-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    fSideOrient[side-firstside] = sideorient;
}

//compute the values of the shape function of the side
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {

    
    DebugStop();

}

template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
{
    if (var == 99) {
        return TPZIntelGen<TSHAPE>::Solution(qsi,var,sol);
    }
    TPZMaterialData data;
	InitMaterialData(data);
	//this->ComputeSolutionHDiv(data);
    this->ComputeRequiredData(data,qsi);
    this->ComputeSolution(qsi,data);
	this->Material()->Solution(data,var,sol);
}



template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data){
    
    TPZBlock<STATE> &bl = this->Mesh()->Block();
    TPZManVector<STATE> multipliers(data.fVecShapeIndex.size());
    int count = 0;
    for (int ic = 0; ic < TSHAPE::NSides; ic++) {
        TPZConnect &c = this->Connect(ic);
        int ndof = c.NDof();
        long seqnum = c.SequenceNumber();
        for (int i=0; i<ndof; i++) {
            multipliers[count++] = bl(seqnum,0,i,0);
        }
    }
    data.sol.resize(1);
    data.sol[0].resize(4);
    data.dsol.resize(1);
    data.dsol[0].Resize(2, 4);
    for (int i=0; i<4; i++) {
        data.sol[0][i]=0.;
        data.dsol[0](0,i)=0.;
        data.dsol[0](1,i)=0.;
    }
    for (int i=0; i<data.fVecShapeIndex.size(); i++) {
        int ishape = data.fVecShapeIndex[i].second;
        int ivec = data.fVecShapeIndex[i].first;
        for (int j=0; j<4; j++) {
            data.sol[0][j] += data.fNormalVec(j,ivec)*data.phi(ishape)*multipliers[i];
            for (int d=0; d<2; d++) {
                data.dsol[0](d,j) += data.fNormalVec(j,ivec)*data.dphix(d,ishape)*multipliers[i];
            }
        }
    }
    
}



template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
	TPZManVector<long,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);

    int connectorder = ConnectOrder(TSHAPE::NSides-1);
    TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,connectorder);
    int i;
    TPZGeoEl *ref = this->Reference();
    for(i=0; i<TSHAPE::NCornerNodes; i++) {
        id[i] = ref->NodePtr(i)->Id();
    }


    int nshape= TSHAPE::NShapeF(ord);
    
    phi.Resize(nshape, 1);
    dphi.Resize(TSHAPE::Dimension, nshape);
    TSHAPE::Shape(pt,id,ord,phi,dphi);

}

/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::FirstShapeIndex(TPZVec<int> &Index) const {
    
    TPZConnect &c = this->Connect(TSHAPE::NSides-1);
    TPZManVector<int> orders(TSHAPE::NSides-TSHAPE::NCornerNodes,c.Order());
    Index[0] = 0;
    for(int iside=0;iside<TSHAPE::NSides;iside++)
    {
        int sideorder = 1;
        if (iside >= TSHAPE::NCornerNodes) {
            sideorder = orders[iside-TSHAPE::NCornerNodes];
        }
        int temp = Index[iside] + TSHAPE::NConnectShapeF(iside,sideorder);
        Index[iside+1] = temp;
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "First  Index " << Index;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
}



template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
    
//    TPZManVector<int,TSHAPE::NSides*TSHAPE::Dimension> normalsidesDG(TSHAPE::Dimension*TSHAPE::NSides);



    this->ComputeShape(qsi,data.x,data.jacobian,data.axes,data.detjac,data.jacinv,data.phi,data.dphi,data.dphix);

    if (data.fNeedsSol) {
        this->ComputeSolution(qsi,data);
    }
}//void

/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
//	if (TSHAPE::Type()==EQuadrilateral) {
//        int maxorder = this->MaxOrder();
//        data.p = maxorder+1;
//    }
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
		{
				LOGPZ_DEBUG(logger,"Initializing MaterialData of TPZCompElSymTensor")
		}
#endif
    
    // FOR A TRIANGLE!!!
    // create the three canonic tensors
    TPZFNMatrix<12,REAL> canonical(4,3,0.);
    canonical(0,0) = 1.; // sigxx
    canonical(3,1) = 1.; // sigyy
    canonical(1,2) = 1.; // tauxy
    canonical(2,2) = 1.; // tauxy
                         // get the corner coordinates
    TPZFNMatrix<6,REAL> cornerco(2,3);
    TPZGeoEl *gel = this->Reference();
    for (int in=0; in < TSHAPE::NCornerNodes; in++) {
        TPZManVector<REAL,3> co(3);
        gel->Node(in).GetCoordinates(co);
#ifdef PZDEBUG
        if (co[2] != 0.) {
            DebugStop();
        }
#endif
        for (int c=0; c<2; c++) {
            cornerco(c,in) = co[c];
        }
    }
    // compute the three tangent vectors
    TPZFNMatrix<6,REAL> tangent(2,3);
    for (int t=0; t<3; t++) {
        TPZManVector<REAL,2> diff(2);
        for (int c=0; c<2; c++) {
            diff[c] = cornerco(c,(t+1)%3)-cornerco(c,t);
        }
        REAL norm = sqrt(diff[0]*diff[0]+diff[1]*diff[1]);
        tangent(0,t) = diff[0]/norm;
        tangent(1,t) = diff[1]/norm;
    }
    // rotate to obtain the normal vectors
    TPZFNMatrix<6,REAL> normalvec(2,3);
    for (int n=0; n<3; n++) {
        normalvec(0,n) = tangent(1,n);
        normalvec(1,n) = -tangent(0,n);
    }
    // compute three tensors associated with each side : txt (txn+nxt) and nxn
    TPZFNMatrix<24,REAL> externalsidetensors(4,6);
    for (int side=0; side<3; side++) {
        TPZFNMatrix<4> nn(2,2),tnt(2,2);
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                nn(i,j) = normalvec(i,side)*normalvec(j,side);
                tnt(i,j) = tangent(i,side)*normalvec(j,side)+tangent(j,side)*normalvec(i,side);
            }
        }
        for(int i=0; i<4; i++)
        {
            externalsidetensors(i,2*side) = nn(i%2,i/2);
            externalsidetensors(i,2*side+1) = tnt(i%2,i/2);
        }
    }
    
    TPZFNMatrix<12,REAL> tangentsidetensors(4,3);
    for (int side=0; side<3; side++) {
        TPZFNMatrix<4> tt(2,2);
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                tt(i,j) = tangent(i,side)*tangent(j,side);
            }
        }
        for(int i=0; i<4; i++)
        {
            tangentsidetensors(i,side) = tt(i%2,i/2);
        }
    }
    // fill the datastructure data.fNormalVec
    data.fNormalVec.Redim(4, 12);
    for (int i=0; i<4; i++) {
        for (int j=0; j<3; j++) {
            data.fNormalVec(i,j) = canonical(i,j);
        }
        for (int j=0; j<6; j++) {
            data.fNormalVec(i,j+3) = externalsidetensors(i,j);
        }
        for (int j=0; j<3; j++) {
            data.fNormalVec(i,j+9) = tangentsidetensors(i,j);
        }
    }
    // associate the shape functions with the tensors
    
    // call firstshapeindex to get a shape function counter.
    TPZManVector<long,9> firstshape(TSHAPE::NSides+1);
    this->FirstShapeIndex(firstshape);
    
    int numfunc = 0;
    for (int ic=0; ic<TSHAPE::NSides; ic++) {
        TPZConnect &c = this->Connect(ic);
        numfunc += c.NDof();
    }
    data.fVecShapeIndex.resize(numfunc);
    // compute the total number of tensor functions 3 x cornernodes + 2 x nlinersides + internal (the internal order may be different!)
    // first the corner shape functions
    int count = 0;
    for (int c=0; c<TSHAPE::NCornerNodes; c++) {
        for (int tensor=0; tensor<3; tensor++) {
            data.fVecShapeIndex[count++] = std::make_pair(tensor, c);
        }
    }
    // two tensors for each side function - use the order of the side connect
    for (int side=0; side<3; side++) {
        TPZConnect &c = this->Connect(side+3);
        int order = c.Order();
        int firstshapeindex = firstshape[3+side];
        for (int shape =0; shape < order-1; shape++) {
            for (int tensor = 0; tensor<2; tensor++) {
                data.fVecShapeIndex[count++] = std::make_pair(tensor+2*side+3, firstshapeindex+shape);
            }
        }
    }
    // internal functions
    // one tensor for each side function - use the internal order
    for (int side=0; side<3; side++) {
        TPZConnect &c = this->Connect(TSHAPE::NSides-1);
        int order = c.Order();
        int firstshapeindex = firstshape[3+side];
        for (int shape =0; shape < order-1; shape++) {
                data.fVecShapeIndex[count++] = std::make_pair(9+side, firstshapeindex+shape);
        }
    }
    // three canonical tensors for each internal function
    for (int sh=firstshape[TSHAPE::NSides-1]; sh < firstshape[TSHAPE::NSides]; sh++) {
        for (int tensor = 0; tensor<3; tensor++) {
            data.fVecShapeIndex[count++] = std::make_pair(tensor, sh);
        }
    }

#ifdef PZDEBUG
    if (count != numfunc) {
        DebugStop();
    }
#endif
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		data.fNormalVec.Print("Normal vector ", sout,EMathematicaInput);
        for (int i=0; i<TSHAPE::NCornerNodes; i++) {
            sout << "Id[" << i << "] = " << this->Reference()->NodePtr(i)->Id() << " ";
        }
        sout << std::endl;
		sout << "NormalVector/Shape indexes \n";
        for (int i=0; i<data.fVecShapeIndex.size(); i++) {
            sout << i << '|' << data.fVecShapeIndex[i] << " ";
        }
        sout << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif    
    
}


// Save the element data to a stream
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
	TPZInterpolatedElement::Write(buf,withclassid);
	TPZManVector<int,3> order(3,0);
	this->fIntRule.GetOrder(order);
	this->WriteObjects(buf,order);
	buf.Write(this->fConnectIndexes.begin(),TSHAPE::NSides);
	buf.Write(&this->fPreferredOrder,1);
    this->WriteObjects(buf,fSideOrient);
	int classid = this->ClassId();
	buf.Write ( &classid, 1 );
}


// Read the element data from a stream
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZInterpolatedElement::Read(buf,context);
	TPZManVector<int,3> order;
	this-> ReadObjects(buf,order);
	this-> fIntRule.SetOrder(order);
    TPZManVector<int, TSHAPE::NFaces> SideOrient;
    this-> ReadObjects(buf,SideOrient);
    fSideOrient = SideOrient;
	buf.Read(this->fConnectIndexes.begin(),TSHAPE::NSides);
	buf.Read(&this->fPreferredOrder,1);
    this->ReadObjects(buf,fSideOrient);
	int classid = -1;
	buf.Read( &classid, 1 );
	if ( classid != this->ClassId() )
	{
		std::stringstream sout;
		sout << "ERROR - " << __PRETTY_FUNCTION__
        << " trying to restore an object id " << this->ClassId() << " and classid read = " << classid;
		LOGPZ_ERROR ( logger, sout.str().c_str() );
	}
}
//refinamento
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::PRefine(int order)
{
    this->SetPreferredOrder(order);
    int side;
    int icon;
    int ncon=TPZIntelGen<TSHAPE>::NConnects();
    int nnodes = this->Reference()->NNodes();
    for(icon=0; icon<ncon; icon++)
    {//somente para os conects de fluxo
//        TPZConnect &con = this->Connect(icon);
//        con.SetOrder(order);
        side= icon;
        
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
                std::stringstream sout;
                sout << "side " << side << " order " << this->PreferredSideOrder(side)<<std::endl;
                LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        this->IdentifySideOrder(side);
    }
		// conect da pressao
    
		
}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    out << "Side orientation " << fSideOrient << std::endl;
    
    TPZIntelGen<TSHAPE>::Print(out);

    
}

template<class TSHAPE>
int TPZCompElSymTensor<TSHAPE>::MaxOrder(){
    
    int maxorder = TPZInterpolationSpace::MaxOrder();
    return maxorder;
}

#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt2dmapped.h"

using namespace pztopology;

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"

#include "pzmeshid.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;

//template<>
//void TPZCompElSymTensor<TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
//	if(dimension == 0) std::cout << "A point element has no graphical representation\n";
//}

template<class TSHAPE>
void TPZCompElSymTensor<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}

//template<>
//int TPZCompElSymTensor<TPZShapePoint>::ClassId() const
//{
//	return TPZHDIVPOINTID;
//}

template<>
int TPZCompElSymTensor<TPZShapeTriang>::ClassId() const
{
	return TPZSYMTENSORTRIANGLEID;
}

//template class
//TPZRestoreClass< TPZCompElSymTensor<TPZShapePoint>, TPZHDIVPOINTID>;

template class
TPZRestoreClass< TPZCompElSymTensor<TPZShapeTriang>, TPZSYMTENSORTRIANGLEID>;



template class TPZCompElSymTensor<TPZShapeTriang>;


#include "TPZCompElSymTensorBound.h"

TPZCompEl * CreateSyMTensorBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElSymTensorBound< TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateSymTensorTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElSymTensor< TPZShapeTriang >(mesh,gel,index);
}

