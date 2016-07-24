/**
 * @file
 * @brief Contains the implementation of the TPZCompElSymTensorBound methods.
 */


#include "TPZCompElSymTensorBound.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzelchdiv.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElSymTensorBound"));
#endif

template<class TSHAPE>
TPZCompElSymTensorBound<TSHAPE>::TPZCompElSymTensorBound(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) :
TPZIntelGen<TSHAPE>(mesh,gel,index), fSideOrient(1){
		
	//int i;
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	//for(i=0; i<TSHAPE::NSides; i++) this->fConnectIndexes[i]=-1;
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
		{
				std::stringstream sout;
				sout << "After creating boundary flux connect " << this->fConnectIndexes[0] << std::endl;
				//	this->Print(sout);
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
	//TPZGeoElSide myInnerSide(gel,gel->NSides()-1);
//	TPZGeoElSide neigh = myInnerSide.Neighbour();
//	while(!neigh.Reference())
//	{
//		neigh = neigh.Neighbour();
//	}
//	if(neigh == myInnerSide)
//	{
//		/**
//		 O codigo pressupoe que os elementos computacionais 2D sao criados antes dos 1D.
//		 Quando serao criados os elementos computacionais 1D, os respectivos vizinhos 2D sao encontrados.
//		 Situacoes assim ocorrem (neste algoritmo) quando eh realizado refinamento uniforme, pois os primeiros elementos sem descendentes sao os 2D (e depois os descendentes 1D de contorno)
//		 
//		 Ocorreu o problema quando tentou-se realizar o refinamento do quadrilatero em 02 triangulos, em que o quadrilatero apresenta descendentes e as arestas nao.
//		 Neste caso a criacao de elementos computacionais eh iniciada pelos 1D, fazendo com que nao encontrem vizinhos computacionais 2D.
//		 Com isso a variavel int connectIndex0 eh setada com -1, dando o BUG observado.
//		 */
//		std::cout << "Nao foi encontrado elemento 2D com elemento computacional inicializado!!!\n"; 
//		DebugStop();
//	}
//	TPZCompElSide compneigh(neigh.Reference());
//    fneighbour = compneigh;
//	int sideoffset = neigh.Element()->NSides()-neigh.Side();
//	int neighnconnects = compneigh.Element()->NConnects();
//	int connectnumber = neighnconnects-sideoffset;
//	TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (compneigh.Element());
//	connectnumber = intel->SideConnectLocId(0,compneigh.Side());
//	int connectIndex0 = compneigh.Element()->ConnectIndex(connectnumber);
//	
//	this->fConnectIndexes[0] = connectIndex0;
//	mesh.ConnectVec()[connectIndex0].IncrementElConnected();
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		
		sout << std::endl<<" Criando Connects: "<< std::endl;
		for(int j=0; j< NConnects();j++)
		{
			sout<<" "<< this->fConnectIndexes[j];
			
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	int sideorder = SideOrder(TSHAPE::NSides-1);
	sideorder = 2*sideorder;
	if (sideorder > this->fIntRule.GetMaxOrder()) sideorder = this->fIntRule.GetMaxOrder();
	//  TPZManVector<int,3> order(3,2*sideorder+2);
	TPZManVector<int,3> order(3,sideorder);
	//TPZManVector<int,3> order(3,20);
	this->fIntRule.SetOrder(order);

#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	 {
         std::stringstream sout;
         sout << "Finalizando criacao do elemento ";
         this->Print(sout);
         LOGPZ_DEBUG(logger,sout.str())
	 }
#endif
	 
}

template<class TSHAPE>
TPZCompElSymTensorBound<TSHAPE>::TPZCompElSymTensorBound(TPZCompMesh &mesh, const TPZCompElSymTensorBound<TSHAPE> &copy) :
TPZIntelGen<TSHAPE>(mesh,copy), fSideOrient(copy.fSideOrient)
{
//	for(int i=0;i<TSHAPE::NSides;i++)
//	{
//		this-> fConnectIndexes[i] = copy.fConnectIndexes[i];
//	}
    long index = copy.fneighbour.Element()->Index();
    TPZCompEl *cel = this->Mesh()->ElementVec()[index];
    if (!cel) {
        DebugStop();
    }
    fneighbour = TPZCompElSide(cel,copy.fneighbour.Side());
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElSymTensorBound<TSHAPE>::TPZCompElSymTensorBound(TPZCompMesh &mesh,
												 const TPZCompElSymTensorBound<TSHAPE> &copy,
												 std::map<long,long> & gl2lcConMap,
												 std::map<long,long> & gl2lcElMap) :
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fSideOrient(copy.fSideOrient)
{
	
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		long lcIdx = -1;
		long glIdx = copy.fConnectIndexes[i];
		if(glIdx == -1)
		{
			// nothing to clone
			this->fConnectIndexes[i] = -1;
			continue;
		}
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end())
		{
			lcIdx = gl2lcConMap[glIdx];
		}
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
	//   gl2lcElMap[copy.fIndex] = this->Index();
    
    long neiIdx = copy.fneighbour.Element()->Index();
    if(gl2lcElMap.find(neiIdx)==gl2lcElMap.end())
    {
        DebugStop();
    }
    TPZCompEl *cel = mesh.ElementVec()[gl2lcElMap[neiIdx]];
    if (!cel) {
        DebugStop();
    }
    fneighbour = TPZCompElSide(cel,copy.fneighbour.Side());
}

// TESTADO
template<class TSHAPE>
TPZCompElSymTensorBound<TSHAPE>::TPZCompElSymTensorBound() : TPZIntelGen<TSHAPE>(),fneighbour()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
		this-> fConnectIndexes[i] = -1;
	}
}

// TESTADO
template<class TSHAPE>
TPZCompElSymTensorBound<TSHAPE>::~TPZCompElSymTensorBound(){
    TPZGeoEl *gel = this->Reference();
    if (gel->Reference() != this) {
        DebugStop();
    }
    int side = TSHAPE::NSides-1;
    TPZGeoElSide gelside(this->Reference(),side);
    TPZStack<TPZCompElSide> celstack;
    TPZCompElSide largecel = gelside.LowerLevelCompElementList2(0);
    if (largecel) {
        DebugStop();
    }
    gelside.HigherLevelCompElementList3(celstack, 0, 1);
    long ncel = celstack.size();
    if (ncel) {
        DebugStop();
    }
    gel->ResetReference();

}

// NAO TESTADO
template<class TSHAPE>
MElementType TPZCompElSymTensorBound<TSHAPE>::Type() {
	return TSHAPE::Type();
}

// NAO TESTADO
template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::SetSideOrient(int side, int sideorient)
{
    if (side != TSHAPE::NSides - 1) {
        DebugStop();
    }
    fSideOrient = sideorient;
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElSymTensorBound<TSHAPE>::GetSideOrient(int side)
{
    if (side != TSHAPE::NSides - 1) {
        DebugStop();
    }
    return fSideOrient;
}

template<class TSHAPE>
int TPZCompElSymTensorBound<TSHAPE>::NConnects() const {
	
    return TSHAPE::NSides;
}

template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::SetConnectIndex(int i, long connectindex)
{
	this->fConnectIndexes[i] = connectindex;
	
}

template<class TSHAPE>
int TPZCompElSymTensorBound<TSHAPE>::NConnectShapeF(int connect) const
{
    int order = ConnectOrder(TSHAPE::NSides-1);
    return TSHAPE::NConnectShapeF(connect,order);
}



//sets the interpolation order of side to order
template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::SetSideOrder(int side, int order) {
    if (side != TSHAPE::NSides-1) {
        DebugStop();
    }
	int connectaux= side;
	if(connectaux<0 || connectaux > this-> NConnects()) {
		PZError << "TPZCompElHDiv::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef LOG4CXX
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
	}
	TPZConnect &c = this->Connect(connectaux);
    long cindex = this->ConnectIndex(connectaux);
    c.SetOrder(order,cindex);
    long seqnum = c.SequenceNumber();
    int nvar = 2;
    int nshape = NConnectShapeF(connectaux);
    c.SetNShape(nshape);
    c.SetNState(nvar);
    this-> Mesh()->Block().Set(seqnum,nshape*nvar);
    if(connectaux == NConnects()-1)
    {
		this->SetIntegrationRule(2*order);
    }
}

/**
 return the interpolation orderof the polynomial for connect
 **/
template<class TSHAPE>
int TPZCompElSymTensorBound<TSHAPE>::ConnectOrder(int connect) const
{
	
	if (connect < 0 || connect >= this->NConnects())
	{
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Connect index out of range connect " << connect <<
			" nconnects " << NConnects();
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
		DebugStop();
#else
		std::cout << sout.str() << std::endl;
#endif
		return -1;
	}
	const TPZConnect &c = this-> Connect(connect);
	return c.Order();
}

/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
    data.fShapeType = TPZMaterialData::EScalarShape;
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		LOGPZ_DEBUG(logger,"Initializing normal vectors")
	}
#endif
    //	Nos temos que inicializar a estrutura de dados das funcoes vetoriais
    // extrai os coordenados dos nos
    // get the corner coordinates
    TPZFNMatrix<6,REAL> cornerco(2,2);
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

    // calcular o vetor tangential
    // compute the three tangent vectors
    TPZFNMatrix<6,REAL> tangent(2,1);
    for (int t=0; t<1; t++) {
        TPZManVector<REAL,2> diff(2);
        for (int c=0; c<2; c++) {
            diff[c] = cornerco(c,t+1)-cornerco(c,t);
        }
        REAL norm = sqrt(diff[0]*diff[0]+diff[1]*diff[1]);
        tangent(0,t) = diff[0]/norm;
        tangent(1,t) = diff[1]/norm;
    }
    // calcular o vetor normal
    // rotate to obtain the normal vectors
    TPZFNMatrix<6,REAL> normalvec(2,1);
    for (int n=0; n<1; n++) {
        normalvec(0,n) = tangent(1,n);
        normalvec(1,n) = -tangent(0,n);
    }
    // calcular sigman_x e sigman_y para cada tensor de borda
    // primeiro para o tensor canonico
    TPZFNMatrix<6> sigman_canonical(2,3,0.);
    sigman_canonical(0,0) = normalvec(0,0);
    
    sigman_canonical(1,1) = normalvec(1,0);
    
    sigman_canonical(0,2) = normalvec(1,0);
    sigman_canonical(1,2) = normalvec(0,0);
    // para os dois tenores de borda nn e tnt
    TPZFNMatrix<4> sigman_border(2,2);
    sigman_border(0,0) = normalvec(0);
    sigman_border(1,0) = normalvec(1);
    
    sigman_border(0,1) = tangent(0);
    sigman_border(1,1) = tangent(1);
    
    data.fNormalVec.Resize(2, 5);
    for (int i=0; i<2; i++) {
        for (int j=0; j<3; j++) {
            data.fNormalVec(i,j) = sigman_canonical(i,j);
        }
        for (int j=0; j<2; j++) {
            data.fNormalVec(i,j+3) = sigman_border(i,j);
        }
    }
    
    // initializar a estrutura vecandshape
    TPZConnect &c = this->Connect(TSHAPE::NSides-1);
    int nshape = c.NShape();
    
    int numvec = 6 + nshape*2;
    
    data.fVecShapeIndex.Resize(numvec);
    
    int count = 0;
    for (int ish = 0; ish<2; ish++) {
        for (int iv=0; iv<3; iv++) {
            data.fVecShapeIndex[count++] = std::make_pair(iv, ish);
        }
    }
    for (int ish = 0; ish < nshape; ish++) {
        for (int iv=0; iv<2; iv++) {
            data.fVecShapeIndex[count++] = std::make_pair(iv+3, ish+2);
        }
    }
    
#ifdef PZDEBUG
    if (count != numvec) {
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

template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data){
    
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
    
    if (count != multipliers.size()) {
        DebugStop();
    }
    if (count != data.fVecShapeIndex.size()) {
        DebugStop();
    }
    data.sol.resize(1);
    data.sol[0].resize(2);
    data.dsol.resize(1);
    data.dsol[0].Resize(2, 2);
    for (int i=0; i<data.fVecShapeIndex.size(); i++) {
        int ishape = data.fVecShapeIndex[i].second;
        int ivec = data.fVecShapeIndex[i].first;
        for (int j=0; j<2; j++) {
            data.sol[0][j] += data.fNormalVec(j,ivec)*data.phi(ishape)*multipliers[i];
//            for (int d=0; d<2; d++) {
//                data.dsol[0](j,d) += data.fNormalVec(j,ivec)*data.dphix(d,ishape)*multipliers[i];
//            }
        }
    }
    
}

template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<long> &shapeindex) {
	
    DebugStop();
}

/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::FirstShapeIndex(TPZVec<long> &Index){
	
	Index.Resize(TSHAPE::NSides+1);
	Index[0]=0;
    int order = ConnectOrder(0);
		
    for(int iside=0;iside<TSHAPE::NSides;iside++)
    {
        
        if(TSHAPE::Type()==EQuadrilateral){
            Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
        }
        else{
            Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
        }
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << " FirsShapeIndex result " << Index;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	
	
}

template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
    if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension || point.size() != TSHAPE::Dimension ){
		DebugStop() ;
	}
    TPZGeoEl *gel = this->Reference();
    int nc = gel->NCornerNodes();
    TPZManVector<long,8> id(nc);
    for (int ic=0; ic<nc; ic++) {
        id[ic] = gel->Node(ic).Id();
    }
    TPZManVector<int,TSHAPE::NSides> ord;
    this->GetInterpolationOrder(ord);

    TPZFNMatrix<50,REAL> philoc(phi.Rows(),phi.Cols()),dphiloc(dphi.Rows(),dphi.Cols());
    TSHAPE::Shape(point,id,ord,phi,dphi);
    
    
    return;
}

/** Compute the shape function at the integration point */
template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{
    TPZManVector<int,TSHAPE::NSides> ordl(TSHAPE::NSides,0);
    this->GetInterpolationOrder(ordl);
    int nshape = TSHAPE::NShapeF(ordl);
    phi.Resize(nshape, 1);
    dphi.Resize(TSHAPE::Dimension, nshape);
    SideShapeFunction(TSHAPE::NSides-1, pt, phi, dphi);


    if (fSideOrient == -1) {
        phi *= -1.;
        dphi *= -1.;
    }
    
    
    
    return;
}

template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data){
    
    this->Shape(intpoint, data.phi, data.dphi);
    
    TPZGeoEl *ref = this->Reference();
    ref->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
    
}

template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                         REAL &detjac, TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx){
    
    std::cout << "Method not implement call the architec." << std::endl;
    DebugStop();
    
}

/** Read the element data from a stream */
template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZIntelGen<TSHAPE>::Read(buf,context);
    buf.Read(&fSideOrient);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
	TPZIntelGen<TSHAPE>::Write(buf,withclassid);
    buf.Write(&fSideOrient);
}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    out << "Side orientation " << fSideOrient << std::endl;
    if (fRestraint.IsInitialized()) {
        fRestraint.Print(out);
    }
    TPZIntelGen<TSHAPE>::Print(out);
    
    
}

/** Returns the actual interpolation order of the polynomial along the side */
template<class TSHAPE>
int TPZCompElSymTensorBound<TSHAPE>::SideOrder(int side) const
{
	if(side == TSHAPE::NSides-1)
	{
		return ConnectOrder(side);
	}
	else {
		return -1;
	}
}

/** Return a matrix with index shape and vector associate to element */
template<class TSHAPE>
void TPZCompElSymTensorBound<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,long> > & ShapeAndVec){
	
	// VectorSide indicates the side associated with each vector entry
	TPZVec<long> FirstIndex;
	// the first index of the shape functions
	FirstShapeIndex(FirstIndex);
#ifdef LOG4CXX
		{
				std::stringstream sout;
				sout << "FirstIndex of shape functions " << FirstIndex;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	
	int tamanho= this->NShapeF();
	
	ShapeAndVec.Resize(tamanho);
	long count=0;
		//for(int jvec=0;jvec< VectorSide.NElements();jvec++)
	for(long jvec=0;jvec< VectorSide.NElements();jvec++)//coloca-se -1 caso queira reduzir o espaco de fluxo
	{
		int lside=VectorSide[jvec];
		long fshape1= FirstIndex[lside];
		long fshape2= FirstIndex[lside+1];
		for (long ishape=fshape1; ishape<fshape2; ishape++)
		{
			
#ifdef LOG4CXX
			std::stringstream sout;
			sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
			LOGPZ_DEBUG(logger,sout.str())
#endif
			
			
			ShapeAndVec[count++]=std::pair<int,long>(jvec,ishape);
		}
		
	}
	
#ifdef LOG4CXX
    std::stringstream sout;
    sout << "VecShapeIndex " << ShapeAndVec;
    LOGPZ_DEBUG(logger,sout.str())
#endif
	
	
}


#include "pzshapelinear.h"

using namespace pzshape;
template<>
int TPZCompElSymTensorBound<TPZShapeLinear>::ClassId() const
{
	return TPZSYMTENSORBOUNDLINEARID;
}

template class
TPZRestoreClass< TPZCompElSymTensorBound<TPZShapeLinear>, TPZSYMTENSORBOUNDLINEARID>;


template class TPZCompElSymTensorBound<TPZShapeLinear>;
