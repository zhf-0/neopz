/**
 * @file
 * @brief Contains the implementation of the TPZHCurlNedFTriEl methods.
 */

#include "TPZHCurlNedFTriEl.h"
#include "pzshapetriang.h"
#include "pzmaterial.h"

using namespace pzshape;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZHCurlNedFTriEl"));
#endif

TPZHCurlNedFTriEl::TPZHCurlNedFTriEl(TPZCompMesh &mesh, TPZGeoEl *geoEl, long &index) :
TPZInterpolatedElement(mesh,geoEl,index),
fConnectIndexes(TPZShapeTriang::NFaces + 1, -1),
fSideOrient(TPZShapeTriang::NFaces, 1){
    
    geoEl->SetReference(this);
    
    TPZStack<int> facesides;
    std::cout<<TPZShapeTriang::NSides-1<<std::endl;
    TPZShapeTriang::LowerDimensionSides(TPZShapeTriang::NSides-1,facesides,TPZShapeTriang::Dimension-1);//inserts 3 edges
    facesides.Push(TPZShapeTriang::NSides-1);//inserts triangle face
    for(int i=0;i< fConnectIndexes.size() ;i++) {
        int sideIndex = facesides[i];
        fConnectIndexes[i] = CreateMidSideConnect( sideIndex );
        mesh.ConnectVec()[ fConnectIndexes[i] ].IncrementElConnected();
        this->IdentifySideOrder(sideIndex);
    }
    
    AdjustIntegrationRule();
    
    int firstside = TPZShapeTriang::NSides-TPZShapeTriang::NFaces-1;
    for(int side = firstside ; side < TPZShapeTriang::NSides-1; side++ )
    {
        fSideOrient[side-firstside] = this->Reference()->NormalOrientation(side);
    }
}

TPZHCurlNedFTriEl::TPZHCurlNedFTriEl(TPZCompMesh &mesh, const TPZHCurlNedFTriEl &copy) :
TPZInterpolatedElement(mesh,copy), fConnectIndexes(copy.fConnectIndexes),
fSideOrient(copy.fSideOrient){
    
    fPreferredOrder = copy.fPreferredOrder;
}

TPZHCurlNedFTriEl::TPZHCurlNedFTriEl(TPZCompMesh &mesh,
                                     const TPZHCurlNedFTriEl &copy,
                                     std::map<long,long> & gl2lcElMap) :
TPZInterpolatedElement(mesh,copy,gl2lcElMap),
fConnectIndexes(copy.fConnectIndexes),
fSideOrient(copy.fSideOrient)
{
}


TPZHCurlNedFTriEl::TPZHCurlNedFTriEl() :
TPZInterpolatedElement(),
fConnectIndexes(TPZShapeTriang::NSides -1){
    fPreferredOrder = -1;
    for( int i = 0 ; i< fConnectIndexes.size() ; i++) {
        fConnectIndexes[i] = -1;
    }
}

TPZHCurlNedFTriEl::~TPZHCurlNedFTriEl() {
}

TPZHCurlNedFTriEl  * TPZHCurlNedFTriEl::Clone(TPZCompMesh &mesh) const{
    DebugStop();
}

TPZHCurlNedFTriEl  * TPZHCurlNedFTriEl::ClonePatchEl(TPZCompMesh &mesh, std::map<long,long> & gl2lcConMap, std::map<long,long> & gl2lcElMap) const{
    DebugStop();
}

int TPZHCurlNedFTriEl::Dimension() const{
    return TPZShapeTriang::Dimension;
}

int TPZHCurlNedFTriEl::NCornerConnects() const{
    return 0;
}

int TPZHCurlNedFTriEl::NConnects() const{
    return fConnectIndexes.size();
}

long TPZHCurlNedFTriEl::ConnectIndex(int iCon) const{
#ifndef NODEBUG
    if( iCon < 0 || iCon >= this->NConnects()) {
        std::cout << "TPZHCurlNedFTriEl::ConnectIndex wrong parameter connect " << iCon <<
        " NConnects " << this-> NConnects() << std::endl;
        DebugStop();
        return -1;
    }
#endif
    
    return this->fConnectIndexes[iCon];
}

void TPZHCurlNedFTriEl::SetConnectIndex(int i, long connectindex){
#ifdef PZDEBUG
    if(i<0 || i>= this->NConnects()) {
        std::cout << "TPZHCurlNedFTriEl::SetConnectIndex index " << i <<
        " out of range\n";
        DebugStop();
    }
#endif
    this-> fConnectIndexes[i] = connectindex;
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << std::endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
}

int TPZHCurlNedFTriEl::NSideConnects(int side) const{
    if(TPZShapeTriang::SideDimension(side) <= Dimension()-2) return 0;
    if(TPZShapeTriang::SideDimension(side) == Dimension()-1) return 1; //side connects
    if(TPZShapeTriang::SideDimension(side) == Dimension()) return 1; //internal connects
    return -1;
}// TODO: TRANSFER TO LINEAR EL

int TPZHCurlNedFTriEl::SideConnectLocId(int iCon,int iSide) const{
#ifdef PZDEBUG
    if( TPZShapeTriang::SideDimension(iSide) <= TPZShapeTriang::Dimension - 2 || iCon >= NSideConnects(iSide) ){
        PZError << "TPZHCurlNedFTriEl::SideConnectLocId no connect associate " <<  std::endl;
        return -1;
    }
#endif
    return iSide - (TPZShapeTriang::NSides-TPZShapeTriang::NumSides(TPZShapeTriang::Dimension-1)-1);
}

int TPZHCurlNedFTriEl::NConnectShapeF(int icon, int order) const{
    int side = icon + TPZShapeTriang::NSides-TPZShapeTriang::NumSides(TPZShapeTriang::Dimension-1)-1 ;
    if(TPZShapeTriang::SideDimension(side) == Dimension()-1) return EffectiveSideOrder(side);
    if(TPZShapeTriang::SideDimension(side) == Dimension()){
        return (this->fPreferredOrder *(this->fPreferredOrder - 1));
    }
    return -1;
}

void TPZHCurlNedFTriEl::Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi){
    DebugStop();    
}

void TPZHCurlNedFTriEl::SideShapeFunction(int side, TPZVec<REAL> &point, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi){
    DebugStop();
}

void TPZHCurlNedFTriEl::SetIntegrationRule(int ord){
    TPZManVector<int,3> order(TPZShapeTriang::Dimension,ord);
    fIntRule.SetOrder(order);
}
const TPZIntPoints &TPZHCurlNedFTriEl::GetIntegrationRule() const{
    if (this->fIntegrationRule) {
        return *fIntegrationRule;
    }
    else
    {
        return fIntRule;
    }
}

TPZIntPoints &TPZHCurlNedFTriEl::GetIntegrationRule(){
    if (fIntegrationRule) {
        return *fIntegrationRule;
    }
    else
    {
        return fIntRule;
    }
}

void TPZHCurlNedFTriEl::SetPreferredOrder(int order){
    fPreferredOrder = order;
}

void TPZHCurlNedFTriEl::GetInterpolationOrder(TPZVec<int> &ord){
    DebugStop();
}

int TPZHCurlNedFTriEl::PreferredSideOrder(int side){
    if(TPZShapeTriang::SideDimension(side) < Dimension()-1)
    {
        DebugStop();
    }
    int connect= SideConnectLocId(0,side);
    if(connect < 0 || connect > NConnects()) {
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
    if(connect < NConnects()) {
        int order =this->fPreferredOrder;
        return order;
    }
    PZError << "TPZHCurlNedFTriEl::PreferredSideOrder called for connect = " << connect << "\n";
    return 0;
}

int TPZHCurlNedFTriEl::ConnectOrder(int connect) const{
    if (connect < 0 || connect >= this->NConnects()){
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "Connect index out of range connect " << connect <<
            " nconnects " << NConnects();
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        DebugStop();
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
        DebugStop();
        return 0;
    }
    
    TPZConnect &c = this->Connect(connect);
    return c.Order();
}

int TPZHCurlNedFTriEl::EffectiveSideOrder(int side) const{
    
    if(!NSideConnects(side)) return -1;
    int firstElCon = SideConnectLocId(0, side);
    int maxorder = 0;
    int conectaux;
    if(firstElCon >= 0 || firstElCon <= NConnects()) return ConnectOrder(firstElCon);
    
    DebugStop();
    return -1;
}

void TPZHCurlNedFTriEl::SetSideOrder(int side, int order){
    
    int connectaux = SideConnectLocId(0,side);
    if(connectaux < 0 || connectaux > this-> NConnects()) {
#ifdef LOG4CXX
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
        LOGPZ_DEBUG(logger,sout.str())
#endif
        DebugStop();
        return;
    }
    TPZConnect &c = this->Connect(connectaux);
    c.SetOrder(order,this->fConnectIndexes[connectaux]);
    long seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial * mat =this-> Material();
    if(mat) nvar = mat->NStateVariables();
    c.SetNState(nvar);
    int nshape =this->NConnectShapeF(connectaux,order);
    c.SetNShape(nshape);
    this-> Mesh()->Block().Set(seqnum,nshape*nvar);
}

TPZTransform TPZHCurlNedFTriEl::TransformSideToElement(int side){
    DebugStop();
}

TPZCompEl * CreateHCurlNedFLinEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    DebugStop();
    return new TPZHCurlNedFTriEl(mesh,gel,index);
}

TPZCompEl * CreateHCurlNedFTriEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
    return new TPZHCurlNedFTriEl(mesh,gel,index);
}
