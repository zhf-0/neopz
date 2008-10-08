// -*- c++ -*-

// $Id: pzelctemp.cpp,v 1.41 2008-10-08 02:13:33 phil Exp $

#include "pzelctemp.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzintelgen"));
#endif

//template<class TSHAPE>
//class TPZIntelGen : public TPZInterpolatedElement {


template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, int &index) :
  TPZInterpolatedElement(mesh,gel,index) {
  int i;
  for(i=0;i<TSHAPE::NSides-TSHAPE::NCornerNodes;i++) {
		//    fSideOrder[i] = gOrder;
		fPreferredOrder = TPZCompEl::GetgOrder();
  }
  for(i=0; i<TSHAPE::NSides; i++) fConnectIndexes[i]=-1;
  //  RemoveSideRestraintsII(EInsert);
  gel->SetReference(this);
  for(i=0;i<TSHAPE::NCornerNodes;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
  }
  for(;i<TSHAPE::NSides;i++) {
    fConnectIndexes[i] = CreateMidSideConnect(i);
    mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
    IdentifySideOrder(i);
  }

  int sideorder = SideOrder(TSHAPE::NSides-1);
  sideorder = 2*sideorder;
  if (sideorder > fIntRule.GetMaxOrder()) sideorder = fIntRule.GetMaxOrder();
  //  TPZManVector<int,3> order(3,2*sideorder+2);
  TPZManVector<int,3> order(3,sideorder);
  //TPZManVector<int,3> order(3,20);
  fIntRule.SetOrder(order);

}

template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, const TPZIntelGen<TSHAPE> &copy) :
  TPZInterpolatedElement(mesh,copy), fIntRule(copy.fIntRule) {
  fPreferredOrder = copy.fPreferredOrder;
  int i;
  for(i=0;i<TSHAPE::NSides;i++) {
    fConnectIndexes[i] = copy.fConnectIndexes[i];
  }
}


template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh,
                                      const TPZIntelGen<TSHAPE> &copy,
                                      std::map<int,int> & gl2lcConMap,
                                      std::map<int,int> & gl2lcElMap) :
    TPZInterpolatedElement(mesh,copy,gl2lcElMap), fIntRule(copy.fIntRule)
{
//   if (gl2lcElMap.find(copy.fIndex)!= gl2lcElMap.end())
//   {
//     std::stringstream sout;
//     sout << "ERROR in : " << __PRETTY_FUNCTION__
//         << " trying to clone an already cloned element index: " << copy.fIndex;
//     LOGPZ_ERROR(logger, sout.str().c_str());
//     exit(-1);
//   }

  fPreferredOrder = copy.fPreferredOrder;
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
      fConnectIndexes[i] = -1;
      return;
    }
    fConnectIndexes[i] = lcIdx;
  }
//   gl2lcElMap[copy.fIndex] = this->Index();
}


template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen() :
  TPZInterpolatedElement(), fIntRule() {
  fPreferredOrder = -1;
  int i;
  for(i=0;i<TSHAPE::NSides;i++) {
    fConnectIndexes[i] = -1;
  }
}

template<class TSHAPE>
TPZIntelGen<TSHAPE>::~TPZIntelGen(){
  if(Reference()) {
    if(Reference()->Reference()) {
      RemoveSideRestraintsII(EDelete);
    }
    Reference()->ResetReference();
  }
}

template<class TSHAPE>
MElementType TPZIntelGen<TSHAPE>::Type() {
  return TSHAPE::Type();
}


//template<class TSHAPE>
//int TPZIntelGen<TSHAPE>::NConnects() {
//  return TSHAPE::NSides;
//}

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetConnectIndex(int i, int connectindex){
#ifndef NODEBUG
  if(i<0 || i>= TSHAPE::NSides) {
    std::cout << " TPZIntelGen<TSHAPE>::SetConnectIndex index " << i <<
      " out of range\n";
    return;
  }
#endif
  fConnectIndexes[i] = connectindex;
}

template<class TSHAPE>
int TPZIntelGen<TSHAPE>::NConnectShapeF(int connect){
  if(connect < TSHAPE::NCornerNodes) return 1;
  int order = SideOrder(connect);
  if(order < 0) return 0;
  return TSHAPE::NConnectShapeF(connect, order);
}

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetIntegrationRule(int ord) {
  TPZManVector<int,3> order(TSHAPE::Dimension,ord);
  fIntRule.SetOrder(order);
}

//template<class TSHAPE>
//int TPZIntelGen<TSHAPE>::Dimension() {
//    return TSHAPE::Dimension;
//}


//template<class TSHAPE>
//int TPZIntelGen<TSHAPE>::NCornerConnects() {
//  return TSHAPE::NCornerNodes;
//}

template<class TSHAPE>
int TPZIntelGen<TSHAPE>::NSideConnects(int side){
  return TSHAPE::NSideConnects(side);
}

template<class TSHAPE>
int TPZIntelGen<TSHAPE>::SideConnectLocId(int node, int side) {
  return TSHAPE::SideConnectLocId(side,node);
}

/**Sets the interpolation order for the interior of the element*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetInterpolationOrder(int order) {
  fPreferredOrder = order;
}

/**Identifies the interpolation order on the interior of the element*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
  int i;
  for(i=0; i<TSHAPE::NSides-TSHAPE::NCornerNodes; i++) {
    ord[i] = SideOrder(i+TSHAPE::NCornerNodes);
  }
}

/**return the preferred order of the polynomial along side iside*/
template<class TSHAPE>
int TPZIntelGen<TSHAPE>::PreferredSideOrder(int side) {
  if(side < TSHAPE::NCornerNodes) return 0;
  if(side<TSHAPE::NSides) {
	  int order =fPreferredOrder;
	  return AdjustPreferredSideOrder(side,order);
  }
  PZError << "TPZIntelgen::PreferredSideOrder called for side = " << side << "\n";
  return 0;

}

template<class TSHAPE>
int TPZIntelGen<TSHAPE>::ConnectIndex(int con) const{

#ifndef NODEBUG
  if(con<0 || con>= TSHAPE::NSides) {
    std::cout << "TPZIntelgen::ConnectIndex wrong parameter con " << con <<
      " NSides " << TSHAPE::NSides << std::endl;
    return -1;
  }

#endif
  return fConnectIndexes[con];
}



/**Sets the preferred interpolation order along a side
   This method only updates the datastructure of the element
   In order to change the interpolation order of an element, use the method PRefine
*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetPreferredOrder(int order)
{
  fPreferredOrder = order;
}

/**sets the interpolation order of side to order*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetSideOrder(int side, int order) {
  if(side<0 || side >= TSHAPE::NSides || (side >= TSHAPE::NCornerNodes && order <1)) {
    PZError << "TPZIntelGen::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef LOG4CXX
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
    LOGPZ_DEBUG(logger,sout.str())
#endif
    return;
  }
  if(side>= TSHAPE::NCornerNodes) {
    if(fConnectIndexes[side] == -1) return;
    TPZConnect &c = Connect(side);
    c.SetOrder(order);
    int seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZAutoPointer<TPZMaterial> mat = Material();
    if(mat) nvar = mat->NStateVariables();
    Mesh()->Block().Set(seqnum,NConnectShapeF(side)*nvar);
    if(side == TSHAPE::NSides-1) {
      SetIntegrationRule(2*order);
    }
  }
}

/**returns the actual interpolation order of the polynomial along the side*/
template<class TSHAPE>
int TPZIntelGen<TSHAPE>::SideOrder(int side) {
  if(side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) return 0;
  if(fConnectIndexes[side] == -1)
  {

    std::stringstream sout ;
    sout << __PRETTY_FUNCTION__ << " side " << side << std::endl;
    //Print(sout);
#ifdef LOG4CXX
    LOGPZ_ERROR(logger,sout.str());
    DebugStop();
#else
    std::cout << sout.str() << std::endl;
#endif
    return -1;
  }
  TPZConnect &c = Connect(side);
  return c.Order();
}

/**transform a point in the parameter space of the side into a point in the space
   of the master element*/
//template<class TSHAPE>
//void TPZIntelGen<TSHAPE>::SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {

//}

/**transform a point in the parameter space of the master element into a point in the
   space of the side*/
//virtual void ElementToSideParameter(int side, TPZVec<REAL> &point, TPZVec<REAL> &par);


/**compute the values of the shape function of the side*/

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {

  int nc = TSHAPE::NSideConnects(side);
  int nn = TSHAPE::NSideNodes(side);
  TPZManVector<int,27> id(nn),order(nc-nn);
  int n,c;
  TPZGeoEl *ref = Reference();
  for (n=0;n<nn;n++){
    int nodloc = TSHAPE::SideNodeLocId(side,n);
    id [n] = ref->NodePtr(nodloc)->Id();
  }
  for (c=nn;c<nc;c++){
    int conloc = TSHAPE::SideConnectLocId(side,c);
    order[c-nn] = SideOrder(conloc);
  }
  TSHAPE::SideShape(side, point, id, order, phi, dphi);
}

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  TPZManVector<int,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
  TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
  int i;
  TPZGeoEl *ref = Reference();
  for(i=0; i<TSHAPE::NCornerNodes; i++) {
    id[i] = ref->NodePtr(i)->Id();
  }
  for(i=0; i<TSHAPE::NSides-TSHAPE::NCornerNodes; i++) {
    ord[i] = SideOrder(i+TSHAPE::NCornerNodes);
  }
  TSHAPE::Shape(pt,id,ord,phi,dphi);
}


// template<class TSHAPE>
// void TPZIntelGen<TSHAPE>::Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol) {

// }

  /** Jorge 09/06/2001
   * Returns the transformation which transform a point from the side to the interior of the element
   */
template<class TSHAPE>
TPZTransform TPZIntelGen<TSHAPE>::TransformSideToElement(int side){
  return TSHAPE::TransformSideToElement(side);
}

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
/*
template<class TSHAPE>
int TPZIntelGen<TSHAPE>::ClassId() const
{
  std::cout << "TPZIntelGen<TSHAPE>::ClassId() type not specified at " << __FILE__ << ':' << __LINE__ << std::endl;
  return -1;
}
*/
  /**
  Save the element data to a stream
  */
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
  TPZInterpolatedElement::Write(buf,withclassid);
  TPZManVector<int,3> order(3,0);
  fIntRule.GetOrder(order);
  WriteObjects(buf,order);
  buf.Write(fConnectIndexes,TSHAPE::NSides);
  buf.Write(&fPreferredOrder,1);
  int classid = this->ClassId();
  buf.Write ( &classid, 1 );
}

  /**
  Read the element data from a stream
  */
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::Read(TPZStream &buf, void *context)
{
  TPZInterpolatedElement::Read(buf,context);
  TPZManVector<int,3> order;
  ReadObjects(buf,order);
  fIntRule.SetOrder(order);
  buf.Read(fConnectIndexes,TSHAPE::NSides);
  buf.Read(&fPreferredOrder,1);
  int classid = -1;
  buf.Read( &classid, 1 );
  if ( classid != this->ClassId() )
  {
    std::stringstream sout;
    sout << "ERROR - " << __PRETTY_FUNCTION__
        << " trying to restore an object id " << this->ClassId() << " for an package of id = " << classid;
    LOGPZ_ERROR ( logger, sout.str().c_str() );
  }
}

#include "tpzint1point.h"
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
#include "tpzpyramid.h"
/*
void TPZIntelGen<TPZPoint,TPZShapePoint>::Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  phi(0,0) = 1.;
}
*/
using namespace pzgeom;
using namespace pzshape;

template<>
void TPZIntelGen<TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 0) std::cout << "A point element has no graphical representation\n";
}



template<class TSHAPE>
void TPZIntelGen<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == TSHAPE::Dimension && Material()->Id() > 0) {
    new typename TSHAPE::GraphElType(this,&grafgrid);
  }
}

template<>
int TPZIntelGen<TPZShapePoint>::ClassId() const
{
  return TPZINTELPOINTID;
}

#ifndef WIN32
template class
		TPZRestoreClass< TPZIntelGen<TPZShapePoint>, TPZINTELPOINTID>;
#endif


template<>
int TPZIntelGen<TPZShapeLinear>::ClassId() const
{
  return TPZINTELLINEARID;
}

#ifndef WIN32
template class
		TPZRestoreClass< TPZIntelGen<TPZShapeLinear>, TPZINTELLINEARID>;
#endif

template<>
int TPZIntelGen<TPZShapeTriang>::ClassId() const
{
	return TPZINTELTRIANGLEID;
}

#ifndef WIN32
template class
		TPZRestoreClass< TPZIntelGen<TPZShapeTriang>, TPZINTELTRIANGLEID>;
#endif

template<>
int TPZIntelGen<TPZShapeQuad>::ClassId() const
{
	return TPZINTELQUADID;
}

#ifndef WIN32
template class
		TPZRestoreClass< TPZIntelGen<TPZShapeQuad>, TPZINTELQUADID>;
#endif

template<>
int TPZIntelGen<TPZShapeCube>::ClassId() const
{
	return TPZINTELCUBEID;
}

#ifndef WIN32
template class
		TPZRestoreClass< TPZIntelGen<TPZShapeCube>, TPZINTELCUBEID>;
#endif

template<>
int TPZIntelGen<TPZShapeTetra>::ClassId() const
{
	return TPZINTELTETRAID;
}

#ifndef WIN32
template class
		TPZRestoreClass< TPZIntelGen<TPZShapeTetra>, TPZINTELTETRAID>;
#endif

template<>
int TPZIntelGen<TPZShapePrism>::ClassId() const
{
	return TPZINTELPRISMID;
}

#ifndef WIN32
template class
		TPZRestoreClass< TPZIntelGen<TPZShapePrism>, TPZINTELPRISMID>;
#endif

template<>
int TPZIntelGen<TPZShapePiram>::ClassId() const
{
	return TPZINTELPYRAMID;
}

#ifndef WIN32
template class
		TPZRestoreClass< TPZIntelGen<TPZShapePiram>, TPZINTELPYRAMID>;
#endif


template class TPZIntelGen<TPZShapeTriang>;
template class TPZIntelGen<TPZShapePoint>;
template class TPZIntelGen<TPZShapeLinear>;
template class TPZIntelGen<TPZShapeQuad>;
template class TPZIntelGen<TPZShapeTetra>;
template class TPZIntelGen<TPZShapePrism>;
template class TPZIntelGen<TPZShapePiram>;
template class TPZIntelGen<TPZShapeCube>;

