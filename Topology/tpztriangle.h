//
// C++ Interface: tpztriangle
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZTOPOLOGYTPZTRIANGLE_H
#define PZTOPOLOGYTPZTRIANGLE_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzeltype.h"

class TPZIntPoints;
class TPZGraphElTd;
class TPZIntTriang;
class TPZGraphElT2dMapped;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

namespace pztopology {

/**
This class defines the topology of a triangle

@author Philippe R. B. Devloo
*/
	
/// This class defines the topology of the Triangle element
class TPZTriangle{
public:

enum {NSides = 7, NCornerNodes= 3, Dimension = 2};

TPZTriangle(){
};

virtual ~TPZTriangle(){
};

 /**
  * return the number of nodes (not connectivities) associated with a side
  */
static int NSideNodes(int side);
 /**
  * returns the local node number of the node "node" along side "side"
  */
static int SideNodeLocId(int side, int node);
	
/**
 * returns the local node number side "side" is particular method for hdivspace
 */
//static int SideNodeLocId(int side);
	
	
	
 /**
  * return the number of nodes (not connectivities) associated with a side
  */
 /**
  * returns the dimension of the side
  */
static int SideDimension(int side);
	/**
	 * return the number of connects for a set dimension
	 */
static int NumSides(int dimension);
	
 /**
  * returns the transformation which takes a point from the side sidefrom ot
  * the side sideto
  * @param sidefrom side where the point resides
  * @param sideto side whose closure contains sidefrom
  */
static TPZTransform SideToSideTransform(int sidefrom, int sideto);
/**
 * Returns the transformation which transform a point from the interior of the element to the side
 * @param side side to which the point will be tranformed (0<=side<=8)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformElementToSide(int side);
/**
 * Returns the transformation which transform a point from the side to the interior of the element
 * @param side side from which the point will be tranformed (0<=side<=2)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformSideToElement(int side);

/** Verifies if the parametric point pt is in the element parametric domain
 */
static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6);

static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);

 /**
  * returns all sides whose closure contains side
  * @param side smaller dimension side
  * @param high vector which will contain all sides whose closure contain sidefrom
  */
static void HigherDimensionSides(int side, TPZStack<int> &high);
/**
 * returns the barycentric coordinates in the master element space of the original element
 */

 static void CenterPoint(int side, TPZVec<REAL> &center);

 /**volume of the master element*/
static REAL RefElVolume(){return 0.5;}

  /**
   * Create an integration rule 
   * @param order order of the integration rule to be created
   * @param side side to create integration rule
   */
static TPZIntPoints * CreateSideIntegrationRule(int side, int order);

  typedef TPZIntTriang IntruleType;
  typedef TPZGraphElT2dMapped GraphElType;


  /**
   * return the type of the element as specified in file pzeltype.h
   */
 static MElementType Type();// { return ETriangle;}

/**
  * return the type of the element as specified in file pzeltype.h
  */
static MElementType Type(int side) ;


/**
 * Number of connects of the element (7)
 * @return number of connects of the element
 */
static int NumSides();


static int NContainedSides(int side);
	/**
	 * returns the local connect along side "side" especial for hdivspace
	 */
	static int ContainedSideLocId(int side);

	/// function pointer which determines the type of computational element
 /**
	 * function pointer which determines what type of computational element will be created
  * @returns the local connect number of the connect "c" along side "side"
  */
	static TPZCompEl *(*fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index);
static int ContainedSideLocId(int side, int c);
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	static int GetTransformId(TPZVec<int> &id);
	
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	static int GetTransformId(int side, TPZVec<int> &id);
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	static void GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather);
	
	
};

}

#endif
