/**
 * @file
 * @brief Contains the TPZQuadraticLine class which defines a linear geometric element with quadratic map.
 */
#ifndef TPZQUADRATICLINE_H
#define TPZQUADRATICLINE_H

#include "TPZGeoLinear.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

namespace pzgeom {

	/**
	 * @ingroup geometry
	 * @brief Defines a linear geometric element with quadratic map. \ref geometry "Geometry"
	 * @author Paulo Cesar de Alvarenga Lucci
	 * @since 2007
	 */
	class TPZQuadraticLine : public pzgeom::TPZNodeRep<3,pztopology::TPZLine> 
    {		
        public:
		
		enum {NNodes = 3};
		
		bool IsLinearMapping() const {
			return false;
		}
		
		TPZQuadraticLine(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes)
		{
		}
		
		TPZQuadraticLine() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>()
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap)
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
		{
		}
		
		TPZQuadraticLine(const TPZQuadraticLine &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Line";} 
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, \n
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid,
										  int& index);
		
		static void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);
		
		static void X(TPZFMatrix &coord, TPZVec<REAL> &par, TPZVec<REAL> &result);
		
		static void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
	};
    
};

#endif
