/**
 * @file
 * @brief Contains declaration of TPZCompElSymTensor class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCOMPELSYMTENSORH
#define TPZCOMPELSYMTENSORH

#include "pzelctemp.h"
#include "TPZOneShapeRestraint.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElSymTensor : public TPZIntelGen<TSHAPE> {
	
    /// vector which defines whether the normal is outward or not
    TPZManVector<int, TSHAPE::NFaces> fSideOrient;
    
    /** @brief To append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);
public:
	

public:
    
	TPZCompElSymTensor(TPZCompMesh &mesh, TPZGeoEl *gel, long &index);
	
	TPZCompElSymTensor(TPZCompMesh &mesh, const TPZCompElSymTensor<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElSymTensor(TPZCompMesh &mesh,
				  const TPZCompElSymTensor<TSHAPE> &copy,
				  std::map<long,long> & gl2lcConMap,
				  std::map<long,long> & gl2lcElMap);
	
	TPZCompElSymTensor();
	
	virtual ~TPZCompElSymTensor();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZCompElSymTensor<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<long,long> & gl2lcConMap,std::map<long,long>&gl2lcElMap) const
	{
		return new TPZCompElSymTensor<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh){
        mesh->ApproxSpace().SetAllCreateFunctionsSymTensor(TSHAPE::Dimension);
	}
	
    /** @brief Prints the relevant data of the element to the output stream */
	virtual void Print(std::ostream &out = std::cout) const;
	

	
	virtual MElementType Type();
	
	
	virtual void SetConnectIndex(int i, long connectindex);
	
    /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect) const;
	
	virtual int Dimension() const {
		return TSHAPE::Dimension;
	}
	
    /**return the first shape associate to each side*/
    void FirstShapeIndex(TPZVec<int> &Index) const;

	/**
     * @brief return the number of shape for flux(just for flux)
	 **/
	virtual int NFluxShapeF() const;
    
    /**
     * @brief It returns the normal orientation of the reference element by the side.
     * Only side that has dimension larger than zero and smaller than me.
     * @param side: side of the reference elemen
     */
    virtual int GetSideOrient(int side);
    
    /**
     * @brief It set the normal orientation of the element by the side.
     * Only side that has dimension equal to my dimension minus one.
     * @param side: side of the reference elemen
     */
    virtual void SetSideOrient(int side, int sideorient);
    
	
	virtual void SetIntegrationRule(int ord);
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord);
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside);
	
	/*
     * @brief Sets the preferred interpolation order along a side \n
	 * This method only updates the datastructure of the element
	 * In order to change the interpolation order of an element, use the method PRefine
	 */
	virtual void SetPreferredOrder(int order);
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order);
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int SideOrder(int side) const;
	
    /**
     * @brief return the interpolation order of the polynomial for connect
     **/
	virtual int ConnectOrder(int connect) const;
    
	/**
     * @brief return the number of continuous functions 
     **/
	int NShapeContinuous(TPZVec<int> &order);
    
    /// Return the maximum order??
    virtual int MaxOrder();
    
    /// the orientation of the face
    int SideOrient(int face)
    {
#ifdef PZDEBUG
        if (face < 0 || face >= TSHAPE::NFaces) {
            DebugStop();
        }
#endif
        return fSideOrient[face];
    }
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data);
    
	/** @brief Compute and fill data with requested attributes */
	virtual void ComputeRequiredData(TPZMaterialData &data,
									 TPZVec<REAL> &qsi);

	/** @brief Compute the correspondence between the normal vectors and the shape functions */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<long> &shapeindex);
	
	/** 
	 * @brief Returns the vector index  of the first index shape associate to to each side 
	 * Special implementation to Hdiv
	 */
	void FirstShapeIndex(TPZVec<long> &Index) const;
    
	/**
     * @brief Returns a matrix index of the shape and vector  associate to element
     * @param[in] VectorSide Indicates the side associated with each vector
     * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
	 * @param[in] pressureorder Order of the pressure (to select shape functions?)
	 */
	void IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,long> > & IndexVecShape, int pressureorder);
	
	/**
     * @brief Returns a matrix index of the shape and vector  associate to element
     * @param[in] VectorSide Indicates the side associated with each vector
     * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
	 * @param[in] pressureorder Order of the pressure (to select shape functions?)
	 */
	void IndexShapeToVec(TPZVec<int> &VectorSide, TPZVec<int> &bilinear, TPZVec<int> &direction, TPZVec<std::pair<int,long> > & IndexVecShape, int pressureorder);
    void IndexShapeToVec2(TPZVec<int> &VectorSide, TPZVec<int> &bilinear, TPZVec<int> &direction, TPZVec<std::pair<int,long> > & IndexVecShape, int pressureorder);

	/** @brief Computes the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
    
    /** @brief Compute the solution for a given variable */
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol);
	
public:
    
    /** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] qsi point in master element coordinates 
	 * @param[in] data stores all input data
	 */
    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data);
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
	
	
	/** Jorge 09/06/2001
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
    /** @brief Refinement along the element */
    virtual void PRefine(int order);
	
};

/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateSymTensorTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateSyMTensorBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);

/** @} */

#endif
