/**
 * @file TPZMatMFHCurlH1.h
 * @brief Header file for class TPZMatMFHCurlH1.\n
 */

#ifndef TPZMATMFHCURLH1_H
#define TPZMATMFHCURLH1_H

#include "TPZVecL2.h"
#include "pzaxestools.h"
#include "pzvec_extras.h"
#include "../HCurl2D/TPZMatHCurl2D.h"

enum whichMatrix { NDefined = 0 , A = 1 , B = 2};
/**
 * @ingroup material
 * @brief This class implements the weak statement of the model problem from Oden's book, Chapter 1, within the PZ environment
 */
class  TPZMatMFHCurlH1 : public TPZVecL2
{
    
protected:
    
    //COM CERTEZA
    STATE (*fUr)( const TPZVec<REAL>&);
    STATE (*fEr)( const TPZVec<REAL>&);
    REAL fLambda;
    REAL fKz;
    REAL fE0;
    REAL fW;
    REAL fTheta;
    REAL fScale;
    whichMatrix assembling;
    
    
public:
    
    TPZMatMFHCurlH1(int id, REAL lambda, REAL kz , STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &), REAL e0, REAL t, REAL scale);
    
    TPZMatMFHCurlH1(int id);
    
    /** @brief Default constructor */
    TPZMatMFHCurlH1();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatMFHCurlH1(const TPZMatMFHCurlH1 &mat);
    /** @brief Default destructor */
    virtual ~TPZMatMFHCurlH1();
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatMFHCurlH1"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const {return 2;}
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() { return 1;}
    
public:
    /**
     * @brief Sets Matrix A for assembling
     * @details This material is designed for solving the
     * generalised eigenvalue problem stated as Ax = lBx
     * Matrices A and B are assembled separatedly.
     */
    virtual void SetMatrixA(){ assembling = A;};
    /**
     * @brief Sets Matrix B for assembling
     * @details This material is designed for solving the
     * generalised eigenvalue problem stated as Ax = lBx
     * Matrices A and B are assembled separatedly.
     */
    virtual void SetMatrixB(){ assembling = B;};
    
    virtual void ContributeValidateFunctions(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    void ComputeCurl(TPZFMatrix<REAL> gradScalarPhi , TPZFMatrix<REAL> ivecHCurl , TPZFMatrix<REAL> &curlPhi );
    void RotateForHCurl(TPZVec<REAL> normal , TPZFMatrix<REAL> vHdiv , TPZFMatrix<REAL> &vHcurl );
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeForcingRTBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * to multiphysics simulation.
     * @param datavec [in]  stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 18, 2011
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the residual vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * @brief This method defines which parameters need to be initialized in order to compute the contribution of an element
     * @param datavec [out] vector of TPZMaterialData, each position will specifie the requirements for its correspondent state variable
     */
    virtual void FillDataRequirements(TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        data.fNeedsNormal = true;
    }
    
    
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        data.fNeedsNormal = true;
    }
    
    /**
     * @brief This method defines which parameters need to be initialized in order to compute the contribution of an element
     * @param datavec [out] vector of TPZMaterialData, each position will specifie the requirements for its correspondent state variable
     */
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
    {
        int nref = datavec.size();
        for(int iref = 0; iref<nref; iref++){
            datavec[iref].SetAllRequirements(false);
            datavec[iref].fNeedsNormal = true;
        }
    }
    
    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const;
    
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
    {
        // default is no specific data requirements
        int nref = datavec.size();
        for(int iref = 0; iref<nref; iref++){
            datavec[iref].SetAllRequirements(false);
            datavec[iref].fNeedsNormal = true;
        }
    }
    
    virtual int VariableIndex(const std::string &name);
    
    /**
     * @brief Returns the number of variables associated with the variable indexed by var.
     * @param var Index variable into the solution, is obtained by calling VariableIndex
     */
    virtual int NSolutionVariables(int var);
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
};

#endif

