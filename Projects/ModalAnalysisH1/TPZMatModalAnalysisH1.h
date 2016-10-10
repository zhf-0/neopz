/**
 * @file TPZMatModalAnalysisH1.h
 * @brief Header file for class TPZMatModalAnalysisH1.\n
 */

#ifndef TPZMATMODALANALYSISH1_H
#define TPZMATMODALANALYSISH1_H

#include "pzaxestools.h"
#include "pzmaterial.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement of a waveguide problem as stated in Jin's 
 * The Finite Element Method in Electromagnetics (chapter 8 of 3rd edition).
 * It used a 2D Hcurl space for the transversal components of the electric field and an 1D
 * H1 space for the longitudinal component.
 */

enum modeType{ NDefined = -1 , modesTE=0, modesTM=1};

class  TPZMatModalAnalysisH1 : public TPZMaterial
{
    
protected:
    enum whichMatrix { NDefined = 0 , A = 1 , B = 2};
    //COM CERTEZA
    STATE (*fUr)( const TPZVec<REAL>&);
    STATE (*fEr)( const TPZVec<REAL>&);
    REAL fW;
    whichMatrix assembling;
    modeType whichMode;
    STATE fGammaZ;
public:
    
    TPZMatModalAnalysisH1(int id, REAL freq, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &));
    
    TPZMatModalAnalysisH1(int id);
    
    /** @brief Default constructor */
    TPZMatModalAnalysisH1();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatModalAnalysisH1(const TPZMatModalAnalysisH1 &mat);
    /** @brief Default destructor */
    virtual ~TPZMatModalAnalysisH1();
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatModalAnalysisH1"; }
    
    virtual void SetGammaZ(const STATE val){ fGammaZ = val; }
    
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
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
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
    }
    
    
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
    }
    
    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const;
    
    virtual int VariableIndex(const std::string &name);
    
    /**
     * @brief Returns the number of variables associated with the variable indexed by var.
     * @param var Index variable into the solution, is obtained by calling VariableIndex
     */
    virtual int NSolutionVariables(int var);
    
//    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
};

#endif

