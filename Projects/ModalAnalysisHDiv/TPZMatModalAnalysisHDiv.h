/**
 * @file TPZMatModalAnalysisHDiv.h
 * @brief Header file for class TPZMatModalAnalysisHDiv.\n
 */

#ifndef TPZMATMODALANALYSISHDIV_H
#define TPZMATMODALANALYSISHDIV_H

#include "pzaxestools.h"
#include "pzmaterial.h"
#include "TPZVecL2.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement of a waveguide problem as stated in Jin's 
 * The Finite Element Method in Electromagnetics (chapter 8 of 3rd edition).
 * It used a 2D Hcurl space for the transversal components of the electric field and an 1D
 * H1 space for the longitudinal component.
 */

enum modeType{ NDefined = -1 , modesTE=0, modesTM=1};

class  TPZMatModalAnalysisHDiv : public TPZVecL2
{
    
protected:
    enum whichMatrix { NDefined = 0 , A = 1 , B = 2};
    //COM CERTEZA
    STATE (*fUr)( const TPZVec<REAL>&);
    STATE (*fEr)( const TPZVec<REAL>&);
    REAL fW;
    whichMatrix assembling;
    modeType whichMode;
    STATE fKtSquared;
    
    const int l2index = 0;
    const int hdivindex = 1;
public:
    
    TPZMatModalAnalysisHDiv(int id, REAL freq, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &));
    
    TPZMatModalAnalysisHDiv(int id);
    
    /** @brief Default constructor */
    TPZMatModalAnalysisHDiv();
    
    int L2Index(){ return l2index; }
    int HDivIndex(){ return hdivindex; }
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatModalAnalysisHDiv(const TPZMatModalAnalysisHDiv &mat);
    /** @brief Default destructor */
    virtual ~TPZMatModalAnalysisHDiv();
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatModalAnalysisHDiv"; }
    
    virtual void SetKtSquared(const STATE val){ fKtSquared = val; }
    
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
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
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
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
};

#endif