/**
 * @file TPZMatInhomogeneousModalAnalysisH1.h
 * @brief Header file for class TPZMatInhomogeneousModalAnalysisH1.\n
 */

#ifndef TPZMATMODALANALYSISH1H1_H
#define TPZMATMODALANALYSISH1H1_H

#include "pzmaterial.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement of a waveguide problem as stated in Jin's
 * The Finite Element Method in Electromagnetics (chapter 8 of 3rd edition).
 * It used a 2D Hcurl space for the transversal components of the electric field and an 1D
 * H1 space for the longitudinal component.
 */
class  TPZMatInhomogeneousModalAnalysisH1 : public TPZMaterial
{

protected:
    enum whichMatrix { NDefined = 0 , A = 1 , B = 2};
    STATE (*fUr)( const TPZVec<REAL>&);
    STATE (*fEr)( const TPZVec<REAL>&);
    REAL fKzOverK0;
    whichMatrix fAssembling;
    const int fEzIndex = 0;
    const int fHzIndex = 1;
    bool fIsTesting;

public:

    TPZMatInhomogeneousModalAnalysisH1(int id, REAL kzOverk0, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &));


    /** @brief Default constructor */
    TPZMatInhomogeneousModalAnalysisH1();

    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatInhomogeneousModalAnalysisH1(const TPZMatInhomogeneousModalAnalysisH1 &mat);
    /** @brief Default destructor */
    virtual ~TPZMatInhomogeneousModalAnalysisH1();

    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatInhomogeneousModalAnalysisH1"; }

    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const {return 2;}

    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() { return 1;}

    int EzIndex() const { return fEzIndex;}
    int HzIndex() const { return fHzIndex;}

public:
    /**
     * @brief Sets Matrix A for assembling
     * @details This material is designed for solving the
     * generalised eigenvalue problem stated as Ax = lBx
     * Matrices A and B are assembled separatedly.
     */
    virtual void SetMatrixA(){ fAssembling = A;};
    /**
     * @brief Sets Matrix B for assembling
     * @details This material is designed for solving the
     * generalised eigenvalue problem stated as Ax = lBx
     * Matrices A and B are assembled separatedly.
     */
    virtual void SetMatrixB(){ fAssembling = B;};

    virtual void ContributeValidateFunctions(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
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
    }



    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
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

    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
};

#endif
