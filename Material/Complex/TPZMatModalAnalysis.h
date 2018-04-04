/**
 * @file TPZMatModalAnalysis.h
 * @brief Header file for class TPZMatModalAnalysis.\n
 */

#ifndef TPZMATMODALANALYSIS_H
#define TPZMATMODALANALYSIS_H

#include <TPZVecL2.h>

const REAL M_C  (299792458); //velocidade da luz no vacuo
const REAL M_UZERO  (1.256637061435917e-6);//permeabilidade do meio livre
const REAL M_EZERO  (8.85418781762039e-12);//permissividade do meio livre
#ifdef STATE_COMPLEX
const STATE imaginary(0.,1.);//unidade imaginaria
#endif

/**
 * @ingroup material
 * @brief This class implements the weak statement of a waveguide problem as stated in Jin's 
 * The Finite Element Method in Electromagnetics (chapter 8 of 3rd edition).
 * It used a 2D Hcurl space for the transversal components of the electric field and an 1D
 * H1 space for the longitudinal component.
 */
class  TPZMatModalAnalysis : public TPZVecL2
{
    
protected:
    //COM CERTEZA
    const STATE fUr;
    const STATE fEr;
    const REAL fScaleFactor;
    REAL fLambda;
    const int h1meshindex = 1;
    const int hcurlmeshindex = 0;
    bool isTesting;
    
public:
    
    TPZMatModalAnalysis(int id, REAL lambda, const STATE &ur, const STATE &er, const REAL &scale = 1.);
    
    TPZMatModalAnalysis(int id);
    
    /** @brief Default constructor */
    TPZMatModalAnalysis();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatModalAnalysis(const TPZMatModalAnalysis &mat);
    /** @brief Default destructor */
    ~TPZMatModalAnalysis() override;
    
    /** @brief Returns the name of the material */
    std::string Name() override { return "TPZMatModalAnalysis"; }
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const override {return 2;}
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() override { return 1;}
    
    int HCurlIndex() const { return hcurlmeshindex;}
    int H1Index() const { return h1meshindex;}
    
public:
#ifdef PZDEBUG
    void ContributeValidateFunctions(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
#endif
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
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
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the residual vector
     * @since April 16, 2007
     */
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
    /**
     * @brief This method defines which parameters need to be initialized in order to compute the contribution of an element
     * @param datavec [out] vector of TPZMaterialData, each position will specifie the requirements for its correspondent state variable
     */
    void FillDataRequirements(TPZMaterialData &data) override
    {
        data.SetAllRequirements(false);
    }
    
    
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data) override
    {
        data.SetAllRequirements(false);
        data.fNeedsNormal = true;
    }
    
    /**
     * @brief This method defines which parameters need to be initialized in order to compute the contribution of an element
     * @param datavec [out] vector of TPZMaterialData, each position will specifie the requirements for its correspondent state variable
     */
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override
    {
        int nref = datavec.size();
        for(int iref = 0; iref<nref; iref++){
            datavec[iref].SetAllRequirements(false);
        }
    }
    
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override
    {
        // default is no specific data requirements
        int nref = datavec.size();
        for(int iref = 0; iref<nref; iref++){
            datavec[iref].SetAllRequirements(false);
            datavec[iref].fNeedsNormal = true;
        }
    }
    
    int VariableIndex(const std::string &name) override;
    
    /**
     * @brief Returns the number of variables associated with the variable indexed by var.
     * @param var Index variable into the solution, is obtained by calling VariableIndex
     */
    int NSolutionVariables(int var) override;
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
    
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
};

#endif

