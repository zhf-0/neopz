/**
 * @file TPZMatMFHDivRotH1.h
 * @brief Header file for class TPZMatMFHDivRotH1.\n
 */

#ifndef TPZMATMFHDIVROTH1_H
#define TPZMATMFHDIVROTH1_H

#include "TPZMatModalAnalysis.h"
#include "TPZMatHCurlProjection.h"

/**
 * @ingroup material
 * @brief This class implements the weak statement of a waveguide problem as stated in Jin's 
 * The Finite Element Method in Electromagnetics (chapter 8 of 3rd edition).
 * It used a 2D Hcurl space for the transversal components of the electric field and an 1D
 * H1 space for the longitudinal component.
 */
class  TPZMatMFHDivRotH1 : public TPZMatModalAnalysis
{
    
public:
    
    TPZMatMFHDivRotH1(int id, REAL freq, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &));
    
    TPZMatMFHDivRotH1(int id);
    
    /** @brief Default constructor */
    TPZMatMFHDivRotH1();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatMFHDivRotH1(const TPZMatMFHDivRotH1 &mat);
    /** @brief Default destructor */
    virtual ~TPZMatMFHDivRotH1();
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatMFHDivRotH1"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const {return 2;}
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() { return 1;}
    
    int HCurlIndex() const { return hcurlmeshindex;}
    int H1Index() const { return h1meshindex;}
    
public:
    virtual void ContributeValidateFunctions(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
   
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
    
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
};

#endif

