/**
 * @file TPZMatHelmholtz2DHDivRot.h
 * @brief Header file for class TPZMatHelmholtz2DHDivRot.\n
 */

#ifndef TPZMATHELMHOLTZ2DHDIV_H
#define TPZMATHELMHOLTZ2DHDIV_H

#include "TPZMatHelmholtz2D.h"
/**
 * @ingroup material
 * @brief This class implements the inhomogeneous Helmholtz wave equation in 2D.
 */
class TPZMatHelmholtz2DHDivRot : public TPZMatHelmholtz2D {

  public:
    TPZMatHelmholtz2DHDivRot(int id, STATE (&cFunc)(const TPZVec<REAL> &));

    TPZMatHelmholtz2DHDivRot(int id);

    /** @brief Default constructor */
    TPZMatHelmholtz2DHDivRot();

    /** @brief Creates a material object based on the referred object and
     * inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within
     * the vector */
    TPZMatHelmholtz2DHDivRot(const TPZMatHelmholtz2DHDivRot &mat);
    /** @brief Default destructor */
    virtual ~TPZMatHelmholtz2DHDivRot();

    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatHelmholtz2DHDivRot"; }
  public:
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector
     * at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight,
                            TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector
     * at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight,
                              TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                              TPZBndCond &bc);

    /** @brief Returns the solution associated with the var index based on the
     * finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var,
                          TPZVec<STATE> &Solout);
protected:
    void ComputeCurl(TPZFMatrix<REAL> gradScalarPhi , TPZFMatrix<REAL> ivecHCurl , TPZFMatrix<REAL> &curlPhi );
    void RotateForHCurl(TPZVec<REAL> normal , TPZFMatrix<REAL> vHdiv , TPZFMatrix<REAL> &vHcurl );
};

#endif
