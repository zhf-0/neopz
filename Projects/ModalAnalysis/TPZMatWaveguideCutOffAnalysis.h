/**
 * @file TPZMatWaveguideCutOffAnalysis.h
 * @brief Header file for class TPZMatWaveguideCutOffAnalysis.\n
 */

#ifndef TPZMATWAVEGUIDECUTOFFANALYSIS_H
#define TPZMATWAVEGUIDECUTOFFANALYSIS_H

#include "TPZVecL2.h"
#include "pzaxestools.h"
#include "pzvec_extras.h"
#include "TPZMatModalAnalysis.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement of the model problem from Oden's book, Chapter 1, within the PZ environment
 */
class  TPZMatWaveguideCutOffAnalysis : public TPZMatModalAnalysis
{
    
public:
    
    TPZMatWaveguideCutOffAnalysis(int id, REAL f0, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &));
    
    ~TPZMatWaveguideCutOffAnalysis();
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatWaveguideCutOffAnalysis"; }

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
};

#endif

