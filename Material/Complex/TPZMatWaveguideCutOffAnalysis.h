/**
 * @file TPZMatWaveguideCutOffAnalysis.h
 * @brief Header file for class TPZMatWaveguideCutOffAnalysis.\n
 */

#ifndef TPZMATWAVEGUIDECUTOFFANALYSIS_H
#define TPZMATWAVEGUIDECUTOFFANALYSIS_H

#include "TPZMatModalAnalysis.h"

class  TPZMatWaveguideCutOffAnalysis : public TPZMatModalAnalysis
{
    
public:
    
    TPZMatWaveguideCutOffAnalysis(int id, REAL f0, const STATE &ur,const STATE &er, const REAL &scale);
    
    ~TPZMatWaveguideCutOffAnalysis();
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatWaveguideCutOffAnalysis"; }

    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
};

#endif

