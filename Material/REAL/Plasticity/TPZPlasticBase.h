/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZPlasticBase.h
 * Author: quinelato
 *
 * Created on July 19, 2017, 4:29 PM
 */

#ifndef TPZPLASTICBASE_H
#define TPZPLASTICBASE_H

#include "TPZTensor.h"
#include "TPZElasticResponse.h"



class TPZPlasticBase {
public:

    virtual ~TPZPlasticBase() {
    };
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) = 0;
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma) = 0;
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) = 0;
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) = 0;
    virtual void SetState(const TPZPlasticState<REAL> &state) = 0;
    virtual TPZPlasticState<REAL> GetState() const = 0;
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const = 0;

    virtual int IntegrationSteps() const {
        return 1;
    }
    virtual void SetElasticResponse(TPZElasticResponse &ER) = 0;
    virtual TPZElasticResponse GetElasticResponse() const = 0;
    //	virtual void SetIntegrTol(REAL integrTol)=0;
    virtual const char * Name()const = 0;
    virtual void Print(std::ostream & out)const = 0;
    //virtual void Write(TPZStream &buf) const = 0;
    virtual void Read(TPZStream &buf) = 0;

};


#endif /* TPZPLASTICBASE_H */

