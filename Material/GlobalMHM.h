#ifndef GLOBALMHMH
#define GLOBALMHMH

#include <iostream>
#include <cstdarg>

#include "pzmaterial.h"
#include "pzmaterialdata.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include <pzdiscgal.h>

class TPZMaterialData;

class GlobalMHM : public TPZDiscontinuousGalerkin
{
public:

    GlobalMHM() : TPZDiscontinuousGalerkin() { }

	GlobalMHM( int i, int d ) : TPZDiscontinuousGalerkin( i ), dim(d)
    {
        // DO nothing else!
    }

    virtual ~GlobalMHM() { }

    virtual std::string Name() { return "GlobalMHM"; }

    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);

    virtual void Contribute(TPZMaterialData &data,
                          REAL weight,
                          TPZFMatrix &ef)
    {
        // DO nothing!
    }

    virtual void ContributeBC(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ef,
                            TPZBndCond &bc)
    {
        // DO nothing!
    }

    virtual void ContributeInterface( TPZMaterialData & data, TPZMaterialData & dataleft, TPZMaterialData & dataright, REAL weight, TPZFMatrix & ek,TPZFMatrix & ef )
    {
        // DO nothing!
    }

    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc)
    {
        // DO nothing!
    }

    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef)
    { 
        // DO nothing!
    }

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc)
    {
        // DO nothing!
    }

    virtual int Dimension() {return dim;}

    virtual int NStateVariables() {return 1;}

private:
    int dim;
};

#endif