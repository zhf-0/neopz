#ifndef LOCALMHMH
#define LOCALMHMH

#include <iostream>
#include <cstdarg>
#include <vector>

#include "pzmaterial.h"
#include "pzmaterialdata.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include "pzdiscgal.h"
#include "pzcompel.h"
#include <pzinterpolationspace.h>

class TPZMaterialData;

class LocalMHM : public TPZDiscontinuousGalerkin
{
public:

    LocalMHM() : TPZDiscontinuousGalerkin() { }

	LocalMHM( int i, int d );

    virtual ~LocalMHM() { }

    virtual std::string Name() { return "LocalMHM"; }

    // ==================
    // Contribute methods
    // ==================

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

    // ==================
    // Contribute methods (for Multiphysics simulations)
    // ==================

    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);

    // ================
    // Solution methods
    // ================

    virtual int Dimension() {return dim;}

    virtual int NStateVariables() {return 1;}

    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);

    virtual int NSolutionVariables(int var);

    virtual int VariableIndex(const std::string &name);

    // ====================
    // MHM specific methods
    // ====================

    void SetProblemParams(int cf, int dof) { curFace = cf; curDof = dof; }

    void SetFaces(std::vector<TPZCompEl*>& f, std::vector<TPZGeoEl*>& g) { compFaces = f; geoFaces = g; }

    void SetLambdaProblem(bool t) { LambdaProblem = t; }

    void TransformX( const TPZManVector<REAL,3>& dataX, TPZVec<REAL>& curX );

private:

    int dim;
    TPZFMatrix fEpsilon;
    int curFace;
    int curDof;
    std::vector<TPZCompEl*> compFaces;
    std::vector<TPZGeoEl*> geoFaces;

    bool LambdaProblem;

    void bcFunc(TPZMaterialData& , TPZVec<REAL>& res);
    bool IsXOnTheFace( TPZGeoEl* geo, TPZVec<REAL>& X );
};

#endif