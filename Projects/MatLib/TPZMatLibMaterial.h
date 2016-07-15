//
//  TPZMatLibMaterial.hpp
//  PZ
//
//  Created by Philippe Devloo on 7/11/16.
//
//

#ifndef TPZMatLibMaterial_hpp
#define TPZMatLibMaterial_hpp

#include <stdio.h>

#ifdef USING_MATLIB
#include "MatLib.h"

#include "pzmatwithmem.h"

struct TPZMatLibMemory
{
    
    zorglib::MaterialState fPrevState, fNextState;
    zorglib::MatLibMatrix fK;
    
    TPZMatLibMemory() : fPrevState(), fNextState(), fK()
    {
        
    }
    
    TPZMatLibMemory(const TPZMatLibMemory &copy) : fPrevState(copy.fPrevState), fNextState(copy.fNextState), fK(copy.fK)
    {
        
    }
    
    TPZMatLibMemory &operator=(const TPZMatLibMemory &copy)
    {
        fPrevState = copy.fPrevState;
        fNextState = copy.fNextState;
        fK = copy.fK;
        return *this;
    }
    
    void Write(TPZStream &buf, int withclassid)
    {
        DebugStop();
    }
    
    void Read(TPZStream &buf, void *context)
    {
        DebugStop();
    }
    
    void Print(std::ostream &out) const
    {
        out << "Previous state deform " << fPrevState.grad << " sig " << fPrevState.flux << std::endl;
        out << "Next state deform " << fNextState.grad << " sig " << fNextState.flux << std::endl;
        out << "Tangent " << fK << std::endl;
    }

    
};

inline std::ostream &operator<<(std::ostream &out, const TPZMatLibMemory &obj)
{
    obj.Print(out);
    return out;
}

class TPZMatLibMaterial : public TPZMatWithMem<TPZMatLibMemory>
{

private:
    int fDimension;
    
    int fZorgLibDimension;
    
    zorglib::MaterialModel fMatModel;
    struct ZorgLibVar
    {
        int fPZIndex;
        int fZorgIndex;
        int fNumVars;
        bool fIsExternal;
        bool fIsGrad;
        std::string fName;
        
        ZorgLibVar() : fPZIndex(-1), fZorgIndex(-1), fNumVars(-1), fIsExternal(false), fIsGrad(false)
        {
            
        }
        ZorgLibVar(int pzindex, int zorgindex, int numvars, bool isExt, bool isGrad, const std::string &name) : fPZIndex(pzindex), fZorgIndex(zorgindex), fNumVars(numvars), fIsExternal(isExt), fIsGrad(isGrad),
            fName(name)
        {
            
        }
        ZorgLibVar(const ZorgLibVar &cp) : fPZIndex(cp.fPZIndex), fZorgIndex(cp.fZorgIndex),
            fNumVars(cp.fNumVars), fIsExternal(cp.fIsExternal), fIsGrad(cp.fIsGrad), fName(cp.fName)
        {
            
        }
        ZorgLibVar &operator=(const ZorgLibVar &cp)
        {
            fPZIndex = cp.fPZIndex;
            fZorgIndex = cp.fZorgIndex;
            fNumVars = cp.fNumVars;
            fIsExternal = cp.fIsExternal;
            fIsGrad = cp.fIsGrad;
            fName = cp.fName;
            return *this;
        }
    };
    
    TPZStack<ZorgLibVar> fZorgLibVars;
    
    void TransferData(ZorgLibVar &a, zorglib::MaterialState &state, TPZVec<STATE> &solution);
    
    static zorglib::ParameterSet fExternalState;
    
public:
    
    TPZMatLibMaterial(int matid, int dimension, const zorglib::MaterialModel &orig) : TPZMatWithMem<TPZMatLibMemory>(matid), fDimension(dimension), fZorgLibDimension(0), fMatModel(orig)
    {
        zorglib::ConstitutiveModel &model = fMatModel.model();
        int numinternalvars = model.nIntVarBundled();
        int count = 10;
        for (int i=0; i<numinternalvars; i++) {
            zorglib::ConstitutiveModel::VariableType vartype = model.typeIntVar(i);
            int numvar = 1;
            ZorgLibVar a(count,i,numvar,false,false,model.labelIntVar(i));
            count++;
            if (vartype == zorglib::ConstitutiveModel::TYPE_SCALAR) {
            }
            else if(vartype == zorglib::ConstitutiveModel::TYPE_VECTOR)
            {
                a.fNumVars = 3;
            }
            else if(vartype == zorglib::ConstitutiveModel::TYPE_SYM_TENSOR)
            {
                a.fNumVars = 9;
            }
            else
            {
                DebugStop();
            }
            fZorgLibVars.Push(a);
        }
        int numextvars = model.nExtVarBundled();
        for (int i=0; i<numextvars; i++) {
            zorglib::ConstitutiveModel::VariableType vartype = model.typeExtVar(i);
            int numvars = 1;
            ZorgLibVar a(count,i,numvars,true,true,model.labelExtVar(i));
            count++;
            ZorgLibVar b(count,i,numvars,true,false,model.labelExtForce(i));
            count++;

            if (vartype == zorglib::ConstitutiveModel::TYPE_SCALAR) {
            }
            else if(vartype == zorglib::ConstitutiveModel::TYPE_VECTOR)
            {
                a.fNumVars = 3;
                b.fNumVars = 3;
            }
            else if(vartype == zorglib::ConstitutiveModel::TYPE_SYM_TENSOR)
            {
                a.fNumVars = 9;
                b.fNumVars = 9;
            }
            else
            {
                DebugStop();
            }
            fZorgLibVars.Push(a);
            fZorgLibVars.Push(b);

        }
        int matmodelnvar = fMatModel.model().nExtVar();
        fZorgLibDimension = 3;
        if (matmodelnvar == 6) {
            fZorgLibDimension = 3;
        }
        else if(matmodelnvar == 4)
        {
            fZorgLibDimension = 2;
        }
        else if(matmodelnvar == 3)
        {
            fZorgLibDimension = 1;
        }

    }
    
    TPZMatLibMaterial(const TPZMatLibMaterial &copy) : TPZMatWithMem<TPZMatLibMemory>(copy), fDimension(copy.fDimension),
        fZorgLibDimension(copy.fZorgLibDimension),
        fMatModel(copy.fMatModel)
    {
        
    }
    
    TPZMatLibMaterial &operator=(const TPZMatLibMaterial &copy)
    {
        TPZMatWithMem<TPZMatLibMemory>::operator=(copy);
        fDimension = copy.fDimension;
        fZorgLibDimension = copy.fZorgLibDimension;
        fMatModel = copy.fMatModel;
        return *this;
    }
    
    virtual TPZMaterial *NewMaterial()
    {
        return new TPZMatLibMaterial(*this);
    }
    /// Consider the iterative process converged : the nextstate becomes previous state
    // this method will loop over all the integration points and update the plastic state
    void LoadState();
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables()
    {
        return fDimension;
    }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const
    {
        return fDimension;
    }
    
    /// change the dimension of the simulation
    void SetDimension(int dim)
    {
        if (dim <1 || dim > 3) {
            DebugStop();
        }
        fDimension = dim;
    }
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        data.fNeedsSol = true;
    }

    /** To be implemented only in the proper materials. */
    virtual int PushMemItem(int sourceIndex = -1);
    

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
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);
    
    /**
     * @brief Returns the number of variables associated with the variable indexed by var.
     * @param var Index variable into the solution, is obtained by calling VariableIndex
     */
    virtual int NSolutionVariables(int var);
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    

};

#endif
#endif /* TPZMatLibMaterial_hpp */
