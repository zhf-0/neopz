/**
 * \file
 * @brief Contains the TPZLinearConvection class which implements a linear scalar convection equation.
 */
/* Generated by Together */
class TPZFMatrix;
class TPZBndCond;

#ifndef TPZLINEARCONVECTION_H
#define TPZLINEARCONVECTION_H
#include "pzmaterial.h"

/**
 * @ingroup material
 * @brief Implements a linear scalar convection equation with modified SUPG difusion
 */
class TPZLinearConvection : public TPZMaterial {
public:  
    TPZLinearConvection(TPZLinearConvection & copy);
	
    TPZLinearConvection(int id,TPZVec<REAL> &conv) ;
	
	
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() ;
	
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables()  ;
	
    /** @brief Returns the number of components which form the flux function */
    virtual int NFluxes() {return 2;}
	
    virtual void Contribute(TPZMaterialData &data, REAL weight,
							TPZFMatrix &ek,TPZFMatrix &ef);
	
	
    virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
	
    virtual void Contribute(TPZMaterialData &data, REAL weight,
							TPZFMatrix &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix &ef,TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }
	
	/** @brief Print out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout);
	
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);
	
    virtual int NSolutionVariables(int var);
	
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<REAL> &Solout)
	{
        Solution(data.sol,data.dsol,data.axes,var,Solout);
    }
	
    virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux) {}
	
    /** @brief To create another material of the same type */
    virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
    /** @brief Reads data of the material from a istream (file data) */
    virtual void SetData(std::istream &data);
	
private:    
    REAL fConvect[2];

};

#endif //TPZLINEARCONVECTION_H
