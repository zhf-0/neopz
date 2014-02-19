//---------------------------------------------------------------------------

/**
 * @file
 * @brief Contains the TPZMatLaplacian class.
 */
#ifndef TPZMatDualHybridPoissonH
#define TPZMatDualHybridPoissonH


#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief \f$ -Laplac(u) = fXf  \f$
 */
/**
 * \f$ -Laplac(u) = fXf  \f$
 */
class TPZMatDualHybridPoisson : public TPZDiscontinuousGalerkin {

	protected :

	/** @brief Forcing function value */
	STATE fXf;

  STATE fBetaZero;

public:

	TPZMatDualHybridPoisson(int nummat, REAL f, REAL betaZero);

  TPZMatDualHybridPoisson(int matid);

	TPZMatDualHybridPoisson();

	TPZMatDualHybridPoisson(const TPZMatDualHybridPoisson &copy);

	virtual ~TPZMatDualHybridPoisson();

	virtual TPZMaterial * NewMaterial(){
		return new TPZMatDualHybridPoisson(*this);
	}

	int Dimension() { return 2;}

	int NStateVariables(){
    return 1;
  }

	virtual void Print(std::ostream & out);

	virtual std::string Name() { return "TPZMatDualHybridPoisson"; }

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */

  virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight,
									 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	virtual int VariableIndex(const std::string &name);

	virtual int NSolutionVariables(int var);

	virtual int NFluxes(){ return 2;}

  public:

	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);

	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux)
  {
    DebugStop();
  }

	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);

	virtual int NEvalErrors() {return 3;}

  virtual int ClassId() const {
    DebugStop();
    return -1;
  }

	virtual void Write(TPZStream &buf, int withclassid){
    DebugStop();
  }

	virtual void Read(TPZStream &buf, void *context){
    DebugStop();
  }

};

#endif

