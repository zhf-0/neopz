/**
 * @file
 * @brief Contains the TPZConsLawTest class for test. Material as conservation law
 */
#ifndef CONSLAWTESTHPP
#define CONSLAWTESTHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzconslaw.h"


/**
 * @brief Only to test a material as conservation law. It was used for testing purposes
 * @ingroup material
 */
class TPZConsLawTest  : public TPZConservationLaw {
	
	TPZFMatrix fXf;//fonte
	TPZVec<REAL> fB;
	int fArtificialDiffusion;

	/// Integer for integration degree of the initial solution
	int fTest;
	
	REAL fDelta;

	public :
	
	TPZConsLawTest(int nummat, TPZVec<REAL> B,int artdiff,REAL delta_t,int dim,REAL delta,int test=0);
	
	virtual ~TPZConsLawTest();
	
	void SetMaterial(TPZFMatrix &xfin) {
		fXf = xfin;
	}
	
	REAL DeltaOtimo();
	
	REAL CFL(int degree);
	
	REAL B(int i,TPZVec<REAL> &x);
		
	REAL Delta();
	
	REAL T(int jn,TPZVec<REAL> &x);
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZConsLawTest"; }
	
	virtual void Contribute(TPZMaterialData &data,REAL weight,
                            TPZFMatrix &ek,TPZFMatrix &ef);
	virtual void Contribute(TPZMaterialData &data,REAL weight,
							TPZFMatrix &ef) {
		TPZConservationLaw::Contribute(data,weight,ef);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ek,
							  TPZFMatrix &ef,
							  TPZBndCond &bc);
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZConservationLaw::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 1;}
	
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<REAL> &Solout)
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
	
	/** @brief Compute the value of the flux function to be used by ZZ error estimator */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
	
	void ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright);
	void ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft);

};

#endif
