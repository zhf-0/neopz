/**
 * \file
 * @brief Contains the TPZTransientMaterial class which implements an implicit Euler time integrator.
 */

//$Id: pztransientmat.h,v 1.7 2009-05-06 20:13:37 fortiago Exp $


#ifndef TRANSIENTMATH
#define TRANSIENTMATH

#include "pzmaterial.h"
// #include "pzvec.h"
// #include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief Implements an implicit Euler time integrator. \ref material Material 
 */
/**
 * Weak statement is supposed to be Integral[(un+1 - un)/deltaT * v, Omega] + Bilinear Form = Linear Form
 * This class implements only Integral[(un+1 - un)/deltaT * v, Omega]. Bilinear and linear form must be implemented in base class TBASEMAT.
 */
template<class TBASEMAT>
class TPZTransientMaterial : public TBASEMAT {
	
public:
	
	/** @brief Class constructor */
	TPZTransientMaterial(int nummat, int dim, REAL TimeStep);

    /** Quite default quite copy constructor.. Inserted by Frederico (LNCC) **/
	TPZTransientMaterial( TBASEMAT& cp, REAL TimeStep);

	/** @brief Default destructor */
	~TPZTransientMaterial();
    
	/** @brief Copy constructor */
	TPZTransientMaterial(const TPZTransientMaterial &cp);
	
	/** @brief Sets integral scheme as an explicit Euler */
	void SetExplicit();
	
	/** √Sets integral scheme as an implicit Euler */
	void SetImplicit();
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ek,
							  TPZFMatrix &ef,
							  TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                       REAL weight,
                                       TPZFMatrix &ek,
                                       TPZFMatrix &ef,
                                       TPZBndCond &bc);
	
	/** @brief Set material to compute only Integral[- un/deltaT * v, Omega] */
	void SetLastState();
	
	/** @brief Set material to compute Integral[un+1/deltaT * v, Omega] + Bilinear Form = Linear Form  */
	void SetCurrentState();
	
	/** @brief Set material to compute ek = Integral[phi_i phi_j, Omega]/deltaT */
	void SetMassMatrix();
	
	/** @brief Set material to compute ef = Linear Form - Bilinear Form(u) = F -ku */ 
	void SetFluxOnly();
	
	/** @brief Define time step DeltaT */
	void SetTimeStep(REAL TimeStep);
	
	/** @brief Returns time step value. */
	REAL TimeStep();
	
	/** @brief Indicates if the material requires the solution to compute Contribute */
	/** 
	 * By default its value is true, but it can be set as false by derived material classes \n
	 * to increase the performance of method TPZCompEl::CalcStiff
	 */
	virtual bool NeedsSolutionToContribute(){
		return true;
	}
	
	/** @brief Indicates if the material requires the global coordinate X to compute Contribute */
	/**
	 * By default its value is true, but it can be set as false by derived material classes \n
	 * to increase the performance of method TPZCompEl::CalcStiff
	 */
	virtual bool NeedsXCoord(){
		return (this->fStep != ELast);
	}
	
protected:
	
	enum ETemporalScheme{EImplicit = 1, EExplicit = 2};
	
	ETemporalScheme fTemporalIntegrator;
	
	enum STEPS{ENone = -1, ELast = 0, ECurrent = 1, EMassMatrix = 2, EFluxOnly = 3};
	
	STEPS fStep;
	
	REAL fTimeStep;
	
	virtual void ContributeSolutionRhs(TPZVec<REAL> &sol, TPZFMatrix &phi, REAL weight, TPZFMatrix &ef);
	
	virtual void ContributeTangent(TPZVec<REAL> &sol, TPZFMatrix &phi, REAL weight, TPZFMatrix &ek);
};

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetLastState(){
	this->fStep = ELast;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetCurrentState(){
	this->fStep = ECurrent;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetMassMatrix(){
	this->fStep = EMassMatrix;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetFluxOnly(){
	this->fStep = EFluxOnly;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetTimeStep(REAL TimeStep){
	this->fTimeStep = TimeStep;
}

template<class TBASEMAT>
inline REAL TPZTransientMaterial< TBASEMAT >::TimeStep(){
	return this->fTimeStep;
}

/**
 * \file
 * @brief Contains implementations of the TPZTransientMaterial methods.
 */

//$Id: pztransientmat.cpp,v 1.8 2009-05-06 20:22:18 fortiago Exp $

//#include "pztransientmat.h"

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetExplicit(){
	this->fTemporalIntegrator = EExplicit;
}
template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetImplicit(){
	this->fTemporalIntegrator = EImplicit;
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::TPZTransientMaterial(int nummat, int dim, REAL TimeStep):TBASEMAT(nummat, dim){
	this->SetTimeStep(TimeStep);
}

// Inserted by Frederico (LNCC) for TPZMatComposite.
template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::TPZTransientMaterial( TBASEMAT& cp, REAL TimeStep):TBASEMAT(cp){
	this->SetTimeStep(TimeStep);
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::TPZTransientMaterial(const TPZTransientMaterial &cp):TBASEMAT(cp){
	this->fTemporalIntegrator = cp.fTemporalIntegrator;
	this->fStep = cp.fStep;
	this->fTimeStep = cp.fTimeStep;
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::~TPZTransientMaterial(){
	//NOTHING TO BE DONE
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::Contribute(TPZMaterialData &data,
                                                  REAL weight,
                                                  TPZFMatrix &ek,
                                                  TPZFMatrix &ef){
	
	// Mostly for implicit
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
    
	if (this->fStep == ECurrent){
		TBASEMAT::Contribute(data,weight,ek,ef);
		//this->ContributeSolutionRhs(data.sol[0], data.phi, weight, ef); // joao para o projeto lncc2 => comenta essa linha
		this->ContributeTangent(data.sol[0], data.phi, weight, ek);
		return;
	}
	
	if (this->fStep == ELast){
		this->ContributeSolutionRhs(data.sol[0], data.phi, weight, ef);
		return;
	}
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		this->ContributeTangent(data.sol[0], data.phi, weight, ek);
		return;
	}
	
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TBASEMAT::Contribute(data,weight,ek,ef);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::ContributeBC(TPZMaterialData &data,
                                                    REAL weight,
                                                    TPZFMatrix &ek,
                                                    TPZFMatrix &ef,
                                                    TPZBndCond &bc){
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeBC(data,weight,ek,ef,bc);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		TPZFNMatrix<1000> fakeef(ek.Rows(),1,0.);
		TBASEMAT::ContributeBC(data,weight,ek,fakeef,bc);
		return;
	}
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TPZFNMatrix<1000> fakeef(ef.Rows(),ef.Rows(),0.);
		TBASEMAT::ContributeBC(data,weight,ek,ef,bc);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                                           REAL weight,
                                                           TPZFMatrix &ek,
                                                           TPZFMatrix &ef){
	
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeInterface(data,dataleft,dataright, weight, ek, ef);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		return;
	}
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TBASEMAT::ContributeInterface(data,dataleft,dataright, weight, ek, ef);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                                             REAL weight, 
                                                             TPZFMatrix &ek,
                                                             TPZFMatrix &ef,
                                                             TPZBndCond &bc){
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeBCInterface(data,dataleft, weight,ek, ef, bc);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		return;
	}
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TBASEMAT::ContributeBCInterface(data, dataleft, weight,  ek, ef, bc);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::ContributeSolutionRhs(TPZVec<REAL> &sol, TPZFMatrix &phi, REAL weight, TPZFMatrix &ef){
	REAL Mult = +1.; 
	//Last solution is added to residual
	if (this->fStep == ECurrent) Mult = -1.; 
	//Current solution is subtracted from residual
	const int phr = phi.Rows();
	const int nstate = this->NStateVariables();
	const REAL DeltaT = this->TimeStep();
	int i, k;
	for(i = 0; i < phr; i++) {
		for(k = 0; k < nstate; k++){
			ef(i*nstate+k, 0) += Mult * weight * sol[k] * phi(i,0) / DeltaT;
		}//k
	}//i
}//method

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::ContributeTangent(TPZVec<REAL> &sol, TPZFMatrix &phi, REAL weight, TPZFMatrix &ek){
	const int phr = phi.Rows();
	const int nstate = this->NStateVariables();
	const REAL DeltaT = this->TimeStep();
	int i, j, k;
	for(i = 0; i < phr; i++) {
		for(j = 0; j < phr; j++){
			for(k = 0; k < nstate; k++){
				ek(i*nstate+k, j*nstate+k) += weight * phi(i,0) * phi(j,0) / DeltaT;
			}//k
		}//j
	}//i
}//method


#endif
