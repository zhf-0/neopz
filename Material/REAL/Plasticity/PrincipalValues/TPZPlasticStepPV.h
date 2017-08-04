/**
 * @file
 */

#ifndef TPZPlasticStepPV_H
#define TPZPlasticStepPV_H


#include "TPZTensor.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "TPZPlasticState.h"
#include "TPZPlasticIntegrMem.h"
#include "pzlog.h"
#include "pzstepsolver.h"
#include "TPZElasticResponse.h"
#include "TPZPlasticBase.h"

#include <set>
#include <ostream>

// Metodos para deixar o programa mais "encapsulado"
TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat);

/*
 
 enum EElastoPlastic
 {
 EAuto = 0,
 EForceElastic = 1,
 EForcePlastic = 2
 };
 */

/**
 * @brief Classe que efetua avanco de um passo de plastificacao utilizando o metodo de Newton
 */
template <class YieldCriterion_type, class ElasticResponse_type = TPZElasticResponse>
class TPZPlasticStepPV : public TPZPlasticBase
{
public:

    /**
     * @brief Constructor which Initialize the plastic material damage variable only
     *
     * @param[in] alpha damage variable
     */

  TPZPlasticStepPV(REAL alpha=0.):fYieldCriterion(), fElasticResponse(), fResidualTol(1.e-12), fMaxNewtonIterations(30), fPlasticState()
	{ 
        fPlasticState.fAlpha = alpha;
    }

    /**
     * @brief Copy Constructor
     *
     * @param[in] source of copy
     */
	TPZPlasticStepPV(const TPZPlasticStepPV & source)
	{
        fYieldCriterion = source.fYieldCriterion;
        fElasticResponse = source.fElasticResponse;
        fResidualTol = source.fResidualTol;
        fMaxNewtonIterations = source.fMaxNewtonIterations;
        fPlasticState = source.fPlasticState;
    }

    /**
     * @brief Operator =
     *
     * @param[in] source of copy
     */
	TPZPlasticStepPV & operator=(const TPZPlasticStepPV & source)
	{
        fYieldCriterion = source.fYieldCriterion;
        fElasticResponse = source.fElasticResponse;
        fResidualTol = source.fResidualTol;
        fMaxNewtonIterations = source.fMaxNewtonIterations;
        fPlasticState = source.fPlasticState;

        return *this;
    }

    /**
     * @brief Name of the class ina string
     */
	virtual const char * Name() const
	{
        return "TPZPlasticStepPV";
    }

	virtual void Print(std::ostream & out) const
	{
        out << "\n" << this->Name();
        out << "\n YieldCriterion:";
        //fYC.Print(out); FAZER O PRINT
        out << "\n ElasticResponse:";
        fElasticResponse.Print(out);
        out << "\nTPZPlasticStepPV Internal members:";
        out << "\n fResidualTol = " << fResidualTol;
        out << "\n fMaxNewtonIterations = " << fMaxNewtonIterations;
        out << "\n fPlasticState = "; // PlasticState
        fPlasticState.Print(out);
    }

    typedef YieldCriterion_type fNYields;

    /**
     * Imposes the specified strain tensor, evaluating the plastic integration if necessary.
     *
     * @param[in] epsTotal Imposed total strain tensor
     */
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal);

    /**
     * Imposes the specified strain tensor and returns the correspondent stress state.
     *
     * @param[in] epsTotal Imposed total strain tensor
     * @param[out] sigma Resultant stress
     */
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma);

    /**
     * Imposes the specified strain tensor and returns the correspondent
     * stress state and tangent stiffness matrix.
     *
     * @param[in] epsTotal Imposed total strain tensor
     * @param[out] sigma Resultant stress
     * @param[out] Dep Incremental constitutive relation
     */
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep);

    /**
     * Attempts to compute an epsTotal value in order to reach an imposed stress state sigma.
     * This method should be used only for test purposes because it isn't fully robust. Some
     * materials, specially those perfectly plastic and with softening, may fail when applying
     * the Newton Method on ProcessLoad.
     *
     * @param[in] sigma stress tensor
     * @param[out] epsTotal deformation tensor
     */
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal);

    virtual TPZPlasticState<REAL> GetState() const;
    /**
     * @brief Return the value of the yield functions for the given strain
     * @param[in] epsTotal strain tensor (total strain)
     * @param[out] phi vector of yield functions
     */
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const;

    virtual void SetElasticResponse(TPZElasticResponse &ER);

    virtual TPZElasticResponse GetElasticResponse() const
    {
        return fElasticResponse;
    }

    /**
     * @brief Updates the damage values
     * @param[in] state New plastic state
     */
    virtual void SetState(const TPZPlasticState<REAL> &state);


    //void CopyFromFNMatrixToTensor(TPZFNMatrix<6> FNM,TPZTensor<STATE> &copy);
    void CopyFromTensorToFMatrix(TPZTensor<STATE> tensor, TPZFMatrix<STATE> &copy) const;


    //void CopyFromTensorToFNMatrix(TPZTensor<STATE> tensor,TPZFNMatrix<6> &copy);
    void CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM, TPZTensor<STATE> &copy) const;




    virtual void Read(TPZStream &buf);

    virtual void Write(TPZStream &buf) const;


    /**
     * Does the TaylorCheck of the tangent matrix
     *
     * @param[in] epsTotal Imposed total strain tensor
     * @param[out] sigma Resultant stress
     * @param[out] Dep Incremental constitutive relation
     */
    void TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv);


    REAL ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat);

	void SetResidualTolerance(STATE tol)
	{
        fResidualTol = tol;
    }

    void ResetPlasticMem()
    {
        //fPlasticMem.Resize(0);
    }

    //virtual void Write(TPZStream &buf) const;

    //virtual void Read(TPZStream &buf);

    /** @brief Object which represents the yield criterion */
    YieldCriterion_type fYieldCriterion;

    /** @brief Object representing the elastic response */
    ElasticResponse_type fElasticResponse;

    /** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */
    TPZPlasticState<STATE> fPlasticState;

    int fYield;
protected:

    /** @brief Residual tolerance accepted in the plastic loop processes */
    REAL fResidualTol;

    /** @brief Maximum number of Newton interations allowed in the nonlinear solvers */
    int fMaxNewtonIterations; // COLOCAR = 30 (sugestao do erick!)


void ApplyStrain(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma);

};

#endif //TPZPlasticStepPV_H
