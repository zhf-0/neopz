/**
 * @file
 * @brief Contains TPZTransientAnalysis class which implements a simple manner to perform transient simulations.
 */
//$Id: pztransientanalysis.h,v 1.7 2011-04-05 19:32:55 calle Exp $

#ifndef TRANSIENTANALH
#define TRANSIENTANALH

#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>
#include <list>

#include "pztransientmat.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzseqsolver.h"
//#include "checkconv.h"

class TPZCompMesh;
class TPZFMatrix;
class TPZFStructMatrix;

using namespace std;

/** 
 * @note It is associated to a TPZTransientMaterial< TRANSIENTCLASS >
 * @brief Implements a very simple manner to perform transient simulations. \ref analysis "Analysis"
 * @ingroup analysis
 */
/** 
 * It uses an implicit or explicit Euler scheme for time derivative
 */
template<class TRANSIENTCLASS>
class TPZTransientAnalysis : public TPZNonLinearAnalysis {
	
public:
	
	/** @brief Static attribute storing the current time of simulation
	 */
	static double gTime;
	
	TPZAutoPointer<TPZMatrix> fMassMatrix;
	
	/** @brief Method for gTime attribute access
	 */
	double GetgTime(){ return gTime; }
	
	/** @brief Constructor
	 * @param mesh for base class
	 * @param IsLinear for optimizating the process time, linear problems have flux tangent computed only once
	 * @param out for base class
	 */
	TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear = false, std::ostream &out = std::cout);
	
	/** @brief Default destructor
	 */
	~TPZTransientAnalysis();
	
	/**
	 * @brief Assemble flux vector and jacobian matrix
	 */
	virtual void Assemble();
	
	/** 
	 * @brief Executes a Newton's method for the solution of the implicit in time equation 
	 */
	virtual void RunTransient(std::ostream &out = std::cout, bool FromBegining = true, bool linesearch = true);
	
	/** 
	 * @brief Solves a explicit Euler's scheme in time
	 */
	virtual void RunExplicit(std::ostream &out = std::cout, bool FromBegining = true);
	
	/** @brief See base class for comments
	 */  
	virtual void PostProcess(int resolution){ TPZAnalysis::PostProcess(resolution);}
	
	/** @brief See base class for comments
	 */
	virtual void PostProcess(int resolution, int dimension);
	
	/** @brief See base class for comments
	 */
	virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout);
	
	/** 
	 * @brief Defines max number of steps and steady state convergence tolerance.
	 */
	void SetConvergence(int niter, REAL eps, bool ForceAllSteps = true);
	
	/** @brief Defines properties of DX file */
	void SetSaveFrequency(int SaveFrequency, int resolution);
	
	/** @brief Defines to save solution vector with SaveFrequency frequency.
	 *
	 * If not set, no solution is kept in the process.
	 */
	void SetSaveSolution(int SaveFrequency);
	
	/** @brief Access to saved solution. Pair of (solution vec, simulation time)
	 */
	std::list< std::pair<TPZFMatrix, REAL> > & GetSavedSolutions();
	
	/** 
	 * @brief Defines max number of steps and error convergence tolerance for Newton's method.
	 */  
	void SetNewtonConvergence(int niter, REAL eps);
	
	/** @brief Access to time step attribute
	 */
	REAL & TimeStep();
	
	/** @brief Sets problem initial solution
	 */
	void SetInitialSolution(TPZFMatrix & InitialSol);
	
	/** @brief Sets problem initial solution as zero
	 */
	void SetInitialSolutionAsZero();
	
	/** @brief Returns current iteration
	 */
	int GetCurrentIter();
    
protected:
	
	/** @brief Flag indicating whether the problem is linear or not. 
	 * 
	 * Linear problems require the computation and decompostition of tangent
	 * matrix only once.
	 */
	bool fIsLinearProblem;
	
	/** @brief Simulation time step */
	REAL fTimeStep;
	
	/** @brief Current iteration. Variable allowing to restart the simulation. */
	int fCurrentIter;
	
	/** @brief Number of iterations counting from fCurrentIter to fCurrentIter+fNIter */
	int fNIter;
	
	/** @brief Tolerance to consider the problem solution as steady state */
	REAL fSteadyTol;
	
	/** @brief Flag indicating whether all steps must be performed even if tolerance is achieved. */
	bool fForceAllSteps;
	
	/** @brief Frequency which solution must be saved in DX file. */
	int fSaveFrequency;  
	/** @brief Resolution of DX mesh */
	int fDXResolution;
	
	/** @brief Frequency which solution vector must be saved.
	 * 
	 *  Zero (default value) means no solution vector but the current one is saved.
	 */
	int fSaveSolutionVecFrequency;
	
	/** @brief Attribute to store solution vectors during process. Pair of (solution vec, simulation time)
	 *
	 * This attribute is cleaned every time Run method is called
	 */
	std::list< std::pair< TPZFMatrix, REAL> > fSavedSolutionVec;
	
	/** @brief If fSaveSolutionVecFrequency != 0, save current solution vector in fSavedSolutionVec attribute. */
	void SaveCurrentSolutionVec();
	
	/** @brief Max iteration number of Newton's method */
	int fNewtonMaxIter;
	
	/** @brief Tolerance of Newton's method */
	REAL fNewtonTol;
	
	/** @brief Sets all materials in temporal scheme as an implicit Euler */
	void SetImplicit();
	
	/** @brief Sets all materials in temporal scheme as an explicit Euler */
	void SetExplicit();
	
	/** @brief Sets all materials in LastState */
	void SetLastState();
	
	/** @brief Sets all materials in CurrentState */
	void SetCurrentState();
	
	/** @brief Sets all materials to compute the mass matrix - used in the explicit scheme */
	void SetMassMatrix();
	
	/** @brief Sets all materials to compute only the flux contributions - used in the explicit scheme */
	void SetFluxOnly();
	
	/** @brief Sets all materials the time step */
	void SetAllMaterialsDeltaT();
	
	/** @brief Computes linear tangent matrix for linear problems */
	void ComputeLinearTangentMatrix();
	
	/** @brief Computes the mass matrix for the explicit scheme */
	void ComputeMassMatrix();
	
	/** @brief Computes the only the flux contribution for the explicit scheme */
	void ComputeFluxOnly();
	
};

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis< TRANSIENTCLASS >::SetConvergence(int niter, REAL eps, bool ForceAllSteps){
	this->fNIter = niter;
	this->fSteadyTol = eps;
	this->fForceAllSteps = ForceAllSteps;
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis< TRANSIENTCLASS >::SetSaveFrequency(int SaveFrequency, int resolution){
	this->fSaveFrequency = SaveFrequency;
	this->fDXResolution = resolution;
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis< TRANSIENTCLASS >::SetNewtonConvergence(int niter, REAL eps){
	this->fNewtonMaxIter = niter;
	this->fNewtonTol = eps;
}

template<class TRANSIENTCLASS>
inline REAL & TPZTransientAnalysis< TRANSIENTCLASS >::TimeStep(){
	return this->fTimeStep;
}

template<class TRANSIENTCLASS>
inline int TPZTransientAnalysis< TRANSIENTCLASS >::GetCurrentIter(){
	return this->fCurrentIter;
}

template<class TRANSIENTCLASS>
double TPZTransientAnalysis<TRANSIENTCLASS>::gTime = 0.;

template<class TRANSIENTCLASS>
TPZTransientAnalysis<TRANSIENTCLASS>::TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear, std::ostream &out):/*TPZAnalysis*/TPZNonLinearAnalysis(mesh,out), fSavedSolutionVec(){
	this->fTimeStep = 0.;
	this->fCurrentIter = 0;
	this->SetConvergence(0, 0.);
	this->SetNewtonConvergence(0, 0.);
	this->SetInitialSolutionAsZero();
	this->fIsLinearProblem = IsLinear;
	this->SetSaveFrequency(0,0);
	this->fSaveSolutionVecFrequency = 0; 
}

template<class TRANSIENTCLASS>
TPZTransientAnalysis<TRANSIENTCLASS>::~TPZTransientAnalysis(){
	fSavedSolutionVec.clear();
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetInitialSolution(TPZFMatrix & InitialSol){
	const int nrows = this->Mesh()->Solution().Rows();
	const int ncols = this->Mesh()->Solution().Cols();
	if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
		PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
	}
	else{
		this->fSolution = InitialSol;
	}
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetInitialSolutionAsZero(){
	TPZFMatrix & MeshSol = this->Mesh()->Solution();
	this->fSolution.Redim( MeshSol.Rows(), MeshSol.Cols() );
	this->fSolution.Zero();
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::RunTransient(std::ostream &out, bool FromBegining, bool linesearch){
	
	this->SetImplicit();
	
	if (FromBegining){
		this->fCurrentIter = 0;
		this->fSavedSolutionVec.clear();
	}
	
#ifdef _DOCHECKCONV_
	{
		TPZVec<REAL> coefs(1,1.);
		TPZFMatrix cpSol(fSolution);
		TPZFMatrix range(fCompMesh->NEquations(),1,1.);
		this->SetLastState();
		CheckConvergence(*this,cpSol,range,coefs);
		this->SetCurrentState();
		CheckConvergence(*this,cpSol,range,coefs);
	}
#endif
	
	this->LoadSolution(fSolution);
	
	TPZTransientAnalysis::gTime =this->GetCurrentIter() * this->TimeStep();
	this->PostProcess(this->fDXResolution);//plot the initial solution
	
	this->SetAllMaterialsDeltaT();
	
	if (this->fIsLinearProblem)// se o problema é linear, não tem pq usar newton, basta resolver o sistema inear.
	{
		if(!this->fMassMatrix)
		{
			this->ComputeMassMatrix();
		}
		this->ComputeLinearTangentMatrix();  // compute ek = K_difussion +1/(delta t)M_mass and   F_bc + F_load
		TPZFMatrix F= this->fRhs;// F = F_bc + F_load
		for( ; this->fCurrentIter < this->fNIter;)
		{
			this->fCurrentIter++;
			TPZFMatrix Un=this->Solution();
			this->fRhs.Zero();
			this->fMassMatrix->MultAdd(Un,F, this->fRhs,1.,1.,1,1); // ef+=(1/delta t)M Un +F
			this->Solve();
			this->PostProcess(this->fDXResolution);
		}
	}
	else//nonlinear problem
	{	
        TPZFMatrix prevsol, laststate, lastsol;
        for( ; this->fCurrentIter < this->fNIter;){
            
            this->fCurrentIter++;
            TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;
            
            //Computing residual of last state solution
            //     this->fSolution = prevsol;
            this->SetLastState();
            this->Assemble();
            laststate = this->fRhs;
            prevsol = fSolution;
            lastsol = fSolution;
            //Newton's method
            this->SetCurrentState();        
            REAL error = this->fNewtonTol * 2. + 1.;    
            int iter = 0;
            while(error > this->fNewtonTol && iter < this->fNewtonMaxIter) {
                
                fSolution.Redim(0,0);
                this->Assemble();
                this->fRhs += laststate;
                this->Solve();
                
                if (linesearch){
                    TPZFMatrix nextSol;
                    REAL LineSearchTol = 1e-3 * Norm(fSolution);
                    const int niter = 100;
                    this->LineSearch(prevsol, fSolution, nextSol, LineSearchTol, niter);
                    fSolution = nextSol;
                }
                else{
                    fSolution += prevsol;
                }
                
                prevsol -= fSolution;
                REAL norm = Norm(prevsol);
                out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << std::endl;
                
                prevsol = fSolution;
                TPZAnalysis::LoadSolution();
                
                error = norm;
                iter++;
            }//Newton's iterations
            
            prevsol = fSolution;
            
            if (this->fSaveFrequency){
                if (!(this->fCurrentIter % fSaveFrequency)){
                    this->PostProcess(this->fDXResolution);
                }
            }
            this->SaveCurrentSolutionVec();
            
            prevsol -= lastsol;
            REAL steadynorm = Norm(prevsol);
            std::cout << "*********** Steady state error at iteration " << this->fCurrentIter << " = " << steadynorm << "\n\n";
            if (!fForceAllSteps){
                if (steadynorm < this->fSteadyTol){
                    std::cout << "Steady state solution achieved\n\n";
                    this->fNIter = fCurrentIter;
                    break;
                }
            }
            std::cout.flush();   
            
        }//time iterations
	}
    
}//method

//template<class TRANSIENTCLASS>
//void TPZTransientAnalysis<TRANSIENTCLASS>::RunTransient(std::ostream &out, bool FromBegining, bool linesearch){
//	
//	this->SetImplicit();
//	
//	if (FromBegining){
//		this->fCurrentIter = 0;
//		this->fSavedSolutionVec.clear();
//	}
//	
//#ifdef _DOCHECKCONV_
//	{
//		TPZVec<REAL> coefs(1,1.);
//		TPZFMatrix cpSol(fSolution);
//		TPZFMatrix range(fCompMesh->NEquations(),1,1.);
//		this->SetLastState();
//		CheckConvergence(*this,cpSol,range,coefs);
//		this->SetCurrentState();
//		CheckConvergence(*this,cpSol,range,coefs);
//	}
//#endif
//	
//	this->LoadSolution(fSolution);
//	
//	TPZTransientAnalysis::gTime =this->GetCurrentIter() * this->TimeStep();
//	//   this->PostProcess(this->fDXResolution);
//	
//	this->SetAllMaterialsDeltaT();
//	
//	if (this->fIsLinearProblem){
//		this->ComputeLinearTangentMatrix();  
//	}
//	
//	TPZFMatrix prevsol, laststate, lastsol;
//	for( ; this->fCurrentIter < this->fNIter;){
//		
//		this->fCurrentIter++;
//		TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;
//		
//		//Computing residual of last state solution
//		//     this->fSolution = prevsol;
//		this->SetLastState();
//		this->Assemble();
//		laststate = this->fRhs;
//		prevsol = fSolution;
//		lastsol = fSolution;
//		//Newton's method
//		this->SetCurrentState();        
//		REAL error = this->fNewtonTol * 2. + 1.;    
//		int iter = 0;
//		while(error > this->fNewtonTol && iter < this->fNewtonMaxIter) {
//			
//			fSolution.Redim(0,0);
//			this->Assemble();
//			this->fRhs += laststate;
//			this->Solve();
//			
//			if (linesearch){
//				TPZFMatrix nextSol;
//				REAL LineSearchTol = 1e-3 * Norm(fSolution);
//				const int niter = 100;
//				this->LineSearch(prevsol, fSolution, nextSol, LineSearchTol, niter);
//				fSolution = nextSol;
//			}
//			else{
//				fSolution += prevsol;
//			}
//			
//			prevsol -= fSolution;
//			REAL norm = Norm(prevsol);
//			out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << std::endl;
//			
//			prevsol = fSolution;
//			TPZAnalysis::LoadSolution();
//			
//			error = norm;
//			iter++;
//		}//Newton's iterations
//		
//		prevsol = fSolution;
//		
//		if (this->fSaveFrequency){
//			if (!(this->fCurrentIter % fSaveFrequency)){
//				this->PostProcess(this->fDXResolution);
//			}
//		}
//		this->SaveCurrentSolutionVec();
//		
//		prevsol -= lastsol;
//		REAL steadynorm = Norm(prevsol);
//		std::cout << "*********** Steady state error at iteration " << this->fCurrentIter << " = " << steadynorm << "\n\n";
//		if (!fForceAllSteps){
//			if (steadynorm < this->fSteadyTol){
//				std::cout << "Steady state solution achieved\n\n";
//				this->fNIter = fCurrentIter;
//				break;
//			}
//		}
//		std::cout.flush();   
//		
//	}//time iterations
//    
//}//method


template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetImplicit(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
		if (trans){
			trans->SetImplicit();
		}
	}
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetExplicit(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
		if (trans){
			trans->SetExplicit();
		}
	}
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetLastState(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
		if (trans){
			trans->SetLastState();
		}
	}
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetCurrentState(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
		if (trans){
			trans->SetCurrentState();
		}
	}
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetMassMatrix(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
		if (trans){
			trans->SetMassMatrix();
		}
	}
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetFluxOnly(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
		if (trans){
			trans->SetFluxOnly();
		}
	}
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetAllMaterialsDeltaT(){
	TPZCompMesh * mesh = this->Mesh();
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
	for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
		if (trans){
			trans->SetTimeStep(this->TimeStep());
		}
	}
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::PostProcess(int resolution, int dimension){
    REAL T = this->GetCurrentIter() * this->TimeStep();
    this->fTime = T;
    TPZAnalysis::PostProcess(resolution, dimension);
}//method

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::PostProcess(TPZVec<REAL> &loc, std::ostream &out){
    REAL T = this->GetCurrentIter() * this->TimeStep();
    this->gTime = T;
    out << "\nSOLUTION #" << this->GetCurrentIter() << " AT TIME = " << T << std::endl;
    TPZAnalysis::PostProcess(loc, out);
    out << "\n***************************************\n" << std::endl;
}//method

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::Assemble(){
	if(!fCompMesh || !fStructMatrix || !fSolver){
		cout << "TPZTransientAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
		<< " fStructMatrix " << (bool) fStructMatrix << " fSolver " << (bool) fSolver << " at file " 
		<< __FILE__ << " line " << __LINE__ << endl;
		return;
	}
	
	int sz = fCompMesh->NEquations();
	fRhs.Redim(sz,1);
	
	bool exist = false;
	if(fSolver->Matrix()) if (fSolver->Matrix()->Rows()==sz) exist = true;
	TPZAutoPointer<TPZGuiInterface> inter = new TPZGuiInterface;
	if (exist){
		if (fIsLinearProblem){
			//      TPZStructMatrix::Assemble(fRhs, *Mesh());
			fStructMatrix->Assemble(fRhs,inter);
		}
		else{
			fSolver->Matrix()->Zero();
			fStructMatrix->Assemble((TPZMatrix&)fSolver->Matrix(),fRhs,inter);
		}
	}
	else{
		if (this->fIsLinearProblem){
			std::cout << __PRETTY_FUNCTION__ << " @ " << __LINE__ << " Error! StrMatrix must be created using" 
			<< " methodTPZTransientAnalysis::ComputeLinearTangentMatrix()"
			<< " when (this->fIsLinearProblem == true)\n";
		}
		TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs,NULL);
		fSolver->SetMatrix(mat);
	}
	fSolver->UpdateFrom(fSolver->Matrix());
}

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeLinearTangentMatrix(){
	if (!fIsLinearProblem) return;      
	this->SetCurrentState();
	const int sz = this->Mesh()->NEquations();
	fRhs.Redim(sz,1);
	TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs,NULL);
	fSolver->SetMatrix(mat);
}//method

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeMassMatrix(){
	this->SetMassMatrix();
	const int sz = this->Mesh()->NEquations();
	fRhs.Redim(sz,1);
	TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs,NULL);
	fMassMatrix = mat->Clone();
	fSolver->SetMatrix(mat);
}//method

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeFluxOnly(){
	if(!fCompMesh || !fStructMatrix || !fSolver){
		cout << "TPZTransientAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
		<< " fStructMatrix " << (bool) fStructMatrix << " fSolver " << (bool) fSolver << " at file " 
		<< __FILE__ << " line " << __LINE__ << endl;
		return;
	}
	
	this->SetFluxOnly();  
	int sz = fCompMesh->NEquations();
	fRhs.Redim(sz,1);
	if(fSolver->Matrix() && fSolver->Matrix()->Rows()==sz){
		fStructMatrix->Assemble(fRhs,NULL);
		//    TPZStructMatrix::Assemble(fRhs, *Mesh());
	}//if
}//method

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::RunExplicit(std::ostream &out, bool FromBegining){
	
	this->SetExplicit();
	
	if (FromBegining){
		this->fCurrentIter = 0;
		this->fSavedSolutionVec.clear();
	}
	TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;
	this->PostProcess(this->fDXResolution);
	
	this->SetAllMaterialsDeltaT();
	
	TPZFMatrix prevsol;
	for( this->fCurrentIter++ ; this->fCurrentIter < this->fNIter; this->fCurrentIter++){
		
		this->ComputeMassMatrix();
		
		TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;
		
		this->SetFluxOnly();
		
		//Computing residual of last state solution
		prevsol = fSolution;
		TPZAnalysis::LoadSolution();
		this->ComputeFluxOnly();
		
		this->Solve();
		//now fSolution = deltaSol
		fSolution += prevsol;
		
		TPZAnalysis::LoadSolution();
		if (this->fSaveFrequency){
			if (!(this->fCurrentIter % fSaveFrequency)){
				this->PostProcess(this->fDXResolution);
			}
		}
		this->SaveCurrentSolutionVec();
		
		prevsol -= fSolution;
		REAL steadynorm = Norm(prevsol);
		std::cout << "*********** Steady state error at iteration " << (this->fCurrentIter) << " = " << steadynorm << "\n\n";
		if (!fForceAllSteps){
			if (steadynorm < this->fSteadyTol){
				std::cout << "Steady state solution achieved\n\n";
				this->fNIter = fCurrentIter;
				break;
			}
		}
		std::cout.flush();   
		
	}//time iterations
	
}//method

template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SetSaveSolution(int SaveFrequency){
	this->fSaveSolutionVecFrequency = SaveFrequency;
}

template<class TRANSIENTCLASS>
inline std::list< std::pair<TPZFMatrix, REAL> > & TPZTransientAnalysis<TRANSIENTCLASS>::GetSavedSolutions(){
	return this->fSavedSolutionVec;
}

#include <sstream>
template<class TRANSIENTCLASS>
inline void TPZTransientAnalysis<TRANSIENTCLASS>::SaveCurrentSolutionVec(){
	if(!this->fSaveSolutionVecFrequency) return;
	if(this->fCurrentIter % this->fSaveSolutionVecFrequency == 0){
		std::pair< TPZFMatrix, REAL > mypair;
		mypair.first = this->Solution();
		mypair.second = TPZTransientAnalysis::gTime;
		this->fSavedSolutionVec.push_back(mypair);
		
		ofstream file("currentsol.txt");
		stringstream mess; mess << "sol( " << TPZTransientAnalysis::gTime << " ) = ";
		this->Solution().Print(mess.str().c_str(), file);
		
	}
}


#endif
