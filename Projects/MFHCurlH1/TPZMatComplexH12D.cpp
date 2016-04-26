#include "TPZMatComplexH12D.h"
#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif

TPZMatComplexH12D::TPZMatComplexH12D(int id, REAL lambda, REAL kz , STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &), REAL e0, REAL t, REAL scale) : TPZMaterial(id), fUr(ur), fEr(er), fLambda(lambda), fKz (kz) , fE0 (e0) , fTheta(t) , fScale(scale)
{
	fW=2.*M_PI*M_C/fLambda;
}

TPZMatComplexH12D::TPZMatComplexH12D(int id) : TPZMaterial(id), fUr(urDefault),
fEr(erDefault), fLambda(1.55e-9)
{
	fW=2.*M_PI*M_C/fLambda;
	fTheta = 0.;
  fKz = -999;
  fE0 = 1;
  fScale = 1;
}

/** @brief Default constructor */
TPZMatComplexH12D::TPZMatComplexH12D() : TPZMaterial(), fUr(urDefault),
fEr(erDefault), fLambda(1.55e-9)
{
	fW=2.*M_PI*M_C/fLambda;
	fTheta = 0.;
  fKz = -999;
  fE0 = 1;
  fScale = 1;
}


TPZMatComplexH12D::TPZMatComplexH12D(const TPZMatComplexH12D &mat) : TPZMaterial(mat), fUr(mat.fUr),
fEr(mat.fEr)
{
	fLambda = mat.fLambda;
	fTheta = mat.fTheta;
	fE0 = mat.fE0;
	fScale = mat.fScale;
  fKz = mat.fKz;
  fW=2.*M_PI*M_C/fLambda;
}

TPZMatComplexH12D::~TPZMatComplexH12D()
{
	
}

void TPZMatComplexH12D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatComplexH12D::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatComplexH12D::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatComplexH12D::ContributeForcingRTBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
  return;
}

void TPZMatComplexH12D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

void TPZMatComplexH12D::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	DebugStop();
}

void TPZMatComplexH12D::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatComplexH12D::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	DebugStop();
}

int TPZMatComplexH12D::VariableIndex(const std::string &name)
{
	if(name == "absE") {
		return 2;
	}
	else if ( name == "realE")
	{
		return 3;
	}
	else if ( name == "solAnal")
	{
		return 4;
	}
	return TPZMaterial::VariableIndex(name);
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatComplexH12D::NSolutionVariables(int var)
{
	int nVar = 0;
	switch (var) {
		case 2://absE
			nVar = 3;
			break;
		case 3://realE
			nVar = 3;
			break;
		case 4://realE
			nVar = 3;
			break;
		default:
			nVar = TPZMaterial::NSolutionVariables(var);
			break;
	}
	return nVar;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatComplexH12D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
	
	TPZManVector<STATE,3> ax1(3),ax2(3), normal(3);
	for (int i=0; i<3; i++) {
		ax1[i] = data.axes(0,i);
		ax2[i] = data.axes(1,i);
	}
	//ROTATE FOR HCURL
	Cross(ax1, ax2, normal);
	
	Solout.Resize(3);
	Cross(normal, data.sol[0], Solout);
	switch (var) {
		case 2://absE
#ifdef STATE_COMPLEX
			Solout[0] = std::abs(Solout[0]);
			Solout[1] = std::abs(Solout[1]);
			Solout[2] = std::abs(Solout[2]);
#endif
			break;
		case 3://realE
#ifdef STATE_COMPLEX
			Solout[0] = std::real(Solout[0]);
			Solout[1] = std::real(Solout[1]);
			Solout[2] = std::real(Solout[2]);
#endif
			break;
		case 4://solAnal
#ifdef STATE_COMPLEX
			if ( HasfForcingFunctionExact() ) {
				fForcingFunctionExact->Execute(data.x,Solout);
				Solout[0] = std::abs(Solout[0]);
				Solout[1] = std::abs(Solout[1]);
				Solout[2] = std::abs(Solout[2]);
			}
#endif
			break;
		default:
			break;
	}
}









