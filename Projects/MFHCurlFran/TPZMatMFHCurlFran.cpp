#include "TPZMatMFHCurlFran.h"
#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif

TPZMatMFHCurlFran::TPZMatMFHCurlFran(int id, REAL lambda, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &), REAL e0, REAL t, REAL scale) : TPZVecL2(id), fUr(ur), fEr(er), fLambda(lambda), fE0 (e0) , fTheta(t) , fScale(scale)
{
	fW=2.*M_PI*M_C/fLambda;
}

TPZMatMFHCurlFran::TPZMatMFHCurlFran(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault), fLambda(1.55e-9)
{
	fW=2.*M_PI*M_C/fLambda;
	fTheta = 0.;
}

/** @brief Default constructor */
TPZMatMFHCurlFran::TPZMatMFHCurlFran() : TPZVecL2(), fUr(urDefault),
fEr(erDefault), fLambda(1.55e-9)
{
	fW=2.*M_PI*M_C/fLambda;
	fTheta = 0.;
}


TPZMatMFHCurlFran::TPZMatMFHCurlFran(const TPZMatMFHCurlFran &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr)
{
	fLambda = mat.fLambda;
	fTheta = mat.fTheta;
	fE0 = mat.fE0;
	fScale = mat.fScale;
	fW=2.*M_PI*M_C/fLambda;
}

TPZMatMFHCurlFran::~TPZMatMFHCurlFran()
{
	
}

void TPZMatMFHCurlFran::ComputeCurl(TPZFMatrix<REAL> gradScalarPhi , TPZFMatrix<REAL> ivecHCurl , TPZFMatrix<REAL> &curlPhi ){
  int nFunctions = gradScalarPhi.Rows();
  curlPhi.Resize( nFunctions , 3);
  curlPhi.Zero();
  TPZManVector<REAL,3> result(3,0.);
  
  for (int i = 0 ; i < nFunctions; i++) {
    curlPhi(i,0) = gradScalarPhi(i,1)*ivecHCurl(i,2) - ivecHCurl(i,1)*gradScalarPhi(i,2);
    curlPhi(i,1) = gradScalarPhi(i,2)*ivecHCurl(i,0) - ivecHCurl(i,2)*gradScalarPhi(i,0);
    curlPhi(i,2) = gradScalarPhi(i,0)*ivecHCurl(i,1) - ivecHCurl(i,0)*gradScalarPhi(i,1);
  }
}

void TPZMatMFHCurlFran::ContributeValidateFunctions(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}



void TPZMatMFHCurlFran::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatMFHCurlFran::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  TPZFNMatrix<12,REAL> phiSca = datavec[0].phi;
  TPZFNMatrix<36,REAL> dphiSca = datavec[0].dphix;
  
  TPZFNMatrix<12,REAL> phiQ = datavec[1].phi;
  TPZManVector<REAL,3> x = datavec[0].x;

  int phrq = datavec[1].fVecShapeIndex.NElements();
  std::cout<<"x"<<std::endl<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
  
  /*********************CREATE HCURL FUNCTIONS****************************/
  
  TPZFNMatrix< 36 , REAL > phiVecHCurl(phrq , 3 , 0.);
  for (int iq = 0 ; iq < phrq ; iq++) {
    int ivecind = datavec[1].fVecShapeIndex[iq].first;
    int ishapeind = datavec[1].fVecShapeIndex[iq].second;
    
    phiVecHCurl(iq , 0) = phiQ(ishapeind , 0) * datavec[1].fNormalVec(0 , ivecind);
    phiVecHCurl(iq , 1) = phiQ(ishapeind , 0) * datavec[1].fNormalVec(1 , ivecind);
    phiVecHCurl(iq , 2) = phiQ(ishapeind , 0) * datavec[1].fNormalVec(2 , ivecind);
//    std::cout<<"shape n: "<<ishapeind<<std::endl;
//    std::cout<<"vector"<<std::endl
//    <<std::setw(10)<<datavec[1].fNormalVec(0 , ivecind)<<" "
//    <<std::setw(10)<<datavec[1].fNormalVec(1 , ivecind)<<" "
//    <<std::setw(10)<<datavec[1].fNormalVec(2 , ivecind)<<std::endl;
    std::cout<<"function: "<<ishapeind<<std::endl
    <<std::setw(10)<<phiVecHCurl( iq , 0 )<<" "
    <<std::setw(10)<<phiVecHCurl( iq , 1 )<<" "
    <<std::setw(10)<<phiVecHCurl( iq , 2 )<<std::endl;
  }
  
  /*********************COMPUTE CURL****************************/
  TPZFMatrix<REAL> &dphiQdaxes = datavec[1].dphix;
  TPZFNMatrix<3,REAL> dphiQ;
  TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, datavec[1].axes);
  TPZFNMatrix<3,REAL> gradScalarPhi(phrq , 3 , 0.);
  TPZFNMatrix<3,REAL> ivecHCurl(phrq , 3 , 0.);
  
  TPZManVector<REAL,3> iVecCurl(3,0.);
  for (int iPhi = 0; iPhi < phrq; iPhi++) {
    int ivecind = datavec[1].fVecShapeIndex[iPhi].first;
    int ishapeind = datavec[1].fVecShapeIndex[iPhi].second;
    iVecCurl[0] = datavec[1].fNormalVec(0,ivecind);
    iVecCurl[1] = datavec[1].fNormalVec(1,ivecind);
    iVecCurl[2] = datavec[1].fNormalVec(2,ivecind);
    
    for (int i = 0; i<dphiQ.Rows(); i++) {
      gradScalarPhi(iPhi,i) = dphiQ(i,ishapeind);
      ivecHCurl(iPhi,i) = iVecCurl[i];
    }
  }
  TPZFNMatrix<40,REAL> curlPhi;
  ComputeCurl(gradScalarPhi, ivecHCurl, curlPhi);
  
  const STATE muR =  fUr(x);
  const STATE epsilonR = fEr(x);
  REAL k0 = fW*sqrt(M_EZERO*M_UZERO);
  /*****************ACTUAL COMPUTATION OF CONTRIBUTION*****************/
  
  int nHCurlFunctions  = phrq;
  for (int iq = 0; iq < nHCurlFunctions; iq++ ) {
    ef(iq,0) = 0.;
    for (int jq = 0 ; jq < nHCurlFunctions; jq++) {
      
      
      STATE curlIdotCurlJ = 0.;
      curlIdotCurlJ += curlPhi(iq , 0) * curlPhi(jq , 0);
      curlIdotCurlJ += curlPhi(iq , 1) * curlPhi(jq , 1);
      curlIdotCurlJ += curlPhi(iq , 2) * curlPhi(jq , 2);
      
      STATE curlIXStar = 0., curlJX = 0.;
      //it is needed to set curlE dot x = j*k0*sin(theta)*Ez
      curlIXStar =  -1. * imaginary * k0 * sin(fTheta) * phiVecHCurl(iq , 2);
      curlJX = 1. * imaginary * k0 * sin(fTheta) * phiVecHCurl(jq , 2);
      curlIdotCurlJ += curlIXStar * curlJX;
      
      STATE phiIdotPhiJ = 0.;
      phiIdotPhiJ += phiVecHCurl(iq , 0) * phiVecHCurl(jq , 0);
      phiIdotPhiJ += phiVecHCurl(iq , 1) * phiVecHCurl(jq , 1);
      phiIdotPhiJ += phiVecHCurl(iq , 2) * phiVecHCurl(jq , 2);
      
      ek(iq,jq)+= 1./muR * curlIdotCurlJ * weight;
      ek(iq,jq)+= -1. * k0 * k0 /(fScale * fScale) * epsilonR * phiIdotPhiJ * weight;
    }
  }
}

void TPZMatMFHCurlFran::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatMFHCurlFran::ContributeForcingRTBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
  return;
}

void TPZMatMFHCurlFran::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	// Setting the phis
	TPZFMatrix<REAL> &phiQ = data.phi;
	
	int nshape=phiQ.Rows();
	REAL BIG = TPZMaterial::gBigNumber;

	const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
	const STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
	
	switch ( bc.Type() )
	{
		case 0:
			for(int i = 0 ; i<nshape ; i++)
			{
				const STATE rhs = phiQ(i,0) * BIG  * v2;
				ef(i,0) += rhs*weight;
				for(int j=0;j<nshape;j++)
				{
          const STATE stiff = phiQ(i,0) * phiQ(j,0) * BIG ;
					ek(i,j) += stiff*weight;
				}
			}
			break;
		case 1:
			DebugStop();
			break;
		case 2:
			for(int i = 0 ; i<nshape ; i++)
			{
				STATE rhs = phiQ(i,0);
				rhs *= v2/fScale;
				ef(i,0) += rhs*weight;
				for(int j=0;j<nshape;j++)
				{
					STATE stiff = phiQ(i,0) *  phiQ(j,0);
					stiff *= v1/fScale;
					ek(i,j) += stiff*weight;
				}
			}
			break;
	}
}

void TPZMatMFHCurlFran::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	DebugStop();
}

void TPZMatMFHCurlFran::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatMFHCurlFran::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	DebugStop();
}

int TPZMatMFHCurlFran::IntegrationRuleOrder(int elPMaxOrder) const
{
	return 2+elPMaxOrder*2;
}


int TPZMatMFHCurlFran::VariableIndex(const std::string &name)
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
int TPZMatMFHCurlFran::NSolutionVariables(int var)
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
void TPZMatMFHCurlFran::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
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


STATE urDefault( const TPZVec<REAL> &x )
{
	return 1.0;
}

STATE erDefault( const TPZVec<REAL> &x )
{
	return 1.0;
}









