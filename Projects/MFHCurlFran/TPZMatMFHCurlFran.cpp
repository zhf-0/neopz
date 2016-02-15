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

void TPZMatMFHCurlFran::RotateForHCurl(TPZVec<REAL> normal , TPZFMatrix<REAL> vHdiv , TPZFMatrix<REAL> &vHcurl ){
  int nFunctions = vHdiv.Rows();
  vHcurl.Resize( vHdiv.Rows(), vHdiv.Cols());
  vHcurl.Zero();
  
  for (int i = 0 ; i < nFunctions; i++) {
    
    vHcurl(i,0) = normal[1]*vHdiv(i,2) - vHdiv(i,1)*normal[2];
    vHcurl(i,1) = normal[2]*vHdiv(i,0) - vHdiv(i,2)*normal[0];
    vHcurl(i,2) = normal[0]*vHdiv(i,1) - vHdiv(i,0)*normal[1];
  }
}
void TPZMatMFHCurlFran::ContributeForcingRT(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  TPZFNMatrix<12,REAL> phiQ = data.phi;
  TPZManVector<REAL,3> x = data.x;
  
  int phrq = data.fVecShapeIndex.NElements();
  
  
  /*********************CREATE HDIV FUNCTIONS****************************/
  TPZFNMatrix< 36 , REAL > phiVecHDiv(12 , 3 , 0.);
  for (int iq = 0 ; iq < phrq ; iq++) {
    int ivecind = data.fVecShapeIndex[iq].first;
    int ishapeind = data.fVecShapeIndex[iq].second;
    
    phiVecHDiv(iq , 0) = phiQ(ishapeind , 0) * data.fNormalVec(0 , ivecind);
    phiVecHDiv(iq , 1) = phiQ(ishapeind , 0) * data.fNormalVec(1 , ivecind);
    phiVecHDiv(iq , 2) = phiQ(ishapeind , 0) * data.fNormalVec(2 , ivecind);
  }
  /*********************CREATE RT FUNCTIONS****************************/
  TPZFNMatrix< 12 , REAL > phiVecHDivRT(4 , 3 , 0.);
  for (int iq = 0;  iq < 4; iq++) {
    phiVecHDivRT(iq, 0) = phiVecHDiv(2*iq,0) + phiVecHDiv(2*iq + 1,0);
    phiVecHDivRT(iq, 1) = phiVecHDiv(2*iq,1) + phiVecHDiv(2*iq + 1,1);
    phiVecHDivRT(iq, 2) = phiVecHDiv(2*iq,2) + phiVecHDiv(2*iq + 1,2);
  }
  
  /*********************ROTATE FOR HCURL****************************/
  TPZManVector<REAL,3> ax1(3),ax2(3), elNormal(3);
  for (int i=0; i<3; i++) {
    ax1[i] = data.axes(0,i);//ELEMENTO DEFORMADO
    ax2[i] = data.axes(1,i);//ELEMENTO DEFORMADO
  }
  Cross(ax1, ax2, elNormal);
  
  TPZFNMatrix< 12 , REAL > phiVecHCurlNed(4 , 3 , 0.);
  RotateForHCurl(elNormal , phiVecHDivRT , phiVecHCurlNed);

  int nHCurlFunctions  = phiVecHCurlNed.Rows();
  for (int iq = 0; iq < nHCurlFunctions; iq++ ) {
    for (int jq = 0 ; jq < nHCurlFunctions; jq++) {
      
      STATE phiIdotPhiJ = 0.;
      phiIdotPhiJ += phiVecHCurlNed(iq , 0) * phiVecHCurlNed(jq , 0);
      phiIdotPhiJ += phiVecHCurlNed(iq , 1) * phiVecHCurlNed(jq , 1);
      phiIdotPhiJ += phiVecHCurlNed(iq , 2) * phiVecHCurlNed(jq , 2);
      ek(iq,jq)+= phiIdotPhiJ * weight;
    }
  }
  
}

void TPZMatMFHCurlFran::ContributeValidateFunctions(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  TPZFNMatrix<12,REAL> phiQ = data.phi;
  TPZManVector<REAL,3> x = data.x;
  //  for (int i=0; i<3; i++) {
  //    x[i] /= fScale;
  //  }
  int phrq = data.fVecShapeIndex.NElements();
  phiQ.Print("scalar functions");
  
  
  /*********************CREATE HDIV FUNCTIONS****************************/
  
  TPZFNMatrix< 36 , REAL > phiVecHDiv(phrq , 3 , 0.);
  for (int iq = 0 ; iq < phrq ; iq++) {
    int ivecind = data.fVecShapeIndex[iq].first;
    int ishapeind = data.fVecShapeIndex[iq].second;
    
    phiVecHDiv(iq , 0) = phiQ(ishapeind , 0) * data.fNormalVec(0 , ivecind);
    phiVecHDiv(iq , 1) = phiQ(ishapeind , 0) * data.fNormalVec(1 , ivecind);
    phiVecHDiv(iq , 2) = phiQ(ishapeind , 0) * data.fNormalVec(2 , ivecind);
    std::cout<<"shape n: "<<ishapeind<<std::endl;
    std::cout<<"vector"<<std::endl
    <<std::setw(10)<<data.fNormalVec(0 , ivecind)<<" "
    <<std::setw(10)<<data.fNormalVec(1 , ivecind)<<" "
    <<std::setw(10)<<data.fNormalVec(2 , ivecind)<<std::endl;
  }
  
  /*********************ROTATE FOR HCURL****************************/
  TPZManVector<REAL,3> ax1(3),ax2(3), elNormal(3);
  for (int i=0; i<3; i++) {
    ax1[i] = data.axes(0,i);//ELEMENTO DEFORMADO
    ax2[i] = data.axes(1,i);//ELEMENTO DEFORMADO
  }
  Cross(ax1, ax2, elNormal);
  std::cout<<"determinante da jacobiana"<<std::endl;
  std::cout<<data.detjac<<std::endl;
  std::cout<<"normal :"<<sqrt(elNormal[0] * elNormal[0] + elNormal[1] * elNormal[1] + elNormal[2] * elNormal[2])<<std::endl;
  TPZFNMatrix< 12 , REAL > phiVecHCurl(phrq , 3 , 0.);
  RotateForHCurl(elNormal , phiVecHDiv , phiVecHCurl);
  std::cout<<"integration point:"<<std::endl;
  std::cout<<data.xParametric<<std::endl;
  std::cout<<data.x<<std::endl;
	for (int i = 0 ; i < phrq/2 ; i++) {
		std::cout<<"phiNed("<<i<<"):"
		<<std::setw(10)<<phiVecHCurl(2 * i , 0)+phiVecHCurl(2 * i + 1 , 0)<<" "
		<<std::setw(10)<<phiVecHCurl(2 * i , 1)+phiVecHCurl(2 * i + 1 , 1)<<" "
		<<std::setw(10)<<phiVecHCurl(2 * i , 2)+phiVecHCurl(2 * i + 1 , 2)<<std::endl;
	}
  /*********************COMPUTE CURL****************************/
  TPZFMatrix<REAL> &dphiQdaxes = data.dphix;
  TPZFNMatrix<3,REAL> dphiQ;
  TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);
  TPZFNMatrix<3,REAL> gradScalarPhi(phrq , 3 , 0.);
  TPZFNMatrix<3,REAL> ivecHCurl(phrq , 3 , 0.);
  
  TPZManVector<REAL,3> iVecHDiv(3,0.), ivecForCurl(3,0.);
  for (int iPhi = 0; iPhi < phrq; iPhi++) {
    int ivecind = data.fVecShapeIndex[iPhi].first;
    int ishapeind = data.fVecShapeIndex[iPhi].second;
    iVecHDiv[0] = data.fNormalVec(0,ivecind);
    iVecHDiv[1] = data.fNormalVec(1,ivecind);
    iVecHDiv[2] = data.fNormalVec(2,ivecind);

    Cross(elNormal, iVecHDiv, ivecForCurl);
    for (int i = 0; i<dphiQ.Rows(); i++) {
      gradScalarPhi(iPhi,i) = dphiQ(i,ishapeind);
      ivecHCurl(iPhi,i) = ivecForCurl[i];
    }
  }
  TPZFNMatrix<40,REAL> curlPhi;
  ComputeCurl(gradScalarPhi, ivecHCurl, curlPhi);
  
  /*****************ACTUAL COMPUTATION OF CONTRIBUTION*****************/
  
  int nHCurlFunctions  = phiVecHCurl.Rows();
  for (int iq = 0; iq < nHCurlFunctions; iq++ ) {
    for (int jq = 0 ; jq < nHCurlFunctions; jq++) {
      
      
      STATE curlIdotCurlJ = 0.;
      curlIdotCurlJ += curlPhi(iq , 0) * curlPhi(jq , 0);
      curlIdotCurlJ += curlPhi(iq , 1) * curlPhi(jq , 1);
      curlIdotCurlJ += curlPhi(iq , 2) * curlPhi(jq , 2);

      
      STATE phiIdotPhiJ = 0.;
      phiIdotPhiJ += phiVecHCurl(iq , 0) * phiVecHCurl(jq , 0);
      phiIdotPhiJ += phiVecHCurl(iq , 1) * phiVecHCurl(jq , 1);
      phiIdotPhiJ += phiVecHCurl(iq , 2) * phiVecHCurl(jq , 2);
      
      ek(iq,jq)+= phiIdotPhiJ * weight;
//      ek(iq,jq)+= curlIdotCurlJ * weight;
    }
  }
}



void TPZMatMFHCurlFran::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
//  ContributeForcingRT(data, weight, ek, ef);
//  return

  ContributeValidateFunctions(data, weight, ek, ef);
  return;
	TPZFNMatrix<12,REAL> phiQ = data.phi;
	TPZManVector<REAL,3> x = data.x;
	//  for (int i=0; i<3; i++) {
	//    x[i] /= fScale;
	//  }
	int phrq = data.fVecShapeIndex.NElements();
	phiQ.Print("scalar functions");
	
	
	/*********************CREATE HDIV FUNCTIONS****************************/
	
	TPZFNMatrix< 36 , REAL > phiVecHDiv(phrq , 3 , 0.);
	for (int iq = 0 ; iq < phrq ; iq++) {
		int ivecind = data.fVecShapeIndex[iq].first;
		int ishapeind = data.fVecShapeIndex[iq].second;
		
		phiVecHDiv(iq , 0) = phiQ(ishapeind , 0) * data.fNormalVec(0 , ivecind);
		phiVecHDiv(iq , 1) = phiQ(ishapeind , 0) * data.fNormalVec(1 , ivecind);
		phiVecHDiv(iq , 2) = phiQ(ishapeind , 0) * data.fNormalVec(2 , ivecind);
		std::cout<<"shape n: "<<ishapeind<<std::endl;
		std::cout<<"vector"<<std::endl
		<<std::setw(10)<<data.fNormalVec(0 , ivecind)<<" "
		<<std::setw(10)<<data.fNormalVec(1 , ivecind)<<" "
		<<std::setw(10)<<data.fNormalVec(2 , ivecind)<<std::endl;
	}
	
	/*********************ROTATE FOR HCURL****************************/
	TPZManVector<REAL,3> ax1(3),ax2(3), elNormal(3);
	for (int i=0; i<3; i++) {
		ax1[i] = data.axes(0,i);//ELEMENTO DEFORMADO
		ax2[i] = data.axes(1,i);//ELEMENTO DEFORMADO
	}
	Cross(ax1, ax2, elNormal);
	std::cout<<"determinante da jacobiana"<<std::endl;
	std::cout<<data.detjac<<std::endl;
	std::cout<<"normal :"<<sqrt(elNormal[0] * elNormal[0] + elNormal[1] * elNormal[1] + elNormal[2] * elNormal[2])<<std::endl;
	TPZFNMatrix< 12 , REAL > phiVecHCurl(phrq , 3 , 0.);
	RotateForHCurl(elNormal , phiVecHDiv , phiVecHCurl);
	std::cout<<"integration point:"<<std::endl;
	std::cout<<data.xParametric<<std::endl;
	std::cout<<data.x<<std::endl;
	for (int i = 0 ; i < phrq/2 ; i++) {
		std::cout<<"phiNed("<<i<<"):"
		<<std::setw(10)<<phiVecHCurl(2 * i , 0)+phiVecHCurl(2 * i + 1 , 0)<<" "
		<<std::setw(10)<<phiVecHCurl(2 * i , 1)+phiVecHCurl(2 * i + 1 , 1)<<" "
		<<std::setw(10)<<phiVecHCurl(2 * i , 2)+phiVecHCurl(2 * i + 1 , 2)<<std::endl;
	}
	/*********************COMPUTE CURL****************************/
	TPZFMatrix<REAL> &dphiQdaxes = data.dphix;
	TPZFNMatrix<3,REAL> dphiQ;
	TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);
	TPZFNMatrix<3,REAL> gradScalarPhi(phrq , 3 , 0.);
	TPZFNMatrix<3,REAL> ivecHCurl(phrq , 3 , 0.);
	
	TPZManVector<REAL,3> iVecHDiv(3,0.), ivecForCurl(3,0.);
	for (int iPhi = 0; iPhi < phrq; iPhi++) {
		int ivecind = data.fVecShapeIndex[iPhi].first;
		int ishapeind = data.fVecShapeIndex[iPhi].second;
		iVecHDiv[0] = data.fNormalVec(0,ivecind);
		iVecHDiv[1] = data.fNormalVec(1,ivecind);
		iVecHDiv[2] = data.fNormalVec(2,ivecind);
		
		Cross(elNormal, iVecHDiv, ivecForCurl);
		for (int i = 0; i<dphiQ.Rows(); i++) {
			gradScalarPhi(iPhi,i) = dphiQ(i,ishapeind);
			ivecHCurl(iPhi,i) = ivecForCurl[i];
		}
	}
	TPZFNMatrix<40,REAL> curlPhi;
	ComputeCurl(gradScalarPhi, ivecHCurl, curlPhi);
	
	const STATE muR =  fUr(x);
	const STATE epsilonR = fEr(x);
	REAL k0 = fW*sqrt(M_EZERO*M_UZERO);
	/*****************ACTUAL COMPUTATION OF CONTRIBUTION*****************/
	
	int nHCurlFunctions  = phiVecHCurl.Rows();
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

void TPZMatMFHCurlFran::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	DebugStop();
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









