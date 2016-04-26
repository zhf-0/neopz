#include "TPZMatMFHCurlFran.h"

#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif




TPZMatMFHCurlFran::TPZMatMFHCurlFran(int id, REAL lambda, REAL kz , STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &), REAL e0, REAL t, REAL scale) : TPZVecL2(id), fUr(ur), fEr(er), fLambda(lambda), fKz (kz) , fE0 (e0) , fTheta(t) , fScale(scale)
{
	fW=2.*M_PI*M_C/fLambda;
}

TPZMatMFHCurlFran::TPZMatMFHCurlFran(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault), fLambda(1.55e-9)
{
	fW=2.*M_PI*M_C/fLambda;
	fTheta = 0.;
  fKz = -999;
  fE0 = 1;
  fScale = 1;
}

/** @brief Default constructor */
TPZMatMFHCurlFran::TPZMatMFHCurlFran() : TPZVecL2(), fUr(urDefault),
fEr(erDefault), fLambda(1.55e-9)
{
	fW=2.*M_PI*M_C/fLambda;
	fTheta = 0.;
  fKz = -999;
  fE0 = 1;
  fScale = 1;
}


TPZMatMFHCurlFran::TPZMatMFHCurlFran(const TPZMatMFHCurlFran &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr)
{
	fLambda = mat.fLambda;
	fTheta = mat.fTheta;
	fE0 = mat.fE0;
	fScale = mat.fScale;
  fKz = mat.fKz;
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
  /*********************CREATE H1 FUNCTIONS****************************/
  TPZFNMatrix<12,REAL> phiH1 = datavec[0].phi;
  TPZFNMatrix<36,REAL> dphiH1daxes = datavec[0].dphix;
  TPZFNMatrix<3,REAL> dphiH1;
  TPZAxesTools<REAL>::Axes2XYZ(dphiH1daxes, dphiH1, datavec[1].axes);
  TPZFNMatrix<3,REAL> gradPhiH1(phiH1.Rows() , 3 , 0.);
  for ( int iFunc = 0 ; iFunc < phiH1.Rows(); iFunc++ ) {
    gradPhiH1 ( iFunc , 0 ) = dphiH1 ( 0 , iFunc );
    gradPhiH1 ( iFunc , 1 ) = dphiH1 ( 1 , iFunc );
  }

  
  /*********************CREATE HCURL FUNCTIONS****************************/
  TPZFNMatrix<12,REAL> phiScaHCurl = datavec[1].phi;
  TPZManVector<REAL,3> x = datavec[0].x;
  
  int phrq = datavec[1].fVecShapeIndex.NElements();
  std::cout<<"x"<<std::endl<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
  
  TPZFNMatrix< 36 , REAL > phiVecHCurl(phrq , 3 , 0.);
  for (int iq = 0 ; iq < phrq ; iq++) {
    int ivecind = datavec[1].fVecShapeIndex[iq].first;
    int ishapeind = datavec[1].fVecShapeIndex[iq].second;
    
    phiVecHCurl(iq , 0) = phiScaHCurl(ishapeind , 0) * datavec[1].fNormalVec(0 , ivecind);
    phiVecHCurl(iq , 1) = phiScaHCurl(ishapeind , 0) * datavec[1].fNormalVec(1 , ivecind);
    phiVecHCurl(iq , 2) = phiScaHCurl(ishapeind , 0) * datavec[1].fNormalVec(2 , ivecind);
//    std::cout<<"shape n: "<<ishapeind<<std::endl;
//    std::cout<<"vector"<<std::endl
//    <<std::setw(10)<<datavec[1].fNormalVec(0 , ivecind)<<" "
//    <<std::setw(10)<<datavec[1].fNormalVec(1 , ivecind)<<" "
//    <<std::setw(10)<<datavec[1].fNormalVec(2 , ivecind)<<std::endl;
//    std::cout<<"function: "<<ishapeind<<std::endl
//    <<std::setw(10)<<phiVecHCurl( iq , 0 )<<" "
//    <<std::setw(10)<<phiVecHCurl( iq , 1 )<<" "
//    <<std::setw(10)<<phiVecHCurl( iq , 2 )<<std::endl;
  }
  
  /*********************COMPUTE CURL****************************/
  TPZFMatrix<REAL> &dphiQdaxes = datavec[1].dphix;
  TPZFNMatrix<3,REAL> dphiQ;
  TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, datavec[1].axes);
  TPZFNMatrix<3,REAL> gradPhiForHCurl(phrq , 3 , 0.);
  TPZFNMatrix<3,REAL> ivecHCurl(phrq , 3 , 0.);
  
  TPZManVector<REAL,3> iVecCurl(3,0.);
  for (int iPhi = 0; iPhi < phrq; iPhi++) {
    int ivecind = datavec[1].fVecShapeIndex[iPhi].first;
    int ishapeind = datavec[1].fVecShapeIndex[iPhi].second;
    iVecCurl[0] = datavec[1].fNormalVec(0,ivecind);
    iVecCurl[1] = datavec[1].fNormalVec(1,ivecind);
    iVecCurl[2] = datavec[1].fNormalVec(2,ivecind);
    
    for (int i = 0; i<dphiQ.Rows(); i++) {
      gradPhiForHCurl(iPhi,i) = dphiQ(i,ishapeind);
      ivecHCurl(iPhi,i) = iVecCurl[i];
    }
  }
  TPZFNMatrix<40,REAL> curlPhi;
  ComputeCurl(gradPhiForHCurl, ivecHCurl, curlPhi);
  
  const STATE muR =  fUr(x);
  const STATE epsilonR = fEr(x);
  REAL k0 = fW*sqrt(M_EZERO*M_UZERO);
  //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//
  
  const int nHCurlFunctions  = phrq;
  const int nH1Functions  = phiH1.Rows();
  STATE stiff = 0.;
  STATE kz = fKz;
  for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
//    if (iVec == nHCurlFunctions - 1){
//      for (int iPrint = 0; iPrint < nHCurlFunctions; iPrint++) {
//        std::cout<<"function: "<<iPrint<<std::endl
//        <<std::setw(10)<<phiVecHCurl( iPrint , 0 )<<" "
//        <<std::setw(10)<<phiVecHCurl( iPrint , 1 )<<" "
//        <<std::setw(10)<<phiVecHCurl( iPrint , 2 )<<std::endl;
//        std::cout<<"curl: "<<iVec<<std::endl
//        <<std::setw(10)<<curlPhi( iPrint , 0 )<<" "
//        <<std::setw(10)<<curlPhi( iPrint , 1 )<<" "
//        <<std::setw(10)<<curlPhi( iPrint , 2 )<<std::endl;
//      }
//      
//    }
    for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
      STATE curlIdotCurlJ = 0.;
      curlIdotCurlJ += std::conj( curlPhi(iVec , 0) ) * curlPhi(jVec , 0);
      curlIdotCurlJ += std::conj( curlPhi(iVec , 1) )* curlPhi(jVec , 1);
      curlIdotCurlJ += std::conj( curlPhi(iVec , 2) ) * curlPhi(jVec , 2);
      STATE phiIdotPhiJ = 0.;
      phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 0) ) * phiVecHCurl(jVec , 0);
      phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 1) ) * phiVecHCurl(jVec , 1);
      phiIdotPhiJ += std::conj( phiVecHCurl(iVec , 2) ) * phiVecHCurl(jVec , 2);
      
      stiff = 1./muR * curlIdotCurlJ;
      stiff -= k0 * k0 * epsilonR * phiIdotPhiJ;
      stiff += kz * kz * 1./muR * phiIdotPhiJ;
      
      ek( iVec , jVec ) += stiff * datavec[0].detjac * weight ;
      
    }
    for (int jSca = 0; jSca < nH1Functions; jSca++) {
//      std::cout<<"function: "<<jSca<<std::endl
//      <<std::setw(10)<<phiSca( jSca , 0 )<<std::endl;
//      std::cout<<"grad: "<<std::endl
//      <<std::setw(10)<<gradPhiSca( jSca , 0 )<<" "
//      <<std::setw(10)<<gradPhiSca( jSca , 1 )<<" "
//      <<std::setw(10)<<gradPhiSca( jSca , 2 )<<std::endl;
      STATE phiVecDotGradPhiSca = 0.;
      phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 0) ) * gradPhiH1(jSca , 0);
      phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 1) ) * gradPhiH1(jSca , 1);
      phiVecDotGradPhiSca += std::conj( phiVecHCurl(iVec , 2) ) * gradPhiH1(jSca , 2);
      
      stiff = kz * kz * 1./muR * phiVecDotGradPhiSca;
      
      ek( iVec , nHCurlFunctions + jSca ) += stiff * datavec[0].detjac * weight ;
    }
  }
  for (int iSca = 0; iSca < nH1Functions; iSca++) {
    for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
      STATE phiVecDotGradPhiSca = 0.;
      phiVecDotGradPhiSca += phiVecHCurl(jVec , 0) * std::conj( gradPhiH1(iSca , 0) );
      phiVecDotGradPhiSca += phiVecHCurl(jVec , 1) * std::conj( gradPhiH1(iSca , 1) );
      phiVecDotGradPhiSca += phiVecHCurl(jVec , 2) * std::conj( gradPhiH1(iSca , 2) );
      
      stiff = kz * kz * 1./muR * phiVecDotGradPhiSca;

      ek( nHCurlFunctions + iSca , jVec ) += stiff * datavec[0].detjac * weight ;
    }
    for (int jSca = 0; jSca < nH1Functions; jSca++) {
      STATE gradPhiScaDotGradPhiSca = 0.;
      gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 0) ) * gradPhiH1(jSca , 0);
      gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 1) ) * gradPhiH1(jSca , 1);
      gradPhiScaDotGradPhiSca += std::conj( gradPhiH1(iSca , 2) ) * gradPhiH1(jSca , 2);
    
      stiff = kz * kz * 1./muR * gradPhiScaDotGradPhiSca;
      stiff -= kz * kz * k0 * k0 * epsilonR * std::conj( phiH1( iSca , 0 ) ) * phiH1( jSca , 0 );
      ek( nHCurlFunctions + iSca , nHCurlFunctions + jSca ) += stiff * datavec[0].detjac * weight ;
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
  DebugStop();
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











