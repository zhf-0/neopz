#include "TPZMatValidacaoHCurlFran2.h"
#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif

TPZMatValidacaoHCurlFran2::TPZMatValidacaoHCurlFran2(int id, REAL lambda, STATE ( &ur)( const TPZVec<REAL> &),STATE ( &er)( const TPZVec<REAL> &), REAL e0, REAL t, REAL scale) : TPZVecL2(id), fUr(ur), fEr(er), fLambda(lambda), fE0 (e0) , fTheta(t) , fScale(scale)
{
	fW=2.*M_PI*M_C/fLambda;
}

TPZMatValidacaoHCurlFran2::TPZMatValidacaoHCurlFran2(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault), fLambda(1.55e-9)
{
	fW=2.*M_PI*M_C/fLambda;
	fTheta = 0.;
}

/** @brief Default constructor */
TPZMatValidacaoHCurlFran2::TPZMatValidacaoHCurlFran2() : TPZVecL2(), fUr(urDefault),
fEr(erDefault), fLambda(1.55e-9)
{
	fW=2.*M_PI*M_C/fLambda;
	fTheta = 0.;
}


TPZMatValidacaoHCurlFran2::TPZMatValidacaoHCurlFran2(const TPZMatValidacaoHCurlFran2 &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr)
{
	fLambda = mat.fLambda;
	fTheta = mat.fTheta;
	fE0 = mat.fE0;
	fScale = mat.fScale;
	fW=2.*M_PI*M_C/fLambda;
}

TPZMatValidacaoHCurlFran2::~TPZMatValidacaoHCurlFran2()
{
	
}
void TPZMatValidacaoHCurlFran2::RotateForHCurl(TPZVec<REAL> normal , TPZFMatrix<REAL> vHdiv , TPZFMatrix<REAL> &vHcurl ){
  int nFunctions = vHdiv.Rows();
  vHcurl.Resize( vHdiv.Rows(), vHdiv.Cols());
  vHcurl.Zero();
  TPZManVector<REAL,3> result(3,0.);
  
  for (int i = 0 ; i < nFunctions; i++) {
    
    vHcurl(i,0) = normal[1]*vHdiv(i,2) - vHdiv(i,1)*normal[2];
    vHcurl(i,1) = normal[2]*vHdiv(i,0) - vHdiv(i,2)*normal[0];//tava errado
    vHcurl(i,2) = normal[0]*vHdiv(i,1) - vHdiv(i,0)*normal[1];
  }
}
void TPZMatValidacaoHCurlFran2::ContributeOmarVampirinho(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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
  
//  for (int iq = 0 ; iq < phrq ; iq++) {
//    int ivecind = data.fVecShapeIndex[iq].first;
//    int ishapeind = data.fVecShapeIndex[iq].second;
//    
//    TPZManVector<REAL,3> ivecHDiv(3);
//    //		norm = 0.;
//    for(int id=0; id<3; id++){
//      ivecHDiv[id] = data.fNormalVec(id,ivecind);//JA EM XYZ
//      //			norm += ivecHDiv[id] * ivecHDiv[id];
//    }
//    
//    
//    
//    STATE dotProduct = 0;
//    for (int iDot = 0; iDot < 3 ; iDot ++) {
//      dotProduct+=ivecHDiv[iDot] * phiQ(ishapeind , 0) * ivecHDiv[iDot] * phiQ(ishapeind , 0);
//      //      dotProduct+=ivecHDiv[iDot] * phiQ(ishapeind , 0) * ivecHDiv[iDot] * phiQ(ishapeind , 0);
//  		}
//    ek(iq,iq) += dotProduct * weight;
//  }
  
}
void TPZMatValidacaoHCurlFran2::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  ContributeOmarVampirinho(data, weight, ek, ef);
  return;
  ////dot product using hdiv functions
//  if(data.detjac != 1){
//    DebugStop();
//  }
//  
//  TPZFNMatrix<12,REAL> phiQ = data.phi;
//  TPZManVector<REAL,3> x = data.x;
//  TPZManVector<REAL,3> xParametric = data.xParametric;
//  if( fabs(xParametric[0] - 0.861136 ) < 0.01 ){
//    std::cout<<xParametric<<std::endl;
//    std::cout<<phiQ<<std::endl;
//  }
//  
//  int phrq = data.fVecShapeIndex.NElements();
//  for (int iq = 0 ; iq < phrq ; iq++) {
//    int ivecind = data.fVecShapeIndex[iq].first;
//    int ishapeind = data.fVecShapeIndex[iq].second;
//    int ivecindPair = data.fVecShapeIndex[iq + 1*((iq%2+1)%2) - 1*(iq%2)].first;
//    int ishapeindPair = data.fVecShapeIndex[iq + 1*((iq%2+1)%2) - 1*(iq%2)].second;
//    
//    TPZManVector<REAL,3> ivecHDiv(3), ivecHDivPair(3);
//    //		norm = 0.;
//    for(int id=0; id<3; id++){
//      ivecHDiv[id] = data.fNormalVec(id,ivecind);//JA EM XYZ
//      ivecHDivPair[id] = data.fNormalVec(id,ivecindPair);//JA EM XYZ
//      //			norm += ivecHDiv[id] * ivecHDiv[id];
//    }
////    std::cout<<"index = "<<iq<<std::endl;
////    std::cout<<"vecIndex = "<<ivecind<<std::endl;
////    std::cout<<"shapeIndex = "<<ishapeind<<std::endl;
////    std::cout<<ivecHDiv<<std::endl;
//    STATE dotProduct = 0;
//    for (int iDot = 0; iDot < 3 ; iDot ++) {
//      dotProduct+=ivecHDiv[iDot] * phiQ(ishapeind , 0) * ivecHDiv[iDot] * phiQ(ishapeind , 0) + ivecHDivPair[iDot] * phiQ(ishapeindPair , 0) * ivecHDiv[iDot] * phiQ(ishapeind , 0);
////      dotProduct+=ivecHDiv[iDot] * phiQ(ishapeind , 0) * ivecHDiv[iDot] * phiQ(ishapeind , 0);
//		}
//    ek(iq,iq) += dotProduct * weight;
//  }

  
//  TPZFNMatrix<220,REAL> &phiQ = data.phi;
//  TPZFNMatrix<660,REAL> &dphiQdaxes = data.dphix;//ELEMENTO DEFORMADO
//  TPZFNMatrix<3,REAL> dphiQ;
//  TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);
//  TPZManVector<REAL,3> x = data.x;
//  for (int i=0; i<3; i++) {
//    		x[i] /= fScale;
//    	}
//  const STATE muR =  fUr(x);
//  const STATE epsilonR = fEr(x);
//  REAL k0 = fW*sqrt(M_EZERO*M_UZERO);
//  TPZManVector<STATE,3> force(3);
//  if(fForcingFunction) {
//    fForcingFunction->Execute(x,force);
//  }
//  
//  TPZManVector<REAL,3> ax1(3),ax2(3), normal(3);
//  for (int i=0; i<3; i++) {
//    ax1[i] = data.axes(0,i);//ELEMENTO DEFORMADO
//    ax2[i] = data.axes(1,i);//ELEMENTO DEFORMADO
//  }
//  Cross(ax1, ax2, normal);
//  
//  
//  int phrq = data.fVecShapeIndex.NElements();
//  for (int iq = 0 ; iq < phrq ; iq++) {
//    int ivecind = data.fVecShapeIndex[iq].first;
//    int ishapeind = data.fVecShapeIndex[iq].second;
//    int ivecindPair = data.fVecShapeIndex[iq + 1*((iq%2+1)%2) - 1*(iq%2)].first;
//    int ishapeindPair = data.fVecShapeIndex[iq + 1*((iq%2+1)%2) - 1*(iq%2)].second;
//    
//    TPZManVector<REAL,3> ivecHDiv(3), ivecHDivPair(3) , ivecHCurl , ivecHCurlPair;
//    //		norm = 0.;
//    for(int id=0; id<3; id++){
//      ivecHDiv[id] = data.fNormalVec(id,ivecind);//JA EM XYZ
//      ivecHDivPair[id] = data.fNormalVec(id,ivecindPair);//JA EM XYZ
//      //			norm += ivecHDiv[id] * ivecHDiv[id];
//    }
//    Cross(normal, ivecHDiv, ivecHCurl);
//    Cross(normal, ivecHDivPair, ivecHCurlPair);
//    
//    
//    TPZManVector<REAL,3> curlI(3), gradPhiI(3) , curlIPair(3), gradPhiIPair(3);
//    for (int i = 0; i<dphiQ.Rows(); i++) {
//      gradPhiI[i] = dphiQ(i,ishapeind);
//      gradPhiIPair[i] = dphiQ(i,ishapeindPair);
//    }
//    Cross(gradPhiI, ivecHCurl, curlI);
//    Cross(gradPhiIPair, ivecHCurlPair, curlIPair);
//    STATE curlIX;
//    
//    curlIX = -1. * imaginary * k0 * sin(fTheta) * phiQ(ishapeind,0) * ivecHCurl[2] *fScale;
//    curlIX += -1. * imaginary * k0 * sin(fTheta) * phiQ(ishapeindPair,0) * ivecHCurlPair[2] *fScale;
//    
//    for (int jq = 0 ; jq < phrq ; jq++) {
//      int jvecind = data.fVecShapeIndex[jq].first;
//      int jshapeind = data.fVecShapeIndex[jq].second;
//      int jvecindPair = data.fVecShapeIndex[jq + 1*((jq%2+1)%2) - 1*(jq%2)].first;
//      int jshapeindPair = data.fVecShapeIndex[iq + 1*((jq%2+1)%2) - 1*(jq%2)].second;
//      
//      TPZManVector<REAL,3> jvecHDiv(3), jvecHDivPair(3) , jvecHCurl , jvecHCurlPair;
//      //		norm = 0.;
//      for(int jd=0; jd<3; jd++){
//        jvecHDiv[jd] = data.fNormalVec(jd,jvecind);//JA EM XYZ
//        jvecHDivPair[jd] = data.fNormalVec(jd,jvecindPair);//JA EM XYZ
//        //			norm += ivecHDiv[id] * ivecHDiv[id];
//      }
//      Cross(normal, jvecHDiv, jvecHCurl);
//      Cross(normal, jvecHDivPair, jvecHCurlPair);
//      
//      
//      TPZManVector<REAL,3> curlJ(3), gradPhiJ(3) , curlJPair(3), gradPhiJPair(3);
//      for (int j = 0; j<dphiQ.Rows(); j++) {
//        gradPhiJ[j] = dphiQ(j,jshapeind);
//        gradPhiJPair[j] = dphiQ(j,jshapeindPair);
//      }
//      Cross(gradPhiJ, jvecHCurl, curlJ);
//      Cross(gradPhiJPair, jvecHCurlPair, curlJPair);
//      STATE curlJX;
//      
//      curlJX = -1. * imaginary * k0 * sin(fTheta) * phiQ(jshapeind,0) * jvecHCurl[2] *fScale;
//      curlJX += -1. * imaginary * k0 * sin(fTheta) * phiQ(jshapeindPair,0) * jvecHCurlPair[2] *fScale;
//      /*****************ACTUAL COMPUTATION*********************/
//      STATE phiIdotPhiJ = 0;
//      for (int iDot = 0; iDot < 3 ; iDot ++) {
//        phiIdotPhiJ += ivecHCurl[iDot] * phiQ(ishapeind , 0) * jvecHCurl[iDot] * phiQ(jshapeind , 0)
//        + ivecHCurlPair[iDot] * phiQ(ishapeindPair , 0) * jvecHCurl[iDot] * phiQ(jshapeind , 0);
//      }
//      ek(iq,iq) += phiIdotPhiJ * weight;
//    }
//  }
  
	TPZManVector<REAL,3> x = data.x;
	for (int i=0; i<3; i++) {
		x[i] /= fScale;
	}
	const STATE muR =  fUr(x);
	const STATE epsilonR = fEr(x);
	REAL k0 = fW*sqrt(M_EZERO*M_UZERO);
	TPZManVector<STATE,3> force(3);
	if(fForcingFunction) {
		fForcingFunction->Execute(x,force);
	}
	//AQUIFRAN
	// Setting the phis
	TPZFMatrix<REAL> &phiQ = data.phi;
	TPZFMatrix<REAL> &dphiQdaxes = data.dphix;//ELEMENTO DEFORMADO
	TPZFNMatrix<3,REAL> dphiQ;
	TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);
	
	TPZManVector<REAL,3> ax1(3),ax2(3), normal(3);
	for (int i=0; i<3; i++) {
		ax1[i] = data.axes(0,i);//ELEMENTO DEFORMADO
		ax2[i] = data.axes(1,i);//ELEMENTO DEFORMADO
	}
	Cross(ax1, ax2, normal);
	
//	REAL norm = 0.;
//	for (int i = 0 ; i < normal.size(); i++) {
//		norm += normal[i] * normal[i];
//	}
//	norm = sqrt(norm);
//	for (int i = 0 ; i < normal.size(); i++) {
//		normal[i] /= norm;
//	}
	int phrq;
	phrq = data.fVecShapeIndex.NElements();
#ifdef LOG4CXX
	if(logger->isDebugEnabled() )
	{
		std::stringstream sout;
		sout<<std::endl;
		sout<<"geo el id:"<<data.gelElId<<std::endl;
		sout<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	for(int iq=0; iq<phrq; iq++)
	{
		//ef(iq, 0) += 0.;
		int ivecind = data.fVecShapeIndex[iq].first;
		int ishapeind = data.fVecShapeIndex[iq].second;
		
		TPZManVector<REAL,3> ivecHDiv(3), ivecHCurl(3);
//		norm = 0.;
		for(int id=0; id<3; id++){
			ivecHDiv[id] = data.fNormalVec(id,ivecind);//JA EM XYZ
//			norm += ivecHDiv[id] * ivecHDiv[id];
		}
//		norm = sqrt(norm);
//		for(int id=0; id<3; id++){
//			ivecHDiv[id] /= norm;
//		}
		//ROTATE FOR HCURL
		Cross(normal, ivecHDiv, ivecHCurl);
		
		STATE ff = 0.;
		for (int i=0; i<3; i++) {
			ff += ivecHCurl[i]*force[i];
		}
		//		if(elId==1){
		//			ef(iq,0) += weight*ff*phiQ(ishapeind,0);
		//		}
		//		else{
		//			ef(iq,0) -= weight*ff*phiQ(ishapeind,0);
		//		}
		ef(iq,0) += weight*ff*phiQ(ishapeind,0);
		TPZManVector<REAL,3> curlI(3), gradPhiI(3);
		for (int i = 0; i<dphiQ.Rows(); i++) {
			gradPhiI[i] = dphiQ(i,ishapeind);
		}
		Cross(gradPhiI, ivecHCurl, curlI);
		if (curlI[0] != 0) {
			DebugStop();
		}
		STATE curlIXStar = 0.;
		//eh necessario forcar curlE dot x = j*k0*sin(theta)*Ez
		curlIXStar = -1. * imaginary * k0 * sin(fTheta) * phiQ(ishapeind,0) * ivecHCurl[2] ;
		
		for (int jq=0; jq<phrq; jq++)
		{
			TPZManVector<REAL,3> jvecHDiv(3), jvecHCurl(3);//JA EM XYZ
			int jvecind = data.fVecShapeIndex[jq].first;
			int jshapeind = data.fVecShapeIndex[jq].second;
			
			//norm = 0.;
			for(int id=0; id<3; id++){
				jvecHDiv[id] = data.fNormalVec(id,jvecind);
				//norm += jvecHDiv[id] * jvecHDiv[id];
			}
			//norm = sqrt(norm);
//			for(int id=0; id<3; id++){
//				jvecHDiv[id] /= norm;
//			}
			//ROTATE FOR HCURL
			Cross(normal, jvecHDiv,jvecHCurl);
			TPZManVector<REAL,3> curlJ(3),  gradPhiJ(3);
			
			for (int i = 0; i<dphiQ.Rows(); i++) {
				gradPhiJ[i] = dphiQ(i,jshapeind);
			}
			Cross(gradPhiJ, jvecHCurl, curlJ);
			//eh necessario forcar curlE dot x = j*k0*sin(theta)*Ez
			STATE curlJX;
			
			curlJX = 1. * imaginary * k0 * sin(fTheta) * phiQ(jshapeind,0) * jvecHCurl[2] *fScale;//AQUIFRAN
			
			REAL phiIdotphiJ;
			phiIdotphiJ = ( phiQ(ishapeind,0) * ivecHCurl[0] * phiQ(jshapeind,0) * jvecHCurl[0])
			+ ( phiQ(ishapeind,0) * ivecHCurl[1] * phiQ(jshapeind,0) * jvecHCurl[1])
			+ ( phiQ(ishapeind,0) * ivecHCurl[2] * phiQ(jshapeind,0) * jvecHCurl[2]);
			
			STATE stiff = 0.;
			stiff += (1./muR) * ( curlJX * curlIXStar + (curlJ[0] * curlI[0] + curlJ[1] * curlI[1] + curlJ[2] * curlI[2]) );
			//stiff += (1./muR) * ( (curlJ[0] * curlI[0] + curlJ[1] * curlI[1] + curlJ[2] * curlI[2]) );
			stiff += -1. * k0 * k0 /(fScale * fScale) * epsilonR * phiIdotphiJ;
      stiff = phiIdotphiJ;//AQUIFRAN
			ek(iq,jq) += weight * stiff;
		}
	}
}

void TPZMatValidacaoHCurlFran2::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatValidacaoHCurlFran2::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatValidacaoHCurlFran2::ContributeOmarVampirinhoBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
  return;
}

void TPZMatValidacaoHCurlFran2::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	// Setting the phis
//	TPZFMatrix<REAL> &phiQ = data.phi;
//	
//	int nshape=phiQ.Rows();
//	REAL BIG = TPZMaterial::gBigNumber;
//	BIG=BIG*BIG;
//	const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
//	const STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
//	
//	switch ( bc.Type() )
//	{
//		case 0:
//			for(int i = 0 ; i<nshape ; i++)
//			{
//				const STATE rhs = phiQ(i,0) * BIG * (1. + imaginary ) * v2;
//				ef(i,0) += rhs*weight;
//				for(int j=0;j<nshape;j++)
//				{
//					const STATE stiff = phiQ(i,0) * phiQ(j,0) * BIG * (1. + imaginary );
//					ek(i,j) += stiff*weight;
//				}
//			}
//			break;
//		case 1:
//			DebugStop();
//			break;
//		case 2:
//			for(int i = 0 ; i<nshape ; i++)
//			{
//				STATE rhs = phiQ(i,0);
//				rhs *= v2/fScale;
//				ef(i,0) += rhs*weight;
//				for(int j=0;j<nshape;j++)
//				{
//					STATE stiff = phiQ(i,0) *  phiQ(j,0);
//					stiff *= v1/fScale;
//					ek(i,j) += stiff*weight;
//				}
//			}
//			break;
//	}
  //AQUIFRAN
}

void TPZMatValidacaoHCurlFran2::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	DebugStop();
}

void TPZMatValidacaoHCurlFran2::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMatValidacaoHCurlFran2::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	DebugStop();
}

int TPZMatValidacaoHCurlFran2::IntegrationRuleOrder(int elPMaxOrder) const
{
	return 2+elPMaxOrder*2;
}


int TPZMatValidacaoHCurlFran2::VariableIndex(const std::string &name)
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
int TPZMatValidacaoHCurlFran2::NSolutionVariables(int var)
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
void TPZMatValidacaoHCurlFran2::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
	
	TPZManVector<STATE,3> ax1(3),ax2(3), normal(3);
	for (int i=0; i<3; i++) {
		ax1[i] = data.axes(0,i);
		ax2[i] = data.axes(1,i);
	}
	//ROTATE FOR HCURL
	Cross(ax1, ax2, normal);
//	STATE norm = 0.;
//	for (int i = 0 ; i < normal.size(); i++) {
//		norm += normal[i] * normal[i];
//	}
//	norm = sqrt(norm);
//	for (int i = 0 ; i < normal.size(); i++) {
//		normal[i] /= norm;
//	}
	
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









