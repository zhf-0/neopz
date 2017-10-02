//$Id: TDiscoFunction.cpp,v 1.3 2008-04-04 13:36:47 fortiago Exp $

#include "TPZDiscoFunction.h"

template<class TVar>
TPZDiscoFunction<TVar>::TPZDiscoFunction() : TPZRegisterClassId(&TPZDiscoFunction::ClassId){

}

template<class TVar>
TPZDiscoFunction<TVar>::~TPZDiscoFunction(){

}
    
template<class TVar>
void  TPZDiscoFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df){
  f.Resize(1);
  df.Resize(8,1);
  df.Zero();
  double Xc = 0.;
  double Yc = 0.;
  double Xp = x[0];
  double Yp = x[1];
  double r = sqrt( (Xp-Xc)*(Xp-Xc) +(Yp-Yc)*(Yp-Yc) );
  f[0] = log(r);
  df(0,0) = (Xp - Xc)/(r*r);
  df(1,0) = (Yp - Yc)/(r*r);
  df(2,0) = 0.; //Laplaciano
  ///as outras derivadas nao serao usadas
}  

template<class TVar>
int  TPZDiscoFunction<TVar>::NFunctions()const {
  return 1;
}
  
template<class TVar>
int  TPZDiscoFunction<TVar>::PolynomialOrder() const{
  return 100;
}

template<class TVar>
int TPZDiscoFunction<TVar>::ClassId() {
    //CLASSIDFRANreturn TPZFunction<TVar>::ClassId()^Hash("TPZDiscoFunction");
  return 666;
}

template class TPZDiscoFunction<float>;
template class TPZDiscoFunction<double>;
template class TPZDiscoFunction<long double>;
