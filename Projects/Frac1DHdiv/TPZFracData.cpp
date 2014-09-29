#include "pzlog.h"
#include "TPZFracData.h"

#include <iostream>
#include <cmath>


TPZFracData::TPZFracData()
{

  fmu = 0.;
  fRho = 0.;
  fKab.Redim(0, 0);
  fPhi = 0.;;
  fDeltaT = 0.;
  fTime = 0.;
  fTtot = 0.;
  fTheta = 0.;
  fLfrac = 0.;
  fHf = 0.;
  fE = 0.;
  fnu = 0.;
  fQ = 0.;
  fSigmaConf = 0.;
  fCl = 0.;
  fPe = 0.;
  fPref = 0.;
  fvsp = 0.;
  fAccumVl = 0.;
  fPorderPressure = 0;
  fPorderFlow = 0;
  fnelFrac = 0;
  felSize = 0.;
  fdwdp = 0.;
  fpostProcessFileName = "DefaultName.vtk";
  fDebugMap.clear();
  this->SetCurrentState();
}


TPZFracData::~TPZFracData()
{
  
}


void TPZFracData::Porosity(REAL p, REAL &porosity, REAL &dPorosityDp) const
{
  const REAL Rockcomp = 1.0*(5.09858e-10);// 5.09858e-10 (Pa)-1 -> 50e-6 (kgf/m2)-1
  const REAL pref = (18.0e6);// 1.0 Mpa
  const bool islinear = true;
    if (islinear) {
        // Linear model
        porosity = fPhi*(1.0+Rockcomp*(p-pref));
        dPorosityDp = fPhi*Rockcomp;
    }else
    {
        porosity = fPhi*exp(Rockcomp*(p-pref));
        dPorosityDp = fPhi*Rockcomp*exp(Rockcomp*(p-pref));
    }
}

void TPZFracData::Density(REAL p, REAL &RhoFluid, REAL &dRhoDp) const
{
    // Using constant density
  const REAL Oilcomp = 0.0*(1.01971e-10);// 1.01971e-9 (Pa)-1 -> 10e-6 (kgf/m2)-1
  const REAL pref = (18.0e6);// 1.0 Mpa
  const bool islinear = true;
    if (islinear) {
        // Linear model
        RhoFluid = fRho*(1+Oilcomp*(p-pref));
        dRhoDp = fRho*Oilcomp;
    }else
    {
      RhoFluid = fRho*exp(Oilcomp*(p-pref));
      dRhoDp = fRho*Oilcomp*exp(Oilcomp*(p-pref));
    }
}


void TPZFracData::Viscosity(REAL p, REAL &FluidViscosity, REAL &dFluidViscosityDp) const
{
  // Constant for a while.
  FluidViscosity = fmu;
  dFluidViscosityDp = 0;
}

TPZFMatrix<STATE> TPZFracData::K() const{
  return fKab;
}

TPZFMatrix<STATE> TPZFracData::Kinv() const{
  TPZFMatrix<STATE> Kinverse(2,2,0.0);
  
  TPZFMatrix<STATE> Kab = fKab;
  REAL Constant = (-1.0*Kab(0,1)*Kab(1,0)+Kab(0,0)*Kab(1,1));
  
  if (fKab.Rows() != 2) {
    std::cout << "Absolute permeability must to be 2 Rows 2 Cols Matrix" << std::endl;
    DebugStop();
  }
  
  Kinverse(0,0) =     1.0*Kab(1,1)/(Constant);
  Kinverse(0,1) = -   1.0*Kab(0,1)/(Constant);
  Kinverse(1,0) = -   1.0*Kab(1,0)/(Constant);
  Kinverse(1,1) =     1.0*Kab(0,0)/(Constant);
  return Kinverse;
    
}

REAL TPZFracData::G() const
{
#ifdef DEBUG
  if (fE < 0. || fnu < 0.) {
    PZError << "ERROR - Elasticity or poisson not set!" << std::endl;
    DebugStop();
  }
#endif
  return fE / (2.*(1+fnu));
}

void TPZFracData::SetNextTime()
{
  if(fTime + fDeltaT > fTtot)
  {
    fDeltaT = fTtot - fTime;
  }
  fTime += fDeltaT;
  std::cout << "\n----------------- Actual Time of Simulation = " << fTime << " -----------------" << std::endl;
}

void TPZFracData::SetPostProcessFileName(std::string &postProcessFileName)
{
  fpostProcessFileName = postProcessFileName;
}

std::string TPZFracData::PostProcessFileName()
{
  return fpostProcessFileName;
}

REAL TPZFracData::VlFtau(REAL pfrac, REAL tau) const
{
  const REAL Cl = this->Cl();
  REAL Pe = this->Pe();
  REAL Pref = this->Pref();
  REAL vsp = this->Vsp();
  
  REAL Pcalc = (pfrac - Pe)/Pref;
  if(Pcalc < 0.)
  {
    Pcalc = 0.;
  }
  
  REAL Clcorr = Cl;// * sqrt(Pcalc); AQUINATHAN
  REAL Vl = 2. * Clcorr * sqrt(tau) + vsp;
  
  return Vl;
}

REAL TPZFracData::FictitiousTime(REAL VlAcum, REAL pfrac) const
{
  REAL Cl = this->Cl();
  REAL Pe = this->Pe();
  REAL Pref = this->Pref();
  REAL vsp = this->Vsp();
  
  REAL tStar = 0.;
  if(VlAcum > vsp)
  {
    REAL Pcalc = (pfrac - Pe)/Pref;
    if(Pcalc < 0.)
    {
      Pcalc = 0.;
    }
    REAL Clcorr = Cl;// * sqrt(Pcalc); AQUINATHAN
    tStar = (VlAcum - vsp)*(VlAcum - vsp)/( (2. * Clcorr) * (2. * Clcorr) );
  }
  
  return tStar;
}

REAL TPZFracData::QlFVl(REAL VlAcum, REAL pfrac) const
{
  REAL deltaT = this->TimeStep();
  
  REAL tStar = FictitiousTime(VlAcum, pfrac);
  REAL Vlnext = VlFtau(pfrac, tStar + deltaT);
  REAL Ql = (Vlnext - VlAcum)/deltaT;
  
  return Ql;
  
}

REAL TPZFracData::dQlFVl(REAL VlAcum, REAL pfrac) const
{
  
  REAL deltaPfrac = fabs(pfrac/10000.);
  if(deltaPfrac < 1.E-10)
  {
    deltaPfrac = 1.E-10;
  }
  else if(deltaPfrac > 1.E-3)
  {
    deltaPfrac = 1.E-3;
  }
  
  REAL deltaT = this->TimeStep();
  /////////////////////////////////////////////////Ql maior
  REAL pfracUP = pfrac + deltaPfrac;
  REAL tStar1 = FictitiousTime(VlAcum, pfracUP);
  REAL Vlnext1 = VlFtau(pfracUP, tStar1 + deltaT);
  REAL Ql1 = (Vlnext1 - VlAcum )/deltaT;
  //...
  
  /////////////////////////////////////////////////Ql menor
  REAL pfracDOWN = pfrac - deltaPfrac;
  REAL tStar0 = FictitiousTime(VlAcum, pfracDOWN);
  REAL Vlnext0 = VlFtau(pfracDOWN, tStar0 + deltaT);
  REAL Ql0 = (Vlnext0 - VlAcum)/deltaT;
  //...
  
  REAL dQldpfrac = (Ql1-Ql0)/(2.*deltaPfrac);
  
  return dQldpfrac;

}

void TPZFracData::PrintDebugMapForMathematica(std::string filename)
{
  std::ofstream out(filename.c_str());
  std::map<REAL,REAL>::const_iterator it = fDebugMap.begin();
  
  out << "DebugMap = ";
  out << "{{" << it->first << "," << it->second << "}";
  it++;
  for (; it != fDebugMap.end(); it++) {
    out << ",{" << it->first << "," << it->second << "}";
  }
  out << "};" << std::endl;
  out << "ListPlot[DebugMap,Joined -> True,PlotMarkers -> Automatic]";
  
  out.close();
}

void TPZFracData::SetDwDp(){
  if (fnu <= 0. || fHf <= 0. || fE <= 0.){
    DebugStop(); // initialize those variables before using this method
  }
  fdwdp = 0.817 * (1. - this->Poisson()) * this->Hf() / this->G();
}

REAL TPZFracData::GetDwDp() const{
  return fdwdp;
}

REAL TPZFracData::GetW(REAL pfrac) const{
  return this->GetDwDp() * (pfrac - this->SigmaConf());
}