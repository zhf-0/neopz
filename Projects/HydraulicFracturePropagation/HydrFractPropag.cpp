#include <iostream>

#include "TPZRefPatternDataBase.h"
#include "TPZPlaneFractureKernel.h"

#include "TPZTimer.h"

using namespace std;



#define PDoreTest1
int main(int argc, char * const argv[])
{
    std::cout << "\e";
    TPZTimer readRef("ReadingRefPatterns");
    readRef.start();
    
    //#define writeAgain
    #ifdef writeAgain
    gRefDBase.InitializeRefPatterns();
    #else
    std::ifstream inRefP("RefPatternsUsed.txt");
    gRefDBase.ReadRefPatternDBase("RefPatternsUsed.txt");
    #endif
    
    readRef.stop();
    std::cout << "DeltaT leitura refpatterns = " << readRef.seconds() << " s" << std::endl;
    
#ifdef PDoreTest1

    //Transient data
    REAL Ttot = 25. * 60.; /** em segundos */
    REAL maxDeltaT = 60.; /** em segundos */
    int nTimes = 1; /** quantidade de divisao do maxDeltaT para definir minDeltaT (minDeltaT = maxDeltaT/nTimes) */
    globTimeControl.SetTimeControl(Ttot, maxDeltaT, nTimes);
    
    //Geometry data
    REAL lengthX = 100.;
    REAL lengthY = 60.;
    REAL Lmax = 2.0;
    
    REAL bulletTVDIni = 2100.;
    REAL bulletTVDFin = 2120.;
    int nstripes = 2;
    
    //Material data
    TPZVec<TPZLayerProperties> layerVec(3);
  
    REAL Young0 = 4.2747495136E10;
    REAL Poisson0 = 0.25;
    REAL SigMax0  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin0  = -35627279.293;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf0 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL Young1 = 4.1368543680E10;
    REAL Poisson1 = 0.15;
    REAL SigMax1  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin1  = -34528254.983;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf1 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL Young2 = 4.2747495136E10;
    REAL Poisson2 = 0.25;
    REAL SigMax2  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin2  = -37324079.06;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf2 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL TVDi0 = 2070.;
    REAL TVDf0 = 2100.;
    REAL TVDi1 = TVDi0;
    REAL TVDf1 = 2120.;
    REAL TVDi2 = TVDi1;
    REAL TVDf2 = 2160.;
    
    //KIc
    REAL KIc0 = 1.19884E6;//AQUICAJU
    REAL KIc1 = 1.09884E6;
    REAL KIc2 = 1.29884E6;//AQUICAJU
    
    //Leakoff
    REAL Cl0 = 0.;
    REAL Pe0 = 0.;//Sempre positivo
    REAL gradPref0 = 1.;
    REAL vsp0 = 0.;
    
    REAL Cl1 = 0.00019674755398733676;
    REAL Pe1 = 0.;//Sempre positivo
    REAL gradPref1 = 1.;
    REAL vsp1 = 0.;
    
    REAL Cl2 = 0.;
    REAL Pe2 = 0.;//Sempre positivo
    REAL gradPref2 = 1.;
    REAL vsp2 = 0.;
    
    layerVec[0] = TPZLayerProperties(Young0, Poisson0, SigMax0, SigMin0, SigConf0, TVDi0, TVDf0, KIc0, Cl0, Pe0, gradPref0, vsp0);
    layerVec[1] = TPZLayerProperties(Young1, Poisson1, SigMax1, SigMin1, SigConf1, TVDi1, TVDf1, KIc1, Cl1, Pe1, gradPref1, vsp1);
    layerVec[2] = TPZLayerProperties(Young2, Poisson2, SigMax2, SigMin2, SigConf2, TVDi2, TVDf2, KIc2, Cl2, Pe2, gradPref2, vsp2);
    
    //Fluid injection data
    REAL QinjWell = 1.*(-0.0533333333333);//m3/s (1.* pois os 80 bpm jah eh no poco e nao 1 wing)
    REAL visc = 200.02E-3;//N.s/m2
    
    //J-Integral data
    REAL Jradius = 1.5;
    
#else
    
    //Transient data
    REAL Ttot = 200. * 60.; /** em segundos */
    REAL maxDeltaT = Ttot/8.; /** em segundos */
    int nTimes = 10; /** quantidade de divisao do maxDeltaT para definir minDeltaT (minDeltaT = maxDeltaT/nTimes) */
    globTimeControl.SetTimeControl(Ttot, maxDeltaT, nTimes);
    
    //Geometry data
    REAL lengthX = 300.;
    REAL lengthY = 20.;
    REAL Lmax = 10.0;
    
    REAL bulletTVDIni = 2795.;
    REAL bulletTVDFin = 2847.;
    int nstripes = 1;
    
    //Material data
    TPZVec<TPZLayerProperties> layerVec(3);
    
    REAL Young0 = 4.48159205E10;
    REAL Poisson0 = 0.30;
    REAL SigMax0  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin0  = -49297512.55;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf0 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL Young1 = 5.8654345E10;
    REAL Poisson1 = 0.21;
    REAL SigMax1  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin1  = -39300114.9;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf1 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL Young2 = 3.79211635E10;
    REAL Poisson2 = 0.29;
    REAL SigMax2  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin2  = -50676463.95;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf2 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL TVDi0 = 2740.;
    REAL TVDf0 = 2795.;
    REAL TVDi1 = TVDi0;
    REAL TVDf1 = 2847.;
    REAL TVDi2 = TVDi1;
    REAL TVDf2 = 2941.;
    
    //KIc
    REAL KIc0 = 2.19769E6;
    REAL KIc1 = 2.19769E6;
    REAL KIc2 = 2.19769E6;
    
    //Leakoff
    REAL Cl0 = 9.83742E-6;
    REAL Pe0 = 0.;//Sempre positivo
    REAL gradPref0 = 1.;
    REAL vsp0 = 0.;
    
    REAL Cl1 = 9.83742E-6;
    REAL Pe1 = 0.;//Sempre positivo
    REAL gradPref1 = 1.;
    REAL vsp1 = 0.;
    
    REAL Cl2 = 9.83742E-6;
    REAL Pe2 = 0.;//Sempre positivo
    REAL gradPref2 = 1.;
    REAL vsp2 = 0.;
    
    layerVec[0] = TPZLayerProperties(Young0, Poisson0, SigMax0, SigMin0, SigConf0, TVDi0, TVDf0, KIc0, Cl0, Pe0, gradPref0, vsp0);
    layerVec[1] = TPZLayerProperties(Young1, Poisson1, SigMax1, SigMin1, SigConf1, TVDi1, TVDf1, KIc1, Cl1, Pe1, gradPref1, vsp1);
    layerVec[2] = TPZLayerProperties(Young2, Poisson2, SigMax2, SigMin2, SigConf2, TVDi2, TVDf2, KIc2, Cl2, Pe2, gradPref2, vsp2);
    
    //Fluid injection data
    REAL QinjWell = 2.*(-0.1324894);//m3/s
    REAL visc = 200.0E-3;//N.s/m2
    
    //J-Integral data
    REAL Jradius = 1.5;
    
#endif

    //Simulation p-order data
    int porder = 1;
    
    TPZPlaneFractureKernel * plfrac = new TPZPlaneFractureKernel(layerVec, bulletTVDIni, bulletTVDFin, lengthX, lengthY, Lmax, nstripes,
                                                                 QinjWell, visc,
                                                                 Jradius,
                                                                 porder);

    plfrac->Run();
    
//    std::ofstream outRefP("RefPatternsUsed.txt");
//    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}
