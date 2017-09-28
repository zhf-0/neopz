/*
 *  pzgradient.cpp
 *  PZ
 *
 *  Created by Agnaldo on 14/01/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string>

#include "pzgradient.h"

TPZGradient::TPZGradient(): TPZRegisterClassId(&TPZGradient::ClassId),
                            TPZFunction<STATE>(){
    
    fCenter.Resize(3,0.);
    fGradient.Resize(3);
    fUc=0.;
    falphaK = 1.;
}

TPZGradient::TPZGradient(const TPZGradient &cp): TPZRegisterClassId(&TPZGradient::ClassId),TPZFunction<STATE>(cp), fCenter(cp.fCenter), fGradient(cp.fGradient){
    
    fUc = cp.fUc;
    falphaK = cp.falphaK;
}

int TPZGradient::ClassId() {
    //CLASSIDFRANreturn TPZFunction<STATE>::ClassId()^Hash("TPZGradient");
    return 666;
}

