//$Id: poroelastoplasticid.h,v 1.4 2009-10-04 05:44:22 erick Exp $
 
// pz id for the classes appended by this project

#include "TPZPlasticStepID.h"

#ifndef POROELASTOPLASTICID_H
#define POROELASTOPLASTICID_H

const int NUMPLASTICMODELS = BASEPLASTICMODEL_ID - LASTPLASTIC_ID + 1; // 399-10+1 = 390

const int TPZELASTOPLASTICMEM_ID      = LASTPLASTIC_ID - 1;                          // 9
const int TPZPOROELASTOPLASTICMEM_ID  = TPZELASTOPLASTICMEM_ID - 1;                  // 8
const int TPZMATWITHMEM_ID            = TPZPOROELASTOPLASTICMEM_ID -1;               // 7
const int TPZMATELASTOPLASTIC_ID      = TPZMATWITHMEM_ID - 1;                        // 6
const int TPZMATPOROELASTOPLASTIC_ID  = TPZMATELASTOPLASTIC_ID - NUMPLASTICMODELS;   // 6-390=-384
//const int TPZPOSTPROCMAT_ID           = TPZMATPOROELASTOPLASTIC_ID - NUMPLASTICMODELS;
const int TPZMATELASTOPLASTIC2D_ID    = TPZMATWITHMEM_ID - 1;                        // 6

#endif
