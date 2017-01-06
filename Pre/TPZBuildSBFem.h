//
//  TPZBuildSBFem.hpp
//  PZ
//
//  Created by Philippe Devloo on 06/01/17.
//
//

#ifndef TPZBuildSBFem_hpp
#define TPZBuildSBFem_hpp

#include <stdio.h>
#include "pzmanvector.h"
#include <map>

class TPZBuildSBFem
{
    /// The volumetric elements with Group Id will spawn SBFemVolume elements with MatId
    std::map<int,int> fMatIdTranslation;
    
    /// Material Id associated with the skeleton elements
    int fSkeletonMatId;
    
    /// partition to which each sbfem element belongs
    
    /// partition to which each group element belongs
};

#endif /* TPZBuildSBFem_hpp */
