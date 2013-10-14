/*
 *  TPZCompElMHM.h
 *  PZ
 *
 *  Created by Frederico on 30/01/13.
 *  Copyright 2013 LNCC. All rights reserved.
 *
 */

#ifndef TPZCOMELMHMH
#define TPZCOMELMHMH

#include <iostream>
#include <vector>

#include "pzelctemp.h"
#include "pzmanvector.h"
#include "pzvec.h"
#include "pzcompel.h"

#include <pzshapelinear.h>
#include <pzshapequad.h>
#include <pzshapetriang.h>
#include <pzshapetetra.h>
#include <pzshapecube.h>
#include <pzshapepiram.h>
#include <pzshapepoint.h>
#include <pzshapeprism.h>

#include <pzgeopoint.h>
#include <pzgeotriangle.h>
#include <TPZGeoLinear.h>
#include <pzgeoquad.h>
#include <pzgeotetrahedra.h>
#include <pzgeoprism.h>
#include <pzgeopyramid.h>
#include <TPZGeoCube.h>

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"

#include "LocalMHM.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"

template <class TSHAPE, class TGEO>
class TPZCompElMHM : public TPZIntelGen<TSHAPE>
{
private:
    bool LocalBasisComputed;
    int LocalMaterial;
    int LocalInterpolationOrder;
    int FaceInterpolationOrder;
    int LocalRefinement;
    int NumberOfLocalProblems;
    std::vector<TPZCompEl*> compElFaces;
    std::vector<TPZGeoEl*> geoElFaces;

    TPZVec<TPZCompMesh*> LocalCompMesh;
    TPZGeoMesh* LocalGeoMesh;

    void RemoveInterfaces(TPZCompMesh*);
    void ConfigureLangrangeCompMesh(TPZCompMesh*, TPZVec<TPZCompMesh*> );
    void ResizeAndSolveSystem(TPZAnalysis &);
    bool IsXOnTheFace( TPZGeoEl* geo, TPZVec<REAL>& X );

public:
	TPZCompElMHM();
	TPZCompElMHM(TPZCompMesh &mesh, TPZGeoEl *ref, int &index);
	virtual ~TPZCompElMHM();

    virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix &jacobian, TPZFMatrix &axes,
							  REAL &detjac, TPZFMatrix &jacinv, TPZFMatrix &phi, TPZFMatrix &dphix);

    virtual void ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi);

    virtual void ComputeSolution(TPZVec<REAL> &X, TPZSolVec &sol, TPZGradSolVec &dsol, TPZFMatrix& axes);
};

TPZCompEl * CreateMHMPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl * CreateMHMLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl * CreateMHMTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl * CreateMHMQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl * CreateMHMCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl * CreateMHMPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl * CreateMHMTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl * CreateMHMPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif  // TPZCOMELMHMH