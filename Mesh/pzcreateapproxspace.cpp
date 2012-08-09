/**
 * @file
 * @brief Contains the implementation of the functions to creates computational elements.
 */

#include "pzcreateapproxspace.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzelctemp.h"
#include "pzcondensedcompel.h"
#include "pzinterpolationspace.h" 

#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"

#include "pzgeopoint.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzgeotriangle.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"

#include "pzcompelwithmem.h"

using namespace pzshape;

TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapePoint>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
    {
		TPZCompEl *result = new TPZIntelGen<TPZShapeLinear>(mesh,gel,index);
        return result;//new TPZCondensedCompel(result);
    }
	return NULL;
}
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeQuad>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeTriang>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeCube>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapePrism>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapePiram>(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeTetra>(mesh,gel,index);
	return NULL;
}


// with mem
TPZCompEl *CreatePointElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapePoint> >(mesh,gel,index) ;
	return NULL;
}
TPZCompEl *CreateLinearElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapeLinear> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateQuadElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapeQuad> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateTriangleElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZIntelGen<TPZShapeTriang> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateCubeElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapeCube> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreatePrismElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapePrism> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreatePyramElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapePiram> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *CreateTetraElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapeTetra> >(mesh,gel,index);
	return NULL;
}

/**
 * @brief Creates the computational elements, and the degree of freedom nodes
 */ 
/** Only element of material id in the set<int> will be created */
void TPZCreateApproximationSpace::BuildMesh(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs){
    cmesh.AutoBuild(MaterialIDs);
}

/** @brief Creates the computational elements, and the degree of freedom nodes */
void TPZCreateApproximationSpace::AutoBuild(TPZCompMesh &cmesh){
    cmesh.AutoBuild();
}

void TPZCreateApproximationSpace::CreateInterfaces(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs){
    for(int el = 0; el < cmesh.ElementVec().NElements(); el++)
	{
		TPZCompEl * compEl = cmesh.ElementVec()[el];
		if(!compEl) continue;
        int matid = compEl->Reference()->MaterialId();
        
        //checking if matid is in MaterialIDs
        std::set<int>::iterator it;
        it = MaterialIDs.find(matid);
        if(it!=MaterialIDs.end())
        {
            TPZInterpolationSpace * interpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
            if(!interpEl) continue;
            interpEl->CreateInterfaces();
        }
    }
    
}

void TPZCreateApproximationSpace::CreateInterfaces(TPZCompMesh &cmesh){
    for(int el = 0; el < cmesh.ElementVec().NElements(); el++)
	{
		TPZCompEl * compEl = cmesh.ElementVec()[el];
		if(!compEl) continue;
					
        TPZInterpolationSpace * interpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
        if(!interpEl) continue;
        interpEl->CreateInterfaces();
    }
    
}

void TPZCreateApproximationSpace::SetAllCreateFunctions(TPZCompEl &cel, TPZCompMesh *mesh){
	cel.SetCreateFunctions(mesh);
}



void TPZCreateApproximationSpace::SetAllCreateFunctionsDiscontinuous(){
	
    fp[EPoint] = TPZCompElDisc::CreateDisc;
    fp[EOned] = TPZCompElDisc::CreateDisc;
    fp[ETriangle] = TPZCompElDisc::CreateDisc;
    fp[EQuadrilateral] = TPZCompElDisc::CreateDisc;
    fp[ETetraedro] = TPZCompElDisc::CreateDisc;
    fp[EPiramide] = TPZCompElDisc::CreateDisc;
    fp[EPrisma] = TPZCompElDisc::CreateDisc;
    fp[ECube] = TPZCompElDisc::CreateDisc;
    
    /*
	pzgeom::TPZGeoPoint::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoLinear::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoQuad::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoTriangle::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoPrism::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoTetrahedra::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoPyramid::fp = TPZCompElDisc::CreateDisc;
	pzgeom::TPZGeoCube::fp = TPZCompElDisc::CreateDisc;
     */
}


void TPZCreateApproximationSpace::SetAllCreateFunctionsContinuous(){
    fp[EPoint] = CreatePointEl;
    fp[EOned] = CreateLinearEl;
    fp[ETriangle] = CreateTriangleEl;
    fp[EQuadrilateral] = CreateQuadEl;
    fp[ETetraedro] = CreateTetraEl;
    fp[EPiramide] = CreatePyramEl;
    fp[EPrisma] = CreatePrismEl;
    fp[ECube] = CreateCubeEl;
    
    /*
	pzgeom::TPZGeoPoint::fp =  CreatePointEl;
	pzgeom::TPZGeoLinear::fp =  CreateLinearEl;
	pzgeom::TPZGeoQuad::fp = CreateQuadEl;
	pzgeom::TPZGeoTriangle::fp =  CreateTriangleEl;
	pzgeom::TPZGeoPrism::fp = CreatePrismEl;
	pzgeom::TPZGeoTetrahedra::fp = CreateTetraEl;
	pzgeom::TPZGeoPyramid::fp = CreatePyramEl;
	pzgeom::TPZGeoCube::fp = CreateCubeEl;
	*/
	
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsContinuousWithMem(int dimension)
{
	// These ones will be always continuous for viscoelasticity
	fp[EPoint] = CreatePointEl;
	fp[EOned] = CreateLinearEl;
	fp[ETriangle] = CreateTriangleEl;
	fp[EQuadrilateral] = CreateQuadEl;
	
	if (dimension < 3 && dimension > 0) 	// If the mesh is dimension < 3, it doesnt matter which 3d elements create for viscoelasticity
	{
    fp[ETetraedro] = CreateTetraEl;
    fp[EPiramide] = CreatePyramEl;
    fp[EPrisma] = CreatePrismEl;
    fp[ECube] = CreateCubeEl;
	}
	else if (dimension == 3) // Otherwise it creates 3d elements with memory
	{
		fp[ETetraedro] = CreateTetraElWithMem;
    fp[EPiramide] = CreatePyramElWithMem;
    fp[EPrisma] = CreatePrismElWithMem;
    fp[ECube] = CreateCubeElWithMem;
	}
	else
	{
		PZError << "Invalid dimension value in TPZCreateApproximationSpace::SetAllCreateFunctionsContinuousWithMem(int dimension) \n";
		DebugStop();
	}

    /*
	pzgeom::TPZGeoPoint::fp =  CreatePointElWithMem;
	pzgeom::TPZGeoLinear::fp =  CreateLinearElWithMem;
	pzgeom::TPZGeoQuad::fp = CreateQuadElWithMem;
	pzgeom::TPZGeoTriangle::fp =  CreateTriangleElWithMem;
	pzgeom::TPZGeoPrism::fp = CreatePrismElWithMem;
	pzgeom::TPZGeoTetrahedra::fp = CreateTetraElWithMem;
	pzgeom::TPZGeoPyramid::fp = CreatePyramElWithMem;
	pzgeom::TPZGeoCube::fp = CreateCubeElWithMem;
     */
}


#include "pzelchdiv.h"

void TPZCreateApproximationSpace::SetAllCreateFunctionsHDiv(){
	
    fp[EPoint] = CreateHDivPointEl;
    fp[EOned] = CreateHDivLinearEl;
    fp[ETriangle] = CreateHDivTriangleEl;
    fp[EQuadrilateral] = CreateHDivQuadEl;
    fp[ETetraedro] = CreateHDivTetraEl;
    fp[EPiramide] = CreateHDivPyramEl;
    fp[EPrisma] = CreateHDivPrismEl;
    fp[ECube] = CreateHDivCubeEl;
    
    /*
    pzgeom::TPZGeoPoint::fp =  CreateHDivPointEl;
	pzgeom::TPZGeoLinear::fp =  CreateHDivLinearEl;
	pzgeom::TPZGeoQuad::fp = CreateHDivQuadEl;
	pzgeom::TPZGeoTriangle::fp =  CreateHDivTriangleEl;
	pzgeom::TPZGeoPrism::fp = CreateHDivPrismEl;
	pzgeom::TPZGeoTetrahedra::fp = CreateHDivTetraEl;
	pzgeom::TPZGeoPyramid::fp = CreateHDivPyramEl;
	pzgeom::TPZGeoCube::fp = CreateHDivCubeEl;
     */
}

#ifndef STATE_COMPLEX
#include "pzhdivpressure.h"

void TPZCreateApproximationSpace::SetAllCreateFunctionsHDivPressure(){
		
    fp[EPoint] = CreateHDivPressurePointEl;
    fp[EOned] = CreateHDivPressureLinearEl;
    fp[ETriangle] = CreateHDivPressureTriangleEl;
    fp[EQuadrilateral] = CreateHDivPressureQuadEl;
    fp[ETetraedro] = CreateHDivPressureTetraEl;
    fp[EPiramide] = CreateHDivPressurePyramEl;
    fp[EPrisma] = CreateHDivPressurePrismEl;
    fp[ECube] = CreateHDivPressureCubeEl;
    
 }
#endif


#include "pzreferredcompel.h"
#include "pzelctemp.h"
void TPZCreateApproximationSpace::SetAllCreateFunctionsDiscontinuousReferred(){
	
    fp[EPoint] = CreateReferredDisc;
    fp[EOned] = CreateReferredDisc;
    fp[ETriangle] = CreateReferredDisc;
    fp[EQuadrilateral] = CreateReferredDisc;
    fp[ETetraedro] = CreateReferredDisc;
    fp[EPiramide] = CreateReferredDisc;
    fp[EPrisma] = CreateReferredDisc;
    fp[ECube] = CreateReferredDisc;
    
    /*
	pzgeom::TPZGeoPoint::fp =  CreateReferredDisc;
	pzgeom::TPZGeoLinear::fp =  CreateReferredDisc;
	pzgeom::TPZGeoQuad::fp = CreateReferredDisc;
	pzgeom::TPZGeoTriangle::fp =  CreateReferredDisc;
	pzgeom::TPZGeoPrism::fp = CreateReferredDisc;
	pzgeom::TPZGeoTetrahedra::fp = CreateReferredDisc;
	pzgeom::TPZGeoPyramid::fp = CreateReferredDisc;
	pzgeom::TPZGeoCube::fp = CreateReferredDisc;
     */
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsContinuousReferred(){
	
    fp[EPoint] = CreateReferredPointEl;
    fp[EOned] = CreateReferredLinearEl;
    fp[ETriangle] = CreateReferredTriangleEl;
    fp[EQuadrilateral] = CreateReferredQuadEl;
    fp[ETetraedro] = CreateReferredTetraEl;
    fp[EPiramide] = CreateReferredPyramEl;
    fp[EPrisma] = CreateReferredPrismEl;
    fp[ECube] = CreateReferredCubeEl;
    
    /*
	pzgeom::TPZGeoPoint::fp =  CreateReferredPointEl;
	pzgeom::TPZGeoLinear::fp =  CreateReferredLinearEl;
	pzgeom::TPZGeoQuad::fp = CreateReferredQuadEl;
	pzgeom::TPZGeoTriangle::fp =  CreateReferredTriangleEl;
	pzgeom::TPZGeoPrism::fp = CreateReferredPrismEl;
	pzgeom::TPZGeoTetrahedra::fp = CreateReferredTetraEl;
	pzgeom::TPZGeoPyramid::fp = CreateReferredPyramEl;
	pzgeom::TPZGeoCube::fp = CreateReferredCubeEl;
     */
	
}

#include "pzmultiphysicscompel.h"
void TPZCreateApproximationSpace::SetAllCreateFunctionsMultiphysicElem(){
	
    fp[EPoint] = CreateMultiphysicsPointEl;
    fp[EOned] = CreateMultiphysicsLinearEl;
    fp[ETriangle] = CreateMultiphysicsTriangleEl;
    fp[EQuadrilateral] = CreateMultiphysicsQuadEl;
    fp[ETetraedro] = CreateMultiphysicsTetraEl;
    fp[EPiramide] = CreateMultiphysicsPyramEl;
    fp[EPrisma] = CreateMultiphysicsPrismEl;
    fp[ECube] = CreateMultiphysicsCubeEl;
    
    /*
	pzgeom::TPZGeoPoint::fp =  CreateMultiphysicsPointEl;
	pzgeom::TPZGeoLinear::fp =  CreateMultiphysicsLinearEl;
	pzgeom::TPZGeoTriangle::fp =  CreateMultiphysicsTriangleEl;
	pzgeom::TPZGeoQuad::fp = CreateMultiphysicsQuadEl;
	pzgeom::TPZGeoCube::fp = CreateMultiphysicsCubeEl;
	pzgeom::TPZGeoPrism::fp = CreateMultiphysicsPrismEl;
	pzgeom::TPZGeoTetrahedra::fp = CreateMultiphysicsTetraEl;
	pzgeom::TPZGeoPyramid::fp = CreateMultiphysicsPyramEl;
     */
}

/*
 * @brief Create a computational element using the function pointer for the topology
 */
TPZCompEl *TPZCreateApproximationSpace::CreateCompEl(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
    switch (gel->Type()) {
        case EPoint:
            return fp[EPoint](gel,mesh,index);
            break;
        case EOned:
            return fp[EOned](gel,mesh,index);
            break;
        case EQuadrilateral:
            return fp[EQuadrilateral](gel,mesh,index);
            break;
        case ETriangle:
            return fp[ETriangle](gel,mesh,index);
            break;
        case EPiramide:
            return fp[EPiramide](gel,mesh,index);
            break;
        case EPrisma:
            return fp[EPrisma](gel,mesh,index);
            break;
        case ETetraedro:
            return fp[ETetraedro](gel,mesh,index);
            break;
        case ECube:
            return fp[ECube](gel,mesh,index);
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

/**
 * @brief Set custom function pointers
 */

void TPZCreateApproximationSpace::SetCreateFunctions(TPZVec<TCreateFunction> &createfuncs)
{
    fp[EPoint] = createfuncs[EPoint];
    fp[EOned] = createfuncs[EOned];
    fp[ETriangle] = createfuncs[ETriangle];
    fp[EQuadrilateral] = createfuncs[EQuadrilateral];
    fp[ETetraedro] = createfuncs[ETetraedro];
    fp[EPiramide] = createfuncs[EPiramide];
    fp[EPrisma] = createfuncs[EPrisma];
    fp[ECube] = createfuncs[ECube];
    
}

#include "pzcondensedcompel.h"

/**
 * @brief Encapsulate the elements in condensed computational elements
 */
void TPZCreateApproximationSpace::CondenseLocalEquations(TPZCompMesh &cmesh)
{
    int nel = cmesh.NElements();
    int iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh.ElementVec()[iel];
        if(!cel) {
            continue;
        }
        new TPZCondensedCompEl(cel);
    }

}

/**
 * @brief Undo the encapsulate elements
 */
void TPZCreateApproximationSpace::UndoCondenseLocalEquations(TPZCompMesh &cmesh)
{
    int nel = cmesh.NElements();
    int iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh.ElementVec()[iel];
        TPZCondensedCompEl *condel = dynamic_cast<TPZCondensedCompEl *>(cel);
        if(!condel) {
            continue;
        }
        condel->Unwrap();
    }    
}

/**
 * @brief transform in low order Raviar Tomas
 */
void TPZCreateApproximationSpace::MakeRaviartThomas(TPZCompMesh &cmesh)
{
    int numcell = cmesh.NElements();
    int el;
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            continue;
        }
        intel->SetPreferredOrder(1);
    }
    cmesh.ExpandSolution();
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        int geldim = gel->Dimension();
        int is;
        int nsides = gel->NSides();
        for (is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != geldim-1) {
                continue;
            }
            int nsconnects = intel->NSideConnects(is);
            // only interested in HDiv elements
            if (nsconnects != 1) {
                continue;
            }
            int cindex = intel->SideConnectIndex(0, is);
            TPZConnect &c = intel->Connect(intel->SideConnectLocId(0,is));
            if (c.HasDependency()) {
                continue;
            }
            int nshape = 1;
            int nstate = 1;
            int order = 0;
            int cindex2 = cmesh.AllocateNewConnect(nshape, nstate, order);
//            TPZConnect &c2 = cmesh.ConnectVec()[cindex];
            TPZFNMatrix<2> depmat(2,1,1.);
            c.AddDependency(cindex, cindex2, depmat, 0, 0, 2, 1);
        }
    }
    cmesh.ExpandSolution();
}

/**
 * @brief transform in low order Raviar Tomas
 */
void TPZCreateApproximationSpace::UndoMakeRaviartThomas(TPZCompMesh &cmesh)
{
    int numcell = cmesh.NElements();
    int el;
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        int geldim = gel->Dimension();
        int is;
        int nsides = gel->NSides();
        for (is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != geldim-1) {
                continue;
            }
            int nsconnects = intel->NSideConnects(is);
            // only interested in HDiv elements
            if (nsconnects != 1) {
                continue;
            }
//            int cindex = intel->SideConnectIndex(0, is);
            TPZConnect &c = intel->Connect(intel->SideConnectLocId(0,is));
            if (c.HasDependency()) {
                c.RemoveDepend();
            }
        }
    }
    cmesh.ExpandSolution();
    cmesh.CleanUpUnconnectedNodes();
}

/** @brief Create interface elements between the computational elements */
void TPZCreateApproximationSpace::CreateInterfaceElements(TPZCompMesh *mesh, bool betweencontinuous, bool multiphysics)
{
    TPZChunkVector<TPZCompEl *> compelvec = mesh->ElementVec();
    int nel = compelvec.NElements();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = compelvec[el];
        if (!cel) {
            continue;
        }
        if(!multiphysics)
        {
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (intel) {
                intel->CreateInterfaces(betweencontinuous);
            }
        }
        else {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mfcel) {
                continue;
            }
            mfcel->CreateInterfaces();
        }
        
    }
}
