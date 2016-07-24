//
//  TPZSBFemVolume.cpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#include "TPZSBFemVolume.h"
#include "pzintel.h"
#include "pzmaterial.h"
#include "pzelmat.h"


TPZSBFemVolume::TPZSBFemVolume(TPZCompMesh &mesh, TPZGeoEl *gel,long &index) : TPZCompEl(mesh,gel,index), fElementGroupIndex(-1), fSkeleton(-1)
{
    
}



/// Compute the K matrices
void TPZSBFemVolume::ComputeKMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2)
{
    // do all the computations here
    
    TPZElementMatrix efmat(Mesh(),TPZElementMatrix::EF);
    
    TPZGeoEl *Ref2D = Reference();
    TPZGeoMesh *gmesh = Ref2D->Mesh();
    
    TPZCompMesh *cmesh = Mesh();
    
    TPZInterpolatedElement *CSkeleton = dynamic_cast<TPZInterpolatedElement *>(cmesh->Element(fSkeleton));
    
    CSkeleton->InitializeElementMatrix(E0, efmat);
    CSkeleton->InitializeElementMatrix(E1, efmat);
    CSkeleton->InitializeElementMatrix(E2, efmat);
    
    TPZGeoEl *Ref1D = CSkeleton->Reference();
    
    int matid = Ref2D->MaterialId();
    int dimension = Ref2D->Dimension();
    
    // find the first face side
    int nsides = Ref2D->NSides();
    int is;
    for (is=nsides-1; is>0; is--) {
        if (Ref2D->SideDimension(is) < dimension-1) {
            break;
        }
    }
    int faceside = is+1;
    
    TPZGeoElSide thisside(Ref2D,faceside);
    
    
    TPZMaterial *mat2d = cmesh->FindMaterial(matid);
    
    if(!mat2d) DebugStop();
    
    TPZGeoElSide SkeletonSide(Ref1D,Ref1D->NSides()-1);
    
    TPZTransform tr(2, 1);
    tr = SkeletonSide.NeighbourSideTransform(thisside);
    TPZTransform t2 = Ref2D->SideToSideTransform(thisside.Side(), Ref2D->NSides()-1);
    tr = t2.Multiply(tr);
    // create a one-d integration rule
    TPZIntPoints &intpoints = CSkeleton->GetIntegrationRule();
    
    TPZMaterialData data1d;
    TPZMaterialData data2d;
    CSkeleton->InitMaterialData(data1d);
    CSkeleton->InitMaterialData(data2d);
    int nshape = data2d.phi.Rows();
    data2d.phi.Redim(nshape*2, 1);
    data2d.dphi.Redim(2, 2*nshape);
    data2d.dphix.Redim(2, 2*nshape);
    data2d.dsol[0].Redim(2,2);
    
    TPZFNMatrix<200,STATE> ek(nshape*4,nshape*4,0.), ef(nshape*4,1,0.);
    int npoint = intpoints.NPoints();
    for (int ip = 0; ip<npoint; ip++)
    {
        TPZManVector<REAL,3> xi(1), xiquad(2);
        REAL weight;
        intpoints.Point(ip, xi, weight);
        tr.Apply(xi, xiquad);
        TPZFNMatrix<9,REAL> jacobian(1,1),axes(1,3),jacinv(1,1);
        REAL detjac;
        Ref1D->Jacobian(xi,jacobian,axes,detjac,jacinv);
        Ref2D->Jacobian(xiquad, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
        CSkeleton->ComputeRequiredData(data1d, xi);
        ExtendShapeFunctions(data1d,data2d);
        
        weight *= fabs(data2d.detjac)*2.;
        // compute the contributions to K11 K12 and K22
        mat2d->Contribute(data2d,weight,ek,ef);
    }
    for (int i=0; i<2*nshape; i++) {
        for (int j=0; j<2*nshape; j++) {
            E0.fMat(i,j) = ek(i,j);
            E1.fMat(j,i) = ek(i,j+2*nshape);
            E2.fMat(i,j) = ek(i+2*nshape,j+2*nshape);
        }
    }
}

/// extend the border shape functions for SBFem computations
void TPZSBFemVolume::ExtendShapeFunctions(TPZMaterialData &data1d, TPZMaterialData &data2d)
{
    int dim = Reference()->Dimension();
    long nshape = data2d.phi.Rows()/2;
    for (int ish=0; ish<nshape; ish++) {
        data2d.phi(ish+nshape,0) = data1d.phi(ish,0);
        for (int d=0; d<dim-1; d++) {
            data2d.dphi(d,ish+nshape) = data1d.dphi(d,ish);
            data2d.dphi(d,ish) = 0.;
        }
        data2d.dphi(dim-1,ish) = -data1d.phi(ish)/2.;
        data2d.dphi(dim-1,ish+nshape) = 0.;
    }
    TPZInterpolationSpace::Convert2Axes(data2d.dphi, data2d.jacinv, data2d.dphix);

}

TPZCompEl * CreateSBFemCompEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index)
{
    return new TPZSBFemVolume(mesh,gel,index);
}

