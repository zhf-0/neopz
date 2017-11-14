//
//  TPZBuildSBFem.cpp
//  PZ
//
//  Created by Philippe Devloo on 06/01/17.
//
//

#include "TPZBuildSBFem.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"
#include "pzcompel.h"
#include "tpzgeoblend.h"

#include "tpzgeoelrefpattern.h"


/// standard configuration means each element is a partition and a center node is created
void TPZBuildSBFem::StandardConfiguration()
{
    TPZVec<long> elindices(fGMesh->NElements());
    long nel = elindices.size();
    for (long el=0; el<nel; el++) {
        elindices[el] = el;
    }
    StandardConfiguration(elindices);
}

/// standard configuration means each element is a partition and a center node is created for the indicated elements
void TPZBuildSBFem::StandardConfiguration(TPZVec<long> &elementindices)
{
    CreateElementCenterNodes(elementindices);
    AddSkeletonElements();
}

/// build element groups according to the id of the scaling centers
void TPZBuildSBFem::Configure(TPZVec<long> &scalingcenters)
{
    std::map<long,long> nodetogroup;
//    int maxpartition = 4;
    long count = 0;
    for (long el=0; el<scalingcenters.size(); el++) {
        if (scalingcenters[el] == -1) {
            continue;
        }
        if (nodetogroup.find(scalingcenters[el]) == nodetogroup.end()) {
            nodetogroup[scalingcenters[el]] = count++;
//            if (count >= maxpartition) {
//                break;
//            }
        }
    }
    if(fPartitionCenterNode.size())
    {
        DebugStop();
    }
    int dim = fGMesh->Dimension();
    fPartitionCenterNode.resize(nodetogroup.size());
    for (std::map<long,long>::iterator it = nodetogroup.begin(); it != nodetogroup.end(); it++) {
        fPartitionCenterNode[it->second] = it->first;
    }
    long nel = fGMesh->NElements();
    for (long el=0; el < nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (gel->Dimension() != dim && scalingcenters[el] != -1) {
            DebugStop();
        }
        else if(gel->Dimension() != dim)
        {
            continue;
        }
        if (scalingcenters[el] == -1) {
            continue;
        }
        if (nodetogroup.find(scalingcenters[el]) == nodetogroup.end()) {
            DebugStop();
        }
        long partition = nodetogroup[scalingcenters[el]];
//        if(partition < maxpartition)
        {
            fElementPartition[el] = partition;
        }
    }
    AddSkeletonElements();

}


/// add a partition manually
void TPZBuildSBFem::AddPartition(TPZVec<long> &elindices, long centernodeindex)
{
    long npart = fPartitionCenterNode.size();
    fPartitionCenterNode.resize(npart+1);
    if (fGMesh->NodeVec().NElements() <= centernodeindex) {
        DebugStop();
    }
    fPartitionCenterNode[npart] = centernodeindex;
    long nel = elindices.size();
    for (long el=0; el<nel; el++) {
        long elindex = elindices[el];
        if (fElementPartition[elindex] != -1) {
            DebugStop();
        }
        fElementPartition[elindex] = npart;
    }
}

/// add the sbfem elements to the computational mesh, the material should exist in cmesh
void TPZBuildSBFem::BuildComputationMesh(TPZCompMesh &cmesh)
{
    // create the lower dimensional mesh
    std::set<int> matids;
    int dim = cmesh.Dimension();
    TPZGeoMesh *gmesh = cmesh.Reference();
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        if (gel->Dimension() < dim) {
            matids.insert(gel->MaterialId());
        }
    }
    // create the boundary elements
    cmesh.ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh.AutoBuild(matids);
    CreateVolumetricElements(cmesh);
    CreateElementGroups(cmesh);
    
}


/// create the geometric skeleton elements
void TPZBuildSBFem::AddSkeletonElements()
{
    // create a lower dimension element on each boundary
    int dim = fGMesh->Dimension();
    
    long nel = fGMesh->NElements();
    
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim) {
            continue;
        }
        // the element doesnt belong to any partition, do not create a skeleton element
        if (fElementPartition[el] == -1) {
            continue;
        }
        long elpartition = fElementPartition[el];
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide thisside(gel,is);
            // we do not create skeleton elements on the boundary of the domain
            TPZGeoElSide neighbour = thisside.Neighbour();
            if (neighbour == thisside) {
                continue;
            }
            long neighbourelpartition = -1;
            // look for a neighbour with mesh dimension
            while(neighbour != thisside && neighbour.Element()->Dimension() != dim)
            {
                neighbour = neighbour.Neighbour();
            }
            // we found one!
            if(neighbour != thisside)
            {
                neighbourelpartition = fElementPartition[neighbour.Element()->Index()];
            }
            neighbour = thisside.Neighbour();
            
            // look for a neighbour of dimension dim-1
            while (neighbour != thisside && neighbour.Element()->Dimension() != dim-1) {
                neighbour = neighbour.Neighbour();
            }
            // if we didnt find a lower dimension element and the neighbouring elements of same dimension belong to a different partition
            if (thisside == neighbour && elpartition != neighbourelpartition) {
                gel->CreateBCGeoEl(is,fSkeletonMatId);
            }
        }
    }

}

/// create a geometric node at the center of each partition
void TPZBuildSBFem::CreateElementCenterNodes(TPZVec<long> &elindices)
{
    if(fPartitionCenterNode.size())
    {
        DebugStop();
    }
    int dim = fGMesh->Dimension();
    fPartitionCenterNode.resize(elindices.size());
    long count = 0;
    for (long el=0; el<elindices.size(); el++) {
        TPZGeoEl *gel = fGMesh->Element(elindices[el]);
        if (gel->Dimension() != dim) {
            continue;
        }
        int nsides = gel->NSides();
        TPZManVector<REAL,3> xicenter(dim),xcenter(3);
        gel->CenterPoint(nsides-1,xicenter);
        gel->X(xicenter,xcenter);
        long middlenode = fGMesh->NodeVec().AllocateNewElement();
        fGMesh->NodeVec()[middlenode].Initialize(xcenter,fGMesh);
        fPartitionCenterNode[count] = middlenode;
        fElementPartition[elindices[el]] = count;
        count++;
    }
    fPartitionCenterNode.resize(count);
}

/// create geometric volumetric elements
void TPZBuildSBFem::CreateVolumetricElements(TPZCompMesh &cmesh)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    cmesh.LoadReferences();
    std::set<int> matids, matidstarget;
    for (std::map<int,int>::iterator it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        long mat = it->second;
        if (cmesh.FindMaterial(mat)) {
            matids.insert(it->first);
            matidstarget.insert(it->second);
        }
    }
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->HasSubElement() || gel->Reference()) {
            continue;
        }
        if (gel->Dimension() != dim) {
            continue;
        }
        if (fElementPartition[el] == -1) {
            continue;
        }
        if (matids.find(gel->MaterialId()) == matids.end()) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gel,is);
            int onlyinterpolated = true;
            int removeduplicates = true;
            gelside.EqualorHigherCompElementList2(celstack, onlyinterpolated, removeduplicates);
            int ncelstack = celstack.NElements();
            for (int icel=0; icel<ncelstack; icel++) {
                TPZGeoElSide subgelside = celstack[icel].Reference();
                // we are only interested in faces
                if (subgelside.Dimension() != dim-1) {
                    continue;
                }
                int nnodes = subgelside.NSideNodes();
                if (nnodes != 2) {
                    std::cout << "Please extend the code to higher dimensions\n";
                    DebugStop();
                }
                TPZManVector<long,4> Nodes(nnodes+2,-1);
                int matid = fMatIdTranslation[gel->MaterialId()];
                long index;
                for (int in=0; in<nnodes; in++) {
                    Nodes[in] = subgelside.SideNodeIndex(in);
                }
                int elpartition = fElementPartition[el];
                // totototototo
//                if(elpartition != 1) continue;
                Nodes[nnodes] = fPartitionCenterNode[elpartition];
                Nodes[nnodes+1] = fPartitionCenterNode[elpartition];
                if (subgelside.IsLinearMapping())
                {
                    gmesh->CreateGeoElement(EQuadrilateral, Nodes, matid, index);
                }
                else
                {
                    long elementid = gmesh->NElements()+1;
                    TPZGeoEl *gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (Nodes, matid, *gmesh,index);

                }
                if (index >= fElementPartition.size()) {
                    fElementPartition.resize(index+1);
                }
                fElementPartition[index] = elpartition;
            }
        }
    }
    gmesh->BuildConnectivity();
    cmesh.ApproxSpace().SetAllCreateFunctionsSBFem(dim);
    cmesh.AutoBuild(matidstarget);
}

/// put the sbfem volumetric elements in element groups
void TPZBuildSBFem::CreateElementGroups(TPZCompMesh &cmesh)
{
    long numgroups = fPartitionCenterNode.size();
    long groupelementindices(numgroups);
    
    TPZVec<long> elementgroupindices(numgroups);
    
    for (long el=0; el<numgroups; el++) {
        long index;
        new TPZSBFemElementGroup(cmesh,index);
        elementgroupindices[el] = index;
    }

    
    long nel = cmesh.NElements();
    int dim = cmesh.Dimension();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if (!cel) {
            continue;
        }
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
        if (sbfem) {
            TPZGeoEl *gel = sbfem->Reference();
            long gelindex = gel->Index();
            TPZGeoElSide gelside(gel,4);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside && (neighbour.Element()->Dimension() != dim-1 || !neighbour.Element()->Reference())) {
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                // we are not handling open sides (yet)
                DebugStop();
            }
            long skelindex = neighbour.Element()->Reference()->Index();
            sbfem->SetSkeleton(skelindex);
            
            long gelgroup = fElementPartition[gelindex];
            if (gelgroup == -1) {
                DebugStop();
            }
            long celgroupindex = elementgroupindices[gelgroup];
            TPZCompEl *celgr = cmesh.Element(celgroupindex);
            TPZSBFemElementGroup *sbfemgr = dynamic_cast<TPZSBFemElementGroup *>(celgr);
            if (!sbfemgr) {
                DebugStop();
            }
            sbfemgr->AddElement(sbfem);
//            sbfem->SetElementGroupIndex(celgroupindex);
        }
    }

    for (long el=0; el<numgroups; el++) {
        long index;
        
        index = elementgroupindices[el];
        TPZCompEl *cel = cmesh.Element(index);
        TPZSBFemElementGroup *sbfemgroup = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!sbfemgroup) {
            DebugStop();
        }
        TPZStack<TPZCompEl *, 5> subgr = sbfemgroup->GetElGroup();
        long nsub = subgr.NElements();
        for (long is=0; is<nsub; is++) {
            TPZCompEl *cel = subgr[is];
            TPZSBFemVolume *femvol = dynamic_cast<TPZSBFemVolume *>(cel);
            if (!femvol) {
                DebugStop();
            }
            femvol->SetElementGroupIndex(index);
        }
    }

}

/// Divide de skeleton elements
void TPZBuildSBFem::DivideSkeleton(int nref)
{
    int dim = fGMesh->Dimension();
    for (int ir=0; ir<nref; ir++)
    {
        TPZAdmChunkVector<TPZGeoEl *> elvec = fGMesh->ElementVec();
        long nel = elvec.NElements();
        for (long el=0; el<nel; el++) {
            TPZGeoEl *gel = elvec[el];
            if (!gel || gel->HasSubElement()) {
                continue;
            }
            if (gel->Dimension() != dim-1) {
                continue;
            }
            TPZManVector<TPZGeoEl *,10> subel;
            gel->Divide(subel);
        }
    }
}
