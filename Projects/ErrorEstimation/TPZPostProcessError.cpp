//
//  TPZPostProcessError.cpp
//  PZ
//
//  Created by Philippe Devloo on 6/30/16.
//
//

#include "TPZPostProcessError.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzbuildmultiphysicsmesh.h"


TPZPostProcessError::TPZPostProcessError(TPZVec<TPZCompMesh *> &meshvec)
{
    fMeshVector = meshvec;
    TPZCompMesh *multiphysics = meshvec[1];
    this->fSolution = multiphysics->Solution();
    this->fBlock = multiphysics->Block();
    fConnectSeqNumbers.resize(multiphysics->NConnects());
    for (long i=0; i<fConnectSeqNumbers.size(); i++) {
        fConnectSeqNumbers[i] = multiphysics->ConnectVec()[i].SequenceNumber();
    }
    BuildPatchStructures();
    
#ifdef PZDEBUG
    {
        std::ofstream out("../patchinfo.txt");
        PrintPatchInformation(out);
        std::ofstream out2("../mphysics.txt");
        multiphysics->Print(out2);
    }
#endif
    
}

TPZPatch TPZPostProcessError::BuildPatch(TPZCompElSide &seed)
{
    // connected : all elements that will compose the patch
    TPZStack<TPZCompElSide> connected;
    // build the set of elements which contain the node
    TPZCompMesh *cmesh = fMeshVector[1];
    TPZGeoMesh *gmesh = fMeshVector[1]->Reference();
    int meshdim = gmesh->Dimension();
    
    // the geometric mesh needs to point to the multiphysics mesh
    if (gmesh->Reference() != fMeshVector[1]) {
        DebugStop();
    }
    TPZGeoElSide gelside = seed.Reference();
    if (gelside.Dimension() != 0) {
        DebugStop();
    }
    std::set<TPZCompEl *> patchelelements;
    std::set<long> connectset;
    std::set<long> internalconnectset, boundaryconnectset;
    gelside.ConnectedCompElementList(connected,0,0);
    connected.Push(gelside.Reference());
    while (connected.size()) {
        TPZCompElSide complocside = connected.Pop();
        patchelelements.insert(complocside.Element());
        TPZCompEl *cel = complocside.Element();
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            connectset.insert(cel->ConnectIndex(ic));
        }
        TPZStack<TPZGeoElSide> highsides;
        TPZGeoElSide geolocside = complocside.Reference();
        geolocside.Element()->AllHigherDimensionSides(geolocside.Side(),meshdim-1,highsides);
        int nsides = highsides.size();
        for (int is = 0; is<nsides; is++) {
            // this is typical for HDiv approximation spaces
            if (geolocside.Element()->SideDimension(highsides[is]) != meshdim-1) {
                continue;
            }
            TPZGeoElSide geoloclocside(geolocside.Element(),highsides[is]);
            geoloclocside.HigherLevelCompElementList2(connected,0,0);
        }
    }
    // the connects to be included are those who receive all contributions
    std::map<long,int> nelconnected;
    for (std::set<TPZCompEl *>::iterator it = patchelelements.begin(); it != patchelelements.end(); it++)
    {
        TPZStack<long> nodelist;
        (*it)->BuildConnectList(nodelist);
        int nc = nodelist.size();
        for (int ic=0; ic<nc; ic++) {
            nelconnected[nodelist[ic]]++;
        }
    }
    for (std::map<long,int>::iterator it = nelconnected.begin(); it != nelconnected.end(); it++) {
        long cindex = it->first;
        int nelc = it->second;
        if (nelc == cmesh->ConnectVec()[cindex].NElConnected()) {
            internalconnectset.insert(cindex);
        }
        else
        {
            boundaryconnectset.insert(cindex);
        }
    }
    TPZPatch result;
    result.fConnectIndices.Resize(internalconnectset.size(), -1);
    result.fElIndices.Resize(patchelelements.size(), -1);
    result.fBoundaryConnectIndices.Resize(boundaryconnectset.size(), -1);
    int count = 0;
    for (std::set<long>::iterator it = internalconnectset.begin(); it != internalconnectset.end(); it++) {
        result.fConnectIndices[count++] = *it;
    }
    count = 0;
    for (std::set<TPZCompEl *>::iterator it = patchelelements.begin(); it != patchelelements.end(); it++) {
        result.fElIndices[count++] = (*it)->Index();
    }
    count = 0;
    for (std::set<long>::iterator it = boundaryconnectset.begin(); it != boundaryconnectset.end(); it++) {
        result.fBoundaryConnectIndices[count++] = *it;
    }
    
    return result;
}

// build vector of patches of a same color
void TPZPostProcessError::BuildPatchStructures()
{
    // vector indicating which connect indices of the H1 mesh have been processed
    TPZVec<int> connectprocessed(fMeshVector[0]->NConnects(),0);
    bool connectfailed = true;
    
    // load the references of all elements of the HDiv mesh
    fMeshVector[1]->Reference()->ResetReference();
    fMeshVector[1]->LoadReferences();
    
    // connectfailed will be true if there is an element patch that could not be inserted
    while (connectfailed) {
        connectfailed = false;
        // we are creating a new color
        long numvecpatch = fVecVecPatches.size();
        fVecVecPatches.resize(numvecpatch+1);
        // fillin is a vector which contains 0 if the connect has not been touched by any patch
        TPZManVector<int> fillin(fMeshVector[1]->NConnects(),0);
        // loop over all elements of the H1 mesh
        for (long el=0; el<fMeshVector[0]->NElements(); el++) {
            // take an H1 mesh
            TPZCompEl *cel = fMeshVector[0]->Element(el);
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            TPZGeoEl *gel = cel->Reference();
            int ncorner = gel->NCornerNodes();
            // patches are associated with corner nodes
            for (int i=0; i<ncorner; i++) {
                bool vertexfailed = false;
                TPZGeoElSide gelside(gel,i);
                TPZCompElSide celside(gelside.Reference());
                TPZPatch locpatch = BuildPatch(celside);
                long locconnectindex = intel->ConnectIndex(i);
                if (connectprocessed[locconnectindex] == 1) {
                    continue;
                }
                
                for (long ic = 0; ic < locpatch.fConnectIndices.size(); ic++) {
                    if (fillin[locpatch.fConnectIndices[ic]] != 0) {
                        vertexfailed = true;
                        connectfailed = true;
                        break;
                    }
                }
                for (long ic = 0; ic < locpatch.fBoundaryConnectIndices.size(); ic++) {
                    if (fillin[locpatch.fBoundaryConnectIndices[ic]] != 0) {
                        vertexfailed = true;
                        connectfailed = true;
                        break;
                    }
                }
                if (vertexfailed == false) {
                    // all systems are go !!
                    locpatch.fPartitionConnectIndex = intel->ConnectIndex(i);
                    fVecVecPatches[numvecpatch].Push(locpatch);
                    for (long ic = 0; ic < locpatch.fConnectIndices.size(); ic++) {
                        fillin[locpatch.fConnectIndices[ic]] = 1;
                    }
                    for (long ic = 0; ic < locpatch.fBoundaryConnectIndices.size(); ic++) {
                        fillin[locpatch.fBoundaryConnectIndices[ic]] = 1;
                    }
                    connectprocessed[locconnectindex] = 1;
                }
                if (vertexfailed == true) {
                    break;
                }
            }
        }
    }
}

// print the relevant information of the patches
void TPZPostProcessError::PrintPatchInformation(std::ostream &out)
{
    out << "Number of colors = " << fVecVecPatches.size() << std::endl;
    for (long color = 0; color < fVecVecPatches.size(); color++)
    {
        long numpatches = fVecVecPatches[color].size();
        out << "color number " << color << std::endl;
        for (long p = 0; p < numpatches; p++) {
            out << "patch number " << p << std::endl;
            fVecVecPatches[color][p].Print(out);
        }
    }
}

// compute the estimated H1 seminorm errors
void TPZPostProcessError::ComputeHDivSolution()
{
    TPZCompMesh *meshmixed = fMeshVector[1];
    
    long nval = fMeshVector[4]->Solution().Rows();
    for (long i=0; i<nval; i++) {
        fMeshVector[4]->Solution()(i,0) = 1.;
    }
    
    int ModelDimension = meshmixed->Dimension();
    TPZAnalysis an(meshmixed);
#ifdef PZDEBUG
    int numthreads = 0;
#else
    int numthreads = 8;
#endif
#ifdef USING_MKL2
    TPZSymetricSpStructMatrix strmat(meshmixed);
    strmat.SetNumThreads(numthreads);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(meshmixed);
    strmat.SetNumThreads(numthreads);
    strmat.SetDecomposeType(ELDLt);
    //		TPZSkylineStructMatrix strmat3(cmesh);
    //        strmat3.SetNumThreads(8);
#endif
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;

    an.Run();
    
    /** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    
    {
        std::stringstream sout;
        sout << "../" << "Poisson" << ModelDimension << "HDiv" << ".vtk";
        an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
    }

    
    an.PostProcess(1,ModelDimension);

}

// compute the estimated H1 seminorm errors
void TPZPostProcessError::ComputeElementErrors(TPZVec<STATE> &elementerrors)
{
    
    TPZCompMesh *meshmixed = fMeshVector[1];
    {
        std::ofstream out("../CMeshError.txt");
        meshmixed->Print(out);
        std::ofstream out2("../PressureCMesh.txt");
        fMeshVector[3]->Print(out2);
        std::ofstream out3("../FluxCMesh.txt");
        fMeshVector[2]->Print(out3);
    }
    meshmixed->NEquations();
    TPZCompMesh *meshweight = fMeshVector[4];
    TPZVec<TPZCompEl *> elpointers(meshmixed->NElements());
    long nelem = meshmixed->NElements();
    for (long i=0; i< nelem; i++) {
        elpointers[i] = meshmixed->Element(i);
    }
    long ncon = meshmixed->NConnects();
    long nblocks = meshmixed->Block().NBlocks();
    TPZVec<long> origseqnum(ncon);
    for (long ic=0; ic<ncon; ic++) {
        origseqnum[ic] = meshmixed->ConnectVec()[ic].SequenceNumber();
    }
    
    int ncolors = fVecVecPatches.size();
    for (int color = 0; color < ncolors; color++)
    {
        meshmixed->Solution().Zero();
        for (long el=0; el<nelem; el++) {
            meshmixed->ElementVec()[el] = elpointers[el];
        }

        TPZVec<TPZCompEl *> activel(nelem,0);
        TPZVec<long> permute(nblocks,-1);
        
        meshweight->Solution().Zero();
        int npatch = fVecVecPatches[color].size();
        long seqcount = 0;
        long nintequations = 0;
        long totalequations = 0;
        long numintconnects = 0;
        for (long patch=0; patch<npatch; patch++) {
            numintconnects += fVecVecPatches[color][patch].fConnectIndices.size();
        }
        long extseqcount = numintconnects;
        for (long patch=0; patch<npatch; patch++) {
            {
                long partitionindex = fVecVecPatches[color][patch].fPartitionConnectIndex;
                long seqnum = meshweight->ConnectVec()[partitionindex].SequenceNumber();
                meshweight->Block()(seqnum,0,0,0) = 1.;
            }
            {
                long nel = fVecVecPatches[color][patch].fElIndices.size();
                for (long el=0; el<nel; el++) {
                    long elindex = fVecVecPatches[color][patch].fElIndices[el];
                    activel[elindex] = elpointers[elindex];
                }
            }
            {
                long ncon = fVecVecPatches[color][patch].fConnectIndices.size();
                for (long ic=0; ic<ncon; ic++) {
                    long cindex = fVecVecPatches[color][patch].fConnectIndices[ic];
                    TPZConnect &c = meshmixed->ConnectVec()[cindex];
                    long seqnum = c.SequenceNumber();
                    if (seqnum != -1)
                    {
                        permute[seqnum] = seqcount++;
                        long consize = meshmixed->ConnectVec()[cindex].NShape() * meshmixed->ConnectVec()[cindex].NState();
                        nintequations += consize;
                        totalequations += consize;
                    }
                }
            }
            {
                long ncon = fVecVecPatches[color][patch].fBoundaryConnectIndices.size();
                for (long ic=0; ic<ncon; ic++) {
                    long cindex = fVecVecPatches[color][patch].fBoundaryConnectIndices[ic];
                    TPZConnect &c = meshmixed->ConnectVec()[cindex];
                    long seqnum = c.SequenceNumber();
                    if (permute[seqnum] != -1) {
                        DebugStop();
                    }
                    if (seqnum != -1)
                    {
                        permute[seqnum] = extseqcount++;
                        long consize = meshmixed->ConnectVec()[cindex].NShape() * meshmixed->ConnectVec()[cindex].NState();
                        totalequations += consize;
                    }
                }
            }

        }
        for (long ic = 0; ic<nblocks; ic++) {
            if (permute[ic] == -1) {
                permute[ic] = extseqcount++;
            }
        }
        meshmixed->Permute(permute);
        long nactiveel = 0;
        for (long el=0; el<nelem; el++) {
            meshmixed->ElementVec()[el] = activel[el];
            if (activel[el]) {
                nactiveel++;
            }
        }
        meshmixed->ComputeNodElCon();
        long nequations = meshmixed->NEquations();
        if (nequations != totalequations) {
            DebugStop();
        }
        // Saddle permute would put the pressure equations after the boundary equations (and break the code)
//        meshmixed->SaddlePermute();
        meshmixed->ExpandSolution();

        std::cout << "Number of active elements " << nactiveel << " Number of equations " << nequations
        << "Number of internal equations " << nintequations << std::endl;
        {
            std::ofstream out("../meshmixed.txt");
            meshmixed->Print(out);
        }
//        PrintPartitionDiagnostics(color, std::cout);
        nequations = meshmixed->NEquations();
        TPZAnalysis an(meshmixed,false);
        TPZSkylineStructMatrix strmat(meshmixed);
        int numthreads = 0;
        strmat.SetNumThreads(numthreads);
//        strmat.SetEquationRange(0, nequations);
        
        strmat.SetEquationRange(0, nintequations);
        an.SetStructuralMatrix(strmat);
        //		TPZSkylineStructMatrix strmat3(cmesh);
        //        strmat3.SetNumThreads(8);
        
        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
        direct->SetDirect(ELDLt);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
        
        an.Assemble();
        
        TPZAutoPointer<TPZMatrix<STATE> > globmat = an.Solver().Matrix();
        for (long p=0; p<npatch; p++) {
            TPZPatch &patch = fVecVecPatches[color][p];
            if (!PatchHasBoundary(patch))
            {
                long firstlagrangeequation = patch.FirstLagrangeEquation(meshmixed);
                STATE diag = globmat->GetVal(firstlagrangeequation, firstlagrangeequation);
                diag += 1.;
                globmat->Put(firstlagrangeequation, firstlagrangeequation, diag);
            }
        }
        
        an.Solve();
        
        TPZStepSolver<STATE> &step = dynamic_cast<TPZStepSolver<STATE> &>(an.Solver());
        std::list<long> singular = step.Singular();
        
        if (singular.size())
        {
            std::cout << "The following equations were flagged as singular\n";
            for (std::list<long>::iterator it=singular.begin(); it != singular.end(); it++) {
                std::cout << *it << " ";
            }
            std::cout << std::endl;
        }
        else
        {
            std::cout << "The matrix has no singularity\n";
        }

        // now we have a partial solution
        /** Variable names for post processing */
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("POrder");
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        an.SetStep(color);

        {
            int ModelDimension = meshmixed->Dimension();
            std::stringstream sout;
            sout << "../" << "Poisson" << ModelDimension << "HDiv" << ".vtk";
            an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
        }
        
        an.PostProcess(1,meshmixed->Dimension());

        
        for (long el=0; el<nelem; el++) {
            meshmixed->ElementVec()[el] = elpointers[el];
        }
        meshmixed->ComputeNodElCon();
        TransferAndSumSolution(meshmixed);
        ResetState();

        an.SetStep(color);

        // now draw the full mesh with the summed flux
        {
            int ModelDimension = meshmixed->Dimension();
            std::stringstream sout;
            sout << "../" << "FullPoisson" << ModelDimension << "HDiv" << ".vtk";
            an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
        }
        
        an.PostProcess(1,meshmixed->Dimension());
        fMeshVector[2]->Solution().Zero();
        fMeshVector[3]->Solution().Zero();
        an.Solution().Zero();

    }
}

// print partition diagnostics
void TPZPostProcessError::PrintPartitionDiagnostics(int color, std::ostream &out) const
{
    TPZCompMesh *meshmixed = fMeshVector[1];
    int meshdim = meshmixed->Dimension();
    if (color < 0 || color >= fVecVecPatches.size()) {
        DebugStop();
    }
    TPZVec<TPZPatch> &vecpatch = fVecVecPatches[color];
    long numpatch = vecpatch.size();
    // determine if the patch is a boundary patch or not
    TPZVec<int> IsInternalPatch(numpatch,0);
    for (long p = 0; p<numpatch; p++) {
        IsInternalPatch[p] = !PatchHasBoundary(vecpatch[p]);
    }
    out << "Number of patches " << numpatch << std::endl;
    for (long p = 0; p<numpatch; p++) {
        TPZPatch &patch = vecpatch[p];
        out << "Diagnostics for patch number " << p << " of color " << color << std::endl;
        if(IsInternalPatch[p])
        {
            out << "Patch is internal\n";
        }
        else
        {
            out << "Patch has boundary sides\n";
        }
        
        out << "Internal connects sequence numbers ";
        for (long ci=0; ci < patch.fConnectIndices.size(); ci++) {
            TPZConnect &c = meshmixed->ConnectVec()[patch.fConnectIndices[ci]];
            if (c.HasDependency()) {
                out << "*";
            }
            if (c.LagrangeMultiplier() != 0) {
                out << "L";
            }
            out << c.SequenceNumber() << " ";
        }
        out << std::endl;
        out << "Boundary connects sequence numbers ";
        for (long ci=0; ci < patch.fBoundaryConnectIndices.size(); ci++) {
            TPZConnect &c = meshmixed->ConnectVec()[patch.fBoundaryConnectIndices[ci]];
            if (c.HasDependency()) {
                out << "*ERROR*";
            }
            out << c.SequenceNumber() << " ";
        }
        out << std::endl;
    }
}

// determine if a given patch is boundary or not
bool TPZPostProcessError::PatchHasBoundary(TPZPatch &patch) const
{
    TPZCompMesh *meshmixed = fMeshVector[1];
    int meshdim = meshmixed->Dimension();
    long numel = patch.fElIndices.size();
    bool HasBoundary = false;
    for (long el=0; el<numel; el++) {
        long elindex = patch.fElIndices[el];
        TPZCompEl *cel = meshmixed->Element(elindex);
        if (!cel) {
            DebugStop();
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        if (gel->Dimension() != meshdim) {
            HasBoundary = true;
            break;
        }
    }
    return HasBoundary;

}

// return the first equation associated with a lagrange multiplier
long TPZPatch::FirstLagrangeEquation(TPZCompMesh *cmesh) const
{
    long nconnect = fConnectIndices.size();
    for (long ic=0; ic<nconnect; ic++) {
        long cindex = fConnectIndices[ic];
        TPZConnect &c = cmesh->ConnectVec()[cindex];
        if (c.SequenceNumber() == -1 || c.NDof() == 0 || c.LagrangeMultiplier() == 0) {
            continue;
        }
        long seqnum = c.SequenceNumber();
        long eq = cmesh->Block().Position(seqnum);
        return eq;
    }
    DebugStop();
    return -1;
}

// Sum the solution stored in fSolution of the second mesh to the fSolution vector
void TPZPostProcessError::TransferAndSumSolution(TPZCompMesh *cmesh)
{
    long nconnect = cmesh->NConnects();
    for (long ic=0; ic<nconnect; ic++) {
        TPZConnect &c = cmesh->ConnectVec()[ic];
        if (c.SequenceNumber() == -1) {
            continue;
        }
        long seqnum = c.SequenceNumber();
        int neq = c.NDof();
        long pos = cmesh->Block().Position(seqnum);
        
        long targetseqnum = fConnectSeqNumbers[ic];
        long targetpos = fBlock.Position(targetseqnum);
        
        
        for (int eq=0; eq<neq; eq++) {
            // tototototo
            fSolution(targetpos+eq,0) += cmesh->Solution()(pos+eq,0);
        }
    }
    cmesh->Solution().Zero();
}

// Reset the state of the HDiv mesh to its original structure
void TPZPostProcessError::ResetState()
{
    TPZCompMesh *multiphysics = fMeshVector[1];
    multiphysics->Block() = this->fBlock;
    multiphysics->Solution() = this->fSolution;
    for (long i=0; i<fConnectSeqNumbers.size(); i++) {
        multiphysics->ConnectVec()[i].SetSequenceNumber(fConnectSeqNumbers[i]);
    }
    TPZManVector<TPZCompMesh *> meshvec(2);
    meshvec[0] = fMeshVector[2];
    meshvec[1] = fMeshVector[3];
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, multiphysics);

}
