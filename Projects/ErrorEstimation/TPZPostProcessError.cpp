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
    std::set<long> internalconnectset;
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
    }
    TPZPatch result;
    result.fConnectIndices.Resize(internalconnectset.size(), -1);
    result.fElIndices.Resize(patchelelements.size(), -1);
    result.fAllConnectIndices.Resize(connectset.size(), -1);
    int count = 0;
    for (std::set<long>::iterator it = internalconnectset.begin(); it != internalconnectset.end(); it++) {
        result.fConnectIndices[count++] = *it;
    }
    count = 0;
    for (std::set<TPZCompEl *>::iterator it = patchelelements.begin(); it != patchelelements.end(); it++) {
        result.fElIndices[count++] = (*it)->Index();
    }
    count = 0;
    for (std::set<long>::iterator it = connectset.begin(); it != connectset.end(); it++) {
        result.fAllConnectIndices[count++] = *it;
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
                
                for (long ic = 0; ic < locpatch.fAllConnectIndices.size(); ic++) {
                    if (fillin[locpatch.fAllConnectIndices[ic]] != 0) {
                        vertexfailed = true;
                        connectfailed = true;
                        break;
                    }
                }
                if (vertexfailed == false) {
                    // all systems are go !!
                    locpatch.fPartitionConnectIndex = intel->ConnectIndex(i);
                    fVecVecPatches[numvecpatch].Push(locpatch);
                    for (long ic = 0; ic < locpatch.fAllConnectIndices.size(); ic++) {
                        fillin[locpatch.fAllConnectIndices[ic]] = 1;
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
        TPZVec<TPZCompEl *> activel(nelem,0);
        TPZVec<long> permute(nblocks,-1);
        
        meshweight->Solution().Zero();
        int npatch = fVecVecPatches[color].size();
        long seqcount = 0;
        long nequations = 0;
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
                long ncon = fVecVecPatches[color][patch].fAllConnectIndices.size();
                for (long ic=0; ic<ncon; ic++) {
                    long cindex = fVecVecPatches[color][patch].fAllConnectIndices[ic];
                    TPZConnect &c = meshmixed->ConnectVec()[cindex];
                    long seqnum = c.SequenceNumber();
                    if (seqnum != -1)
                    {
                        permute[seqnum] = seqcount++;
                        nequations += meshmixed->ConnectVec()[cindex].NShape() * meshmixed->ConnectVec()[cindex].NState();
                    }
                }
            }
        }
        for (long ic = 0; ic<nblocks; ic++) {
            if (permute[ic] == -1) {
                permute[ic] = seqcount++;
            }
        }
        meshmixed->Permute(permute);
        meshmixed->SaddlePermute();
        meshmixed->ExpandSolution();
        long nactiveel = 0;
        for (long el=0; el<nelem; el++) {
            meshmixed->ElementVec()[el] = activel[el];
            if (activel[el]) {
                nactiveel++;
            }
        }

        std::cout << "Number of active elements " << nactiveel << " Number of equations " << nequations << std::endl;
        meshmixed->ComputeNodElCon();
        {
            std::ofstream out("../meshmixed.txt");
            meshmixed->Print(out);
        }
        long neq = meshmixed->NEquations();
        TPZAnalysis an(meshmixed,false);
        an.SetStep(color);
        TPZSkylineStructMatrix strmat(meshmixed);
        int numthreads = 0;
        strmat.SetNumThreads(numthreads);
        strmat.SetEquationRange(0, nequations);
        //		TPZSkylineStructMatrix strmat3(cmesh);
        //        strmat3.SetNumThreads(8);
        
        TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
        direct->SetDirect(ELDLt);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
        
        an.Run();
        
        // now we have a partial solution
        /** Variable names for post processing */
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("POrder");
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        
        {
            int ModelDimension = meshmixed->Dimension();
            std::stringstream sout;
            sout << "../" << "Poisson" << ModelDimension << "HDiv" << ".vtk";
            an.DefineGraphMesh(ModelDimension,scalnames,vecnames,sout.str());
        }

        an.PostProcess(1,meshmixed->Dimension());
    }
}
