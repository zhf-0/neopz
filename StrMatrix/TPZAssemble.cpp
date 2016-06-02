//
//  TPZAssemble.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/31/16.
//
//

#include "TPZAssemble.h"

#include "run_stats_table.h"
#include "TPZTimer.h"

#include "pzlog.h"
#include "pzsubcmesh.h"
#include "pzmaterial.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.TPZAsemble"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
static LoggerPtr loggerel2(Logger::getLogger("pz.strmatrix.elementinterface"));
static LoggerPtr loggerelmat(Logger::getLogger("pz.strmatrix.elementmat"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.strmatrix.checkconsistency"));
static LoggerPtr loggerGlobStiff(Logger::getLogger("pz.strmatrix.globalstiffness"));
#endif


RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZAssemble::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs)
{
    ass_stiff.start();
    if (fConfig->fEquationFilter.IsActive()) {
        long neqcondense = fConfig->fEquationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense,rhs.Cols(),0.);
        if(this->fConfig->fNumThreads){
            this->MultiThread_Assemble(stiffness,rhsloc);
        }
        else{
            this->Serial_Assemble(stiffness,rhsloc);
        }
        
        fConfig->fEquationFilter.Scatter(rhsloc, rhs);
    }
    else
    {
        if(this->fConfig->fNumThreads){
            this->MultiThread_Assemble(stiffness,rhs);
        }
        else{
            this->Serial_Assemble(stiffness,rhs);
        }
    }
    ass_stiff.stop();
}

void TPZAssemble::Assemble(TPZFMatrix<STATE> & rhs){
    ass_rhs.start();
    if(fConfig->fEquationFilter.IsActive())
    {
        long neqcondense = fConfig->fEquationFilter.NActiveEquations();
        long neqexpand = fConfig->fEquationFilter.NEqExpand();
        if(rhs.Rows() != neqexpand || Norm(rhs) != 0.)
        {
            DebugStop();
        }
        TPZFMatrix<STATE> rhsloc(neqcondense,1,0.);
        if(this->fConfig->fNumThreads)
        {
            this->MultiThread_Assemble(rhsloc);
        }
        else
        {
            this->Serial_Assemble(rhsloc);
        }
        fConfig->fEquationFilter.Scatter(rhsloc,rhs);
    }
    else
    {
        if(this->fConfig->fNumThreads){
            this->MultiThread_Assemble(rhs);
        }
        else{
            this->Serial_Assemble(rhs);
        }
    }
    ass_rhs.stop();
}



void TPZAssemble::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs ){
    
    if(!fConfig->fMesh){
        LOGPZ_ERROR(logger,"Serial_Assemble called without mesh")
        DebugStop();
    }
#ifdef LOG4CXX
    if (loggerelmat->isDebugEnabled())
    {
        if(dynamic_cast<TPZSubCompMesh * >(fConfig->fMesh))
        {
            std::stringstream sout;
            sout << "AllEig = {};";
            LOGPZ_DEBUG(loggerelmat,sout.str())
        }
    }
#endif
    // tototototo
    //    TPZFNMatrix<51*51,STATE> GK(51,51,0.),GF(51,1,0.),Sol(51,1,0.),GKSol(51,1,0.);
    //    for (int i=0; i<4; i++) {
    //        Sol(47+i,0) = 1.;
    //        Sol(39+i,0) = 1.;
    //    }
    // tototototo
#ifdef PZDEBUG
    if (rhs.Rows() != fConfig->fEquationFilter.NActiveEquations()) {
        DebugStop();
    }
#endif
    
    long iel;
    long nelem = fConfig->fMesh->NElements();
    TPZElementMatrix ek(fConfig->fMesh, TPZElementMatrix::EK),ef(fConfig->fMesh, TPZElementMatrix::EF);
#ifdef LOG4CXX
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fConfig->fMesh->ElementVec();
    
    long count = 0;
    for(iel=0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        int matid = 0;
        TPZGeoEl *gel = el->Reference();
        if (gel) {
            matid = gel->MaterialId();
        }
        int matidsize = fConfig->fMaterialIds.size();
        if(matidsize){
            TPZMaterial * mat = el->Material();
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
            if (!mat)
            {
                if (!submesh) {
                    continue;
                }
                else if(submesh->NeedsComputing(fConfig->fMaterialIds) == false) continue;
            }
            else
            {
                if (this->ShouldCompute(matid) == false) continue;
            }
        }
        
        count++;
        if(!(count%1000))
        {
            std::cout << '*';
            std::cout.flush();
        }
        if(!(count%20000))
        {
            std::cout << "\n";
        }
        calcstiff.start();
        ek.Reset();
        ef.Reset();
        el->CalcStiff(ek,ef);
        
        
#ifdef LOG4CXX
        if (loggerelmat->isDebugEnabled())
        {
            if(dynamic_cast<TPZSubCompMesh * >(fConfig->fMesh))
            {
                std::stringstream objname;
                objname << "Element" << iel;
                std::string name = objname.str();
                objname << " = ";
                std::stringstream sout;
                ek.fMat.Print(objname.str().c_str(),sout,EMathematicaInput);
                sout << "AppendTo[AllEig,Eigenvalues[" << name << "]];";
                
                LOGPZ_DEBUG(loggerelmat,sout.str())
                /*		  if(iel == 133)
                 {
                 std::stringstream sout2;
                 el->Reference()->Print(sout2);
                 el->Print(sout2);
                 LOGPZ_DEBUG(logger,sout2.str())
                 }
                 */
            }
        }
#endif
        
#ifdef CHECKCONSISTENCY
        //extern TPZCheckConsistency stiffconsist("ElementStiff");
        stiffconsist.SetOverWrite(true);
        bool result;
        result = stiffconsist.CheckObject(ek.fMat);
        if(!result)
        {
            globalresult = false;
            std::stringstream sout;
            sout << "element " << iel << " computed differently";
            LOGPZ_ERROR(loggerCheck,sout.str())
        }
        
#endif
        
        calcstiff.stop();
        assemble.start();
        
        if(!ek.HasDependency()) {
            ek.ComputeDestinationIndices();
            fConfig->fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            //			TPZSFMatrix<STATE> test(stiffness);
            //			TPZFMatrix<STATE> test2(stiffness.Rows(),stiffness.Cols(),0.);
            //			stiffness.Print("before assembly",std::cout,EMathematicaInput);
            stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
            //			rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
            //			test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			test -= stiffness;
            //			test.Print("matriz de rigidez diference",std::cout);
            //			test2.Print("matriz de rigidez interface",std::cout);
            // tototototo
            //            GK.Zero();
            //            GF.Zero();
            //            GK.AddKel(ek.fMat, ek.fSourceIndex,ek.fDestinationIndex);
            //            GF.AddFel(ef.fMat, ek.fSourceIndex,ek.fDestinationIndex);
            
#ifdef LOG4CXX
            if(loggerel->isDebugEnabled())
            {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                if(gel)
                {
                    TPZManVector<REAL> center(gel->Dimension()),xcenter(3,0.);
                    gel->CenterPoint(gel->NSides()-1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for computational element index " << el->Index() << " material id " << gel->MaterialId() << std::endl;
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                }
                else {
                    sout << "Stiffness for computational element without associated geometric element\n";
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            fConfig->fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            // tototototototo
            //            GK.Zero();
            //            GF.Zero();
            //            GK.AddKel(ek.fConstrMat, ek.fSourceIndex,ek.fDestinationIndex);
            //            GF.AddFel(ef.fConstrMat, ek.fSourceIndex,ek.fDestinationIndex);
            
#ifdef LOG4CXX
            if(loggerel->isDebugEnabled() && ! dynamic_cast<TPZSubCompMesh *>(fConfig->fMesh))
            {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                //                el->Print();
                //                int nc = el->NConnects();
                //                for (int ic=0; ic<nc; ic++) {
                //                    std::cout << "Index " << el->ConnectIndex(ic) << " ";
                //                    el->Connect(ic).Print(*fMesh);
                //                    fMesh->ConnectVec()[ic].Print(*fMesh);
                //                }
                if (gel)
                {
                    TPZManVector<REAL> center(gel->Dimension()),xcenter(3,0.);
                    gel->CenterPoint(gel->NSides()-1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                }
                else{
                    sout << "Stiffness for computational element index " << iel << std::endl;
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        }
        // tototototo
        //        GK.Multiply(Sol, GKSol);
        //        GKSol -= GF;
        //        GKSol.Transpose();
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "GKSol" << iel << " = ";
        //            GKSol.Print(str.str().c_str(),sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        stiffness.Multiply(Sol, GKSol);
        //        GKSol -= rhs;
        //        GKSol.Transpose();
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "StiffSol" << iel << " = ";
        //            GKSol.Print(str.str().c_str(),sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "GK" << iel << " = ";
        //            GK.Print(str.str().c_str(),sout,EMathematicaInput);
        //            std::stringstream str2;
        //            str2 << "ST" << iel << " = ";
        //            stiffness.Print(str2.str().c_str(),sout,EMathematicaInput);
        //            sout << "GK-ST\n";
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //
        //        stiffness.Zero();
        //        rhs.Zero();
        assemble.stop();
    }//fim for iel
    if(count > 1000) std::cout << std::endl;
    
#ifdef LOG4CXX
    if(loggerCheck->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
    if (loggerel->isDebugEnabled())
    {
        std::stringstream sout;
        stiffness.Print("GK = ",sout,EMathematicaInput);
        rhs.Print("GR = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerel,sout.str())
    }
    
#endif
    
}

void TPZAssemble::Serial_Assemble(TPZFMatrix<STATE> & rhs){
    
    long iel;
    long nelem = fConfig->fMesh->NElements();
    
    TPZTimer calcresidual("Computing the residual vector");
    TPZTimer assemble("Assembling the residual vector");
    
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fConfig->fMesh->ElementVec();
    
    for(iel=0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        
        TPZMaterial * mat = el->Material();
        if (!mat) continue;
        int matid = mat->Id();
        if (this->ShouldCompute(matid) == false) continue;
        
        TPZElementMatrix ef(fConfig->fMesh, TPZElementMatrix::EF);
        
        calcresidual.start();
        
        el->CalcResidual(ef);
        
        calcresidual.stop();
        
        assemble.start();
        
        if(!ef.HasDependency()) {
            ef.ComputeDestinationIndices();
            fConfig->fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ef.ComputeDestinationIndices();
            fConfig->fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
        }
        
        assemble.stop();
        
    }//fim for iel
#ifdef LOG4CXX
    {
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << calcresidual.processName() << " " << calcresidual << std::endl;
            sout << assemble.processName() << " " << assemble;
            LOGPZ_DEBUG(logger,sout.str().c_str());
        }
    }
#endif
    //std::cout << std::endl;
}

/// filter out the equations which are out of the range
void TPZAssemble::FilterEquations(TPZVec<long> &origindex, TPZVec<long> &destindex) const
{
    //destindex = origindex;
    fConfig->fEquationFilter.Filter(origindex, destindex);
    
}

static bool CanAssemble(TPZStack<long> &connectlist, TPZVec<long> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        if (elContribute[connectlist[i]] >= 0){
            return false;
        }
    }
    return true;
}

static void AssembleColor(int el,TPZStack<long> &connectlist, TPZVec<long> &elContribute)
{
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        elContribute[connectlist[i]] = el;
    }
}

static int WhoBlockedMe(TPZStack<long> &connectlist, TPZVec<long> &elContribute, TPZVec<long> &elSeqinv)
{
    int el = -1;
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        int elBlocked = elContribute[connectlist[i]];
        if (elBlocked == -1) continue;
        int elBlockedIndex = elSeqinv[elBlocked];
        if (el == -1) el = elBlockedIndex;
        if (elBlockedIndex < el) el = elBlockedIndex;
    }
    return el;
}

static void RemoveEl(int el,TPZCompMesh *cmesh,TPZVec<long> &elContribute,long elSequence)
{
    TPZCompEl *cel = cmesh->ElementVec()[el];
    if(!cel) DebugStop();
    TPZStack<long> connectlist;
    cel->BuildConnectList(connectlist);
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        int conindex = connectlist[i];
        if (elContribute[conindex] != elSequence){
            DebugStop();
        }
        elContribute[conindex] = -1;
    }
}

static int MinPassIndex(TPZStack<long> &connectlist,TPZVec<long> &elContribute, TPZVec<int> &passIndex)
{
    int minPassIndex = -1;
    for (int i = 0 ; i < connectlist.NElements() ; i++)
    {
        int elcont = elContribute[connectlist[i]];
        int passindex = -1;
        if (elcont != -1){
            passindex = passIndex[elcont];
            if (minPassIndex == -1) minPassIndex = passindex;
        }
        if (minPassIndex < passindex) minPassIndex = passindex;
    }
    return minPassIndex;
}



void TPZAssemble::ElementColoring(TPZCompMesh *cmesh, TPZVec<long> &elSequence, TPZVec<long> &elSequenceColor,
                                        TPZVec<long> &elBlocked)
{
    
    const int nnodes = cmesh->NConnects();
    const int nel = cmesh->ElementVec().NElements();
    
    TPZManVector<long> elContribute(nnodes,-1), elSequenceColorInv(nel,-1);
    TPZManVector<int> passIndex(nel,-1);
    elSequenceColor.Resize(nel);
    elSequenceColor.Fill(-1);
    elBlocked.Resize(nel);
    elBlocked.Fill(-1);
    int nelProcessed = 0;
    int currentEl = 0;
    int currentPassIndex = 0;
    while (nelProcessed < elSequence.NElements()){
        
        int elindex = elSequence[currentEl];
        
        if(elSequenceColorInv[elindex] == -1)
        {
            TPZCompEl *cel = cmesh->ElementVec()[elindex];
            
            
            if(!cel) continue;
            TPZStack<long> connectlist;
            cel->BuildConnectList(connectlist);
            //      std::cout << "elcontribute " << elContribute << std::endl;
            //      std::cout << "connectlist " << connectlist << std::endl;
            int minPass = MinPassIndex(connectlist,elContribute,passIndex);
            if (minPass == -1){
                passIndex[elindex] = currentPassIndex;
                AssembleColor(elindex,connectlist,elContribute);
                elSequenceColor[nelProcessed] = elindex;
                elSequenceColorInv[elindex] = nelProcessed;
                nelProcessed++;
            }
            else if (minPass == currentPassIndex){
            }
            else if (minPass < currentPassIndex){
                while (!CanAssemble(connectlist,elContribute)){
                    const int el = WhoBlockedMe(connectlist,elContribute, elSequenceColorInv);
                    if (elBlocked[el] == -1) elBlocked[el] = nelProcessed;
                    int locindex = elSequenceColor[el];
                    RemoveEl(locindex,cmesh,elContribute,locindex);
                    //          std::cout << "elcontribute " << elContribute << std::endl;
                }
                passIndex[elindex] = currentPassIndex;
                AssembleColor(elindex,connectlist,elContribute);
                elSequenceColor[nelProcessed] = elindex;
                elSequenceColorInv[elindex] = nelProcessed;
                nelProcessed++;
            }
            else{
                DebugStop();
            }
        }
        currentEl++;
        if (currentEl == elSequence.NElements()){
            currentEl = 0;
            currentPassIndex++;
        }
    }
    
    //std::cout << "sequence: " << elSequence << std::endl;
    //std::cout << "color: " << elSequenceColorInv << std::endl;
    
    
    //    exit(101);
    /*
     std::ofstream toto("c:\\Temp\\output\\ColorMeshDebug.txt");
     toto << "elSequence\n" << elSequence << std::endl;
     toto << "elSequenceColor\n" << elSequenceColor << std::endl;
     toto << "elSequenceColorInv\n" << elSequenceColorInv << std::endl;
     toto << "elBlocked\n" << elBlocked << std::endl;
     toto << "elContribute\n" << elContribute << std::endl;
     toto << "passIndex\n" << passIndex << std::endl;
     toto.close();
     */
}

void TPZAssemble::OrderElement(TPZCompMesh *cmesh, TPZVec<long> &ElementOrder)
{
    
    int numelconnected = 0;
    int nconnect = cmesh->ConnectVec().NElements();
    int ic;
    //firstelconnect contains the first element index in the elconnect vector
    TPZVec<int> firstelconnect(nconnect+1);
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        numelconnected += cmesh->ConnectVec()[ic].NElConnected();
        firstelconnect[ic+1] = firstelconnect[ic]+cmesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "numelconnected " << numelconnected << endl;
    //cout << "firstelconnect ";
    //  for(ic=0; ic<nconnect; ic++) cout << firstelconnect[ic] << ' ';
    TPZVec<int> elconnect(numelconnected,-1);
    int el;
    TPZCompEl *cel;
    for(el=0; el<cmesh->ElementVec().NElements(); el++) {
        cel = cmesh->ElementVec()[el];
        if(!cel) continue;
        TPZStack<long> connectlist;
        cel->BuildConnectList(connectlist);
        int nc = connectlist.NElements();
        int ic;
        for(ic=0; ic<nc; ic++) {
            int cindex = connectlist[ic];
            elconnect[firstelconnect[cindex]] = el;
            firstelconnect[cindex]++;
        }
    }
    //  for(ic=0; ic<numelconnected; ic++) cout << elconnect[ic] << endl;
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        firstelconnect[ic+1] = firstelconnect[ic]+cmesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "elconnect\n";
    //  int no;
    //  for(no=0; no< fMesh->ConnectVec().NElements(); no++) {
    //cout << "no numero " << no << ' ' << " seq num " << fMesh->ConnectVec()[no].SequenceNumber() << ' ';
    //       for(ic=firstelconnect[no]; ic<firstelconnect[no+1];ic++) cout << elconnect[ic] << ' ';
    //cout << endl;
    //  }
    
    ElementOrder.Resize(cmesh->ElementVec().NElements(),-1);
    ElementOrder.Fill(-1);
    TPZVec<int> nodeorder(cmesh->ConnectVec().NElements(),-1);
    firstelconnect[0] = 0;
    for(ic=0; ic<nconnect; ic++) {
        int seqnum = cmesh->ConnectVec()[ic].SequenceNumber();
        if(seqnum >= 0) nodeorder[seqnum] = ic;
    }
    //  cout << "nodeorder ";
    /*  for(ic=0; ic<fMesh->ConnectVec().NElements(); ic++) cout << nodeorder[ic] << ' ';
     cout << endl;
     cout.flush();*/
    int seq;
    int elsequence = 0;
    TPZVec<int> elorderinv(cmesh->ElementVec().NElements(),-1);
    for(seq=0; seq<nconnect; seq++) {
        ic = nodeorder[seq];
        if(ic == -1) continue;
        int firstind = firstelconnect[ic];
        int lastind = firstelconnect[ic+1];
        int ind;
        for(ind=firstind; ind<lastind; ind++) {
            el = elconnect[ind];
            if(el == -1) {
                continue;
            }
            if(elorderinv[el]==-1) elorderinv[el] = elsequence++;
        }
    }
    //  cout << "elorderinv ";
    //  for(seq=0;seq<fMesh->ElementVec().NElements();seq++) cout << elorderinv[seq] << ' ';
    //  cout << endl;
    elsequence = 0;
    for(seq=0;seq<cmesh->ElementVec().NElements();seq++) {
        if(elorderinv[seq] == -1) continue;
        ElementOrder[elorderinv[seq]] = seq;
    }
    
    for(seq=0;seq<cmesh->ElementVec().NElements();seq++) {
        if(ElementOrder[seq]==-1) break;
    }
    
    ElementOrder.Resize(seq);
}