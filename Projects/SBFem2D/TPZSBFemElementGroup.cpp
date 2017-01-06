//
//  TPZSBFemElementGroup.cpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#include "TPZSBFemElementGroup.hpp"
#include "TPZSBFemVolume.h"
#include <Accelerate/Accelerate.h>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemelementgroup"));
#endif


/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZSBFemElementGroup::ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0)
{
    std::map<long,long> locindex;
    long ncon = fConnectIndexes.size();
    for (long ic=0; ic<ncon ; ic++) {
        locindex[fConnectIndexes[ic]] = ic;
    }
    TPZElementMatrix ef(Mesh(),TPZElementMatrix::EF);
    long nel = fElGroup.size();
    InitializeElementMatrix(E0, ef);
    InitializeElementMatrix(E1, ef);
    InitializeElementMatrix(E2, ef);
    InitializeElementMatrix(M0, ef);
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        TPZElementMatrix E0Loc(Mesh(),TPZElementMatrix::EK);
        TPZElementMatrix E1Loc(Mesh(),TPZElementMatrix::EK);
        TPZElementMatrix E2Loc(Mesh(),TPZElementMatrix::EK);
        TPZElementMatrix M0Loc(Mesh(),TPZElementMatrix::EK);
        sbfem->ComputeKMatrices(E0Loc, E1Loc, E2Loc,M0Loc);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            TPZGeoEl *gel = cel->Reference();
            
            int matid = 0;
            if(gel) matid = gel->MaterialId();
            std::stringstream sout;
            if (gel) {
                sout << "Material id " << matid <<std::endl;
            }
            else
            {
                sout << "No associated geometry\n";
            }
            sout << "Connect indexes ";
            for (int i=0; i<cel->NConnects(); i++) {
                sout << cel->ConnectIndex(i) << " ";
            }
            sout << std::endl;
            sout << "Local indexes ";
            for (int i=0; i<cel->NConnects(); i++) {
                sout << locindex[cel->ConnectIndex(i)] << " ";
            }
            sout << std::endl;
            E0.fMat.Print("Matriz elementar E0",sout);
            E1.fMat.Print("Matriz elementar E1",sout);
            E2.fMat.Print("Matriz elementar E2",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
        
#endif
        int nelcon = E0Loc.NConnects();
        for (int ic=0; ic<nelcon; ic++) {
            int iblsize = E0Loc.fBlock.Size(ic);
            int icindex = E0Loc.fConnect[ic];
            int ibldest = locindex[icindex];
            for (int jc = 0; jc<nelcon; jc++) {
                int jblsize = E0Loc.fBlock.Size(jc);
                int jcindex = E0Loc.fConnect[jc];
                int jbldest = locindex[jcindex];
                for (int idf = 0; idf<iblsize; idf++) {
                    for (int jdf=0; jdf<jblsize; jdf++) {
                        E0.fBlock(ibldest,jbldest,idf,jdf) += E0Loc.fBlock(ic,jc,idf,jdf);
                        E1.fBlock(ibldest,jbldest,idf,jdf) += E1Loc.fBlock(ic,jc,idf,jdf);
                        E2.fBlock(ibldest,jbldest,idf,jdf) += E2Loc.fBlock(ic,jc,idf,jdf);
                        M0.fBlock(ibldest,jbldest,idf,jdf) += M0Loc.fBlock(ic,jc,idf,jdf);
                    }
                }
            }
        }
    }
}

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZSBFemElementGroup::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    TPZElementMatrix E0,E1,E2, M0;
    ComputeMatrices(E0, E1, E2, M0);
    
//    E0.fMat.Print("E0");
//    E1.fMat.Print("E1Check = ",std::cout,EMathematicaInput);
//    E2.fMat.Print("E2");
    
    int n = E0.fMat.Rows();
    
    TPZFMatrix<STATE> E0Inv(E0.fMat);
    TPZVec<int> pivot(E0Inv.Rows(),0);
    int nwork = 4*n*n + 2*n;
    TPZVec<STATE> work(nwork,0.);
    int info=0;
    dgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
    if (info != 0) {
        DebugStop();
    }
    dgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
    if (info != 0) {
        DebugStop();
    }
    E0Inv.Print("E0InvCheck = ",std::cout,EMathematicaInput);
    
    TPZFMatrix<STATE> globmat(2*n,2*n,0.);
    
//    cblas_dgemm(<#const enum CBLAS_ORDER __Order#>, <#const enum CBLAS_TRANSPOSE __TransA#>, <#const enum CBLAS_TRANSPOSE __TransB#>, <#const int __M#>, <#const int __N#>, <#const int __K#>, <#const double __alpha#>, <#const double *__A#>, <#const int __lda#>, <#const double *__B#>, <#const int __ldb#>, <#const double __beta#>, <#double *__C#>, <#const int __ldc#>)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1.fMat(0,0), n, 0., &globmat(0,0), 2*n);

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i,j+n) = -E0Inv(i,j);
        }
    }
    
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1.fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i+n,j) -= E2.fMat(i,j);
        }
    }

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1.fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
    globmat.Print("GlobMatCheck = ",std::cout, EMathematicaInput);

    TPZFMatrix<STATE> globmatkeep(globmat);
    TPZFMatrix< std::complex<double> > eigenVectors;
    TPZManVector<std::complex<double> > eigenvalues;
    
    globmatkeep.SolveEigenProblem(eigenvalues, eigenVectors);
    
    {
        TPZManVector<STATE> eigvalreal(2*n,0.);
        TPZFMatrix<STATE> eigvecreal(2*n,2*n,0.);
        for (int i=0; i<2*n; i++) {
            eigvalreal[i] = eigenvalues[i].real();
            for (int j=0; j<2*n; j++) {
                eigvecreal(i,j) = eigenVectors(i,j).real();
            }
        }
        std::cout << "eigval = {" << eigvalreal << "};\n";
        eigvecreal.Print("eigvec =",std::cout,EMathematicaInput);
    }
    
    TPZFNMatrix<200,STATE> QVectors(n,n,0.);
    fPhi.Resize(n, n);
    TPZManVector<STATE> eigvalsel(n,0);
    int count = 0;
    for (int i=0; i<2*n; i++) {
        if (eigenvalues[i].real() < -1.e-6) {
            for (int j=0; j<n; j++) {
                QVectors(j,count) = eigenVectors(j+n,i).real();
                fPhi(j,count) = eigenVectors(j,i).real();
                eigvalsel[count] = eigenvalues[i].real();
            }
            count++;
        }
    }
    int nstate = Connect(0).NState();
    if (nstate != 2 && nstate != 1) {
        DebugStop();
    }
    if (count != n-nstate) {
        DebugStop();
    }
    int ncon = fConnectIndexes.size();
    int eq=0;
    std::set<long> cornercon;
    BuildCornerConnectList(cornercon);
    for (int ic=0; ic<ncon; ic++) {
        long conindex = ConnectIndex(ic);
        if (cornercon.find(conindex) != cornercon.end())
        {
            fPhi(eq,count) = 1;
            if (nstate == 2)
            {
                fPhi(eq+1,count+1) = 1;
            }
        }
        eq += Connect(ic).NShape()*Connect(ic).NState();
    }
    fEigenvalues = eigvalsel;
//    std::cout << "eigenvalues " << eigvalsel << std::endl;
//    PhiVectors.Print("Phivec ",std::cout);
    TPZFMatrix<STATE> phicopy(fPhi);
    phicopy.Inverse(fPhiInverse, ELU);
//    fPhiInverse.Print("fPhiInverse",std::cout);
//    QVectors.Print("QVectors ", std::cout);
    QVectors.Multiply(fPhiInverse, ek.fMat);
    
    long nel = fElGroup.size();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->SetPhiEigVal(fPhi, fEigenvalues);
    }
    
    ComputeMassMatrix(M0);
    
    ek.fMat.Print("Stiffness",std::cout);
}

void TPZSBFemElementGroup::LoadSolution()
{
    TPZFNMatrix<200,STATE> uh_local(fPhi.Rows(),fMesh->Solution().Cols(),0.);
    int nc = NConnects();
    int count = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = Connect(ic);
        int nshape = c.NShape();
        int nstate = c.NState();
        int blsize = nshape*nstate;
        long seqnum = c.SequenceNumber();
        long pos = fMesh->Block().Position(seqnum);
        for (int seq=0; seq < blsize; seq++) {
            for (int c=0; c<uh_local.Cols(); c++)
            {
                uh_local(count+seq,c) = fMesh->Solution()(pos+seq,c);
            }
        }
        count += blsize;
    }
    fPhiInverse.Multiply(uh_local, fCoef);
    long nel = fElGroup.size();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
#ifdef PZDEBUG
        if (!sbfem) {
            DebugStop();
        }
#endif
        sbfem->LoadCoef(fCoef);
    }

}

/// Compute the mass matrix based on the value of M0 and the eigenvectors
void TPZSBFemElementGroup::ComputeMassMatrix(TPZElementMatrix &M0)
{
    //    M0.fMat.Print("Mass = ",std::cout,EMathematicaInput);
    TPZFMatrix<STATE> temp;
    REAL alpha = 1.;
    REAL beta = 0.;
    int transpose = 1;
    fPhi.MultAdd(M0.fMat, M0.fMat, temp, alpha, beta, transpose);
    //    temp.Print("Temp = ",std::cout,EMathematicaInput);
    temp.Multiply(fPhi,fMassMatrix);
    int nrow = fEigenvalues.size();
    for (int i=0; i<nrow; i++) {
        for (int j=0; j<nrow; j++) {
            fMassMatrix(i,j) /= (2.-fEigenvalues[i]-fEigenvalues[j]);
        }
    }
    fPhiInverse.MultAdd(fMassMatrix, M0.fMat, temp,alpha, beta, transpose);
    temp.Multiply(fPhiInverse,fMassMatrix);
    
    fMassMatrix.Print("Mass Matrix",std::cout);
}


//http://www.netlib.org/lapack/lug/node50.html
//https://software.intel.com/en-us/node/521079
