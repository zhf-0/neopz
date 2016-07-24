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


/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZSBFemElementGroup::ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2)
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
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = fElGroup[el];
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
        TPZElementMatrix E0Loc(Mesh(),TPZElementMatrix::EK);
        TPZElementMatrix E1Loc(Mesh(),TPZElementMatrix::EK);
        TPZElementMatrix E2Loc(Mesh(),TPZElementMatrix::EK);
        sbfem->ComputeKMatrices(E0Loc, E1Loc, E2Loc);
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
    TPZElementMatrix E0,E1,E2;
    ComputeMatrices(E0, E1, E2);
    
    E0.fMat.Print("E0");
    E1.fMat.Print("E1Check = ",std::cout,EMathematicaInput);
    E2.fMat.Print("E2");
    
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

    
//    int dgehrd_(__CLPK_integer *__n, __CLPK_integer *__ilo, __CLPK_integer *__ihi,
//                __CLPK_doublereal *__a, __CLPK_integer *__lda, __CLPK_doublereal *__tau,
//                __CLPK_doublereal *__work, __CLPK_integer *__lwork,
//                __CLPK_integer *__info) __OSX_AVAILABLE_STARTING(__MAC_10_2,
//                                                                 __IPHONE_4_0);

    TPZFMatrix<STATE> globmatkeep(globmat);
    int globsize = 2*n;
    int one = 1;
    TPZVec<STATE> tau(globsize,0.);
    dgehrd_(&globsize, &one, &globsize, &globmat(0,0), &globsize, &tau[0], &work[0], &nwork, &info);
    if (info != 0) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> Qorth(globmat);
    dorghr_(&globsize, &one, &globsize, &Qorth(0,0), &globsize, &tau[0], &work[0], &nwork, &info);
    if (info != 0) {
        DebugStop();
    }
    Qorth.Print("Matriz Q");
    
    TPZFMatrix<STATE> Identity(globsize,globsize,0.);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, globsize, globsize, globsize, 1., &Qorth(0,0), globsize, &Qorth(0,0), globsize, 0., &Identity(0,0), globsize);
    
    Identity.Print("Identity");
    
    TPZFMatrix<STATE> check(globsize,globsize,0.);
    for (int i=0; i<globsize; i++) {
        for (int j=0; j< globsize; j++) {
            if (i-j <= 1) {
                check(i,j) = -globmat(i,j);
            }
            for (int k=0; k<globsize; k++) {
                for (int l=0; l<globsize; l++) {
                    check(i,j) += Qorth(k,i)*globmatkeep(k,l)*Qorth(l,j);
                }
            }
        }
    }
    std::cout << "Norm of check " << Norm(check) << std::endl;
    check.Print("check");
    
    TPZVec<STATE> realeig(globsize,0.), imageig(globsize,0.);
    char S = 'S', V = 'V';
    dhseqr_(&S, &V, &globsize, &one, &globsize, &globmat(0,0), &globsize, &realeig[0], &imageig[0], &Qorth(0,0), &globsize, &work[0], &nwork, &info);

    /// globmat now contains the T matrix
    if (info != 0) {
        DebugStop();
    }
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, globsize, globsize, globsize, 1., &Qorth(0,0), globsize, &Qorth(0,0), globsize, 0., &Identity(0,0), globsize);
    
    Identity.Print("Identity");
    
    char side = 'R', eigsrc = 'Q',initv = 'N', howmny = 'A';
    TPZVec<int> select(globsize,1);
    int mm_required = 0;
    TPZFMatrix<STATE> LeftEigVec(globsize,globsize,0.), RightEigVec(globsize,globsize,0.);
    dtrevc_(<#char *__side#>, <#char *__howmny#>, __CLPK_logical *__select, <#__CLPK_integer *__n#>, <#__CLPK_doublereal *__t#>, <#__CLPK_integer *__ldt#>, <#__CLPK_doublereal *__vl#>, <#__CLPK_integer *__ldvl#>, <#__CLPK_doublereal *__vr#>, <#__CLPK_integer *__ldvr#>, <#__CLPK_integer *__mm#>, <#__CLPK_integer *__m#>, <#__CLPK_doublereal *__work#>, <#__CLPK_integer *__info#>)
    
    TPZVec<int> failL(globsize,0), failR(globsize,0);
    dhsein_(&side, &eigsrc, &initv, &select[0], &globsize, &globmat(0,0), &globsize, &realeig[0], &imageig[0], &LeftEigVec(0,0), &globsize, &RightEigVec(0,0), &globsize, &globsize, &mm_required, &work[0], &failL[0], &failR[0], &info);

    if (info != 0) {
        DebugStop();
    }
    

    std::cout << "Eig real " << realeig << std::endl;
    std::cout << "Eig imag " << imageig << std::endl;
    
    Qorth.Print("EigVecCheck = ",std::cout,EMathematicaInput);
}

//http://www.netlib.org/lapack/lug/node50.html
//https://software.intel.com/en-us/node/521079
