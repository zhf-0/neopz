//
// Created by Francisco Teixeira Orlandini on 11/27/17.
//

#include "TPZLapackWrapper.h"

template<class TVar>
int TPZLapackWrapper<TVar>::SolveEigenProblem(TPZFMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type  > &eigenVectors){
    DebugStop();
}
template<class TVar>
int TPZLapackWrapper<TVar>::SolveEigenProblem(TPZFMatrix<TVar> &A, TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
    DebugStop();
}
template<class TVar>
int TPZLapackWrapper<TVar>::SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &A, TPZFMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w, TPZFMatrix < typename SPZAlwaysComplex<TVar>::type > &eigenVectors){
    DebugStop();
}
template<class TVar>
int TPZLapackWrapper<TVar>::SolveGeneralisedEigenProblem(TPZFMatrix<TVar> &A, TPZFMatrix< TVar > &B , TPZVec < typename SPZAlwaysComplex<TVar>::type > &w){
    DebugStop();
}
/*******************
*    TPZFMATRIX    *
*******************/

#ifdef USING_LAPACK
/** CBlas Math Library */
#ifdef MACOSX
#include <Accelerate/Accelerate.h>
typedef __CLPK_doublecomplex vardoublecomplex;
typedef __CLPK_complex varfloatcomplex;
#elif USING_MKL
#include <mkl.h>
typedef MKL_Complex16 vardoublecomplex;
typedef MKL_Complex8 varfloatcomplex;
#else
#include "cblas.h"
#define BLAS_MULT
#endif
#endif

#ifdef USING_MKL
#include <mkl.h>
typedef MKL_Complex16 vardoublecomplex;
typedef MKL_Complex8 varfloatcomplex;
#endif
template <>
int TPZLapackWrapper<float>::SolveEigenProblem(TPZFMatrix<float> &A,TPZVec < typename SPZAlwaysComplex<float>::type > &eigenvalues)
{
    if (A.Rows() != A.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< float > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<float> realeigen(dim,0.);
    TPZVec<float> imageigen(dim,0.);

    TPZVec<float> work(lwork);
    sgeev_(jobvl, jobvr, &dim, A.fElem, &dim, &realeigen[0], &imageigen[0], VL.fElem, &dim, VR.fElem, &dim, &work[0], &lwork, &info);

    if (info != 0) {
        DebugStop();
    }

    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){

        eigenvalues[i] = realeigen[i] + I*imageigen[i];
    }
    return 1;
}

template <>
int TPZLapackWrapper<double>::SolveEigenProblem(TPZFMatrix<double> &A, TPZVec < typename SPZAlwaysComplex<double>::type > &eigenvalues)
{
    if (A.Rows() != A.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< double > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);

    TPZVec<double> work(lwork);
    dgeev_(jobvl, jobvr, &dim, A.fElem, &dim, &realeigen[0], &imageigen[0], VL.fElem, &dim, VR.fElem, &dim, &work[0], &lwork, &info);

    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = realeigen[i] + I*imageigen[i];
    }
    return 1;
}

template <>
int TPZLapackWrapper<std::complex<float>>::SolveEigenProblem(TPZFMatrix<std::complex<float>> &A, TPZVec <  typename SPZAlwaysComplex<std::complex<float>>::type > &eigenvalues)
{
    if (A.Rows() != A.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< std::complex<float> > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<std::complex<float> > eigen(dim,0.);

    TPZVec<std::complex<float> > work(lwork);
    TPZVec< float > rwork( 2 * dim);

    cgeev_(jobvl, jobvr, &dim, (varfloatcomplex *)A.fElem, &dim, (varfloatcomplex *)&eigen[0], (varfloatcomplex *)VL.fElem, &dim, (varfloatcomplex *)VR.fElem, &dim, (varfloatcomplex *)&work[0], &lwork, &rwork[0], &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = eigen[i];
    }
    return 1;
}

template <>
int TPZLapackWrapper<std::complex<double>>::SolveEigenProblem(TPZFMatrix<std::complex<double>> &A, TPZVec < typename SPZAlwaysComplex<std::complex<double>>::type> &eigenvalues)
{
    if (A.Rows() != A.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< std::complex<double> > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<std::complex<double> > eigen(dim,0.);

    TPZVec<std::complex<double> > work(lwork);
    TPZVec< double > rwork( 2 * dim);

    zgeev_(jobvl, jobvr, &dim, (vardoublecomplex *)A.fElem, &dim, (vardoublecomplex *)&eigen[0], (vardoublecomplex *)VL.fElem, &dim, (vardoublecomplex *)VR.fElem, &dim, (vardoublecomplex *)&work[0], &lwork, &rwork[0], &info);

    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = eigen[i];
    }

    return 1;
}

template <>
int TPZLapackWrapper<float>::SolveEigenProblem(TPZFMatrix<float> &A, TPZVec < SPZAlwaysComplex<float>::type > &eigenvalues, TPZFMatrix < SPZAlwaysComplex<float>::type > &eigenvectors)
{
    if (A.Rows() != A.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< float > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<float> realeigen(dim,0.);
    TPZVec<float> imageigen(dim,0.);
    TPZVec<float> work(lwork);

    sgeev_(jobvl, jobvr, &dim, A.fElem, &dim, &realeigen[0], &imageigen[0], VL.fElem, &dim, VR.fElem, &dim, &work[0], &lwork, &info);

    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = realeigen[i] + I*imageigen[i];
    }
    for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }

    return 1;
}

template <>
int TPZLapackWrapper<double>::SolveEigenProblem(TPZFMatrix<double> &A, TPZVec < SPZAlwaysComplex<double>::type > &eigenvalues, TPZFMatrix < SPZAlwaysComplex<double>::type > &eigenvectors)
{
    if (A.Rows() != A.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< double > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);
    TPZVec<double> work(lwork);
    dgeev_(jobvl, jobvr, &dim, A.fElem, &dim, &realeigen[0], &imageigen[0], VL.fElem, &dim, VR.fElem, &dim, &work[0], &lwork, &info);

    if (info != 0) {
        DebugStop();
    }
    //    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = realeigen[i] + I*imageigen[i];
    }
    for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            double *realVRptr = VR.fElem;
            double *imagVRptr = VR.fElem + dim;
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }

    return 1;
}

template <>
int TPZLapackWrapper<std::complex<double> >::SolveEigenProblem(TPZFMatrix<std::complex<double>> &A, TPZVec < SPZAlwaysComplex<std::complex<double>>::type > &eigenvalues, TPZFMatrix < SPZAlwaysComplex<std::complex<double>>::type > &eigenvectors)
{
    if (A.Rows() != A.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< std::complex<double> > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<std::complex<double> > eigen(dim,0.);
    TPZVec<std::complex<double> > work(lwork);
    TPZVec< double > rwork( 2 * dim);

zgeev_(jobvl, jobvr, &dim, (vardoublecomplex *)A.fElem, &dim, (vardoublecomplex *)&eigen[0], (vardoublecomplex *)VL.fElem, &dim, (vardoublecomplex *)VR.fElem, &dim, (vardoublecomplex *)&work[0], &lwork, &rwork[0], &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = eigen[i];
    }
    for(int i = 0 ; i < dim ; i ++){

        for( int iV = 0 ; iV < dim ; iV++ ){
            eigenvectors(iV,i) = VR(iV,i);
        }

    }

    return 1;
}

template <>
int TPZLapackWrapper<std::complex< float> >::SolveEigenProblem(TPZFMatrix<std::complex<float>> &A, TPZVec < SPZAlwaysComplex<std::complex<float>>::type > &eigenvalues, TPZFMatrix < SPZAlwaysComplex<std::complex<float>>::type > &eigenvectors)
{
    if (A.Rows() != A.Cols()) {
        DebugStop();
    }
    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< std::complex<float> > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<std::complex<float> > eigen(dim,0.);
    TPZVec<std::complex<float> > work(lwork);
    TPZVec< float > rwork( 2 * dim);

    cgeev_(jobvl, jobvr, &dim, (varfloatcomplex *)A.fElem, &dim, (varfloatcomplex *)&eigen[0], (varfloatcomplex *)VL.fElem, &dim, (varfloatcomplex *)VR.fElem, &dim, (varfloatcomplex *)&work[0], &lwork, &rwork[0], &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        eigenvalues[i] = eigen[i];
    }
    for(int i = 0 ; i < dim ; i ++){

        for( int iV = 0 ; iV < dim ; iV++ ){
            eigenvectors(iV,i) = VR(iV,i);
        }

    }

    return 1;
}

template<>
int TPZLapackWrapper<float>::SolveGeneralisedEigenProblem(TPZFMatrix<float> &A, TPZFMatrix<float> &B , TPZVec <SPZAlwaysComplex<float>::type > &eigenvalues, TPZFMatrix < SPZAlwaysComplex<float>::type > &eigenvectors)
{
    if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
    {
        TPZFMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }

    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< float > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<float> realeigen(dim,0.);
    TPZVec<float> imageigen(dim,0.);

    TPZVec<float> beta(dim);
    TPZVec<float> work(lwork);

    sggev_(jobvl, jobvr, &dim, A.fElem, &dim , B.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = (realeigen[i] + I*imageigen[i]) / beta[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }
    return 1;
}


template<>
int
TPZLapackWrapper<float>::SolveGeneralisedEigenProblem(TPZFMatrix<float> &A, TPZFMatrix<float> &B , TPZVec <SPZAlwaysComplex<float>::type > &eigenvalues)
{
    if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
    {
        TPZFMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }

    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< float > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<float> I(0,1.);
    TPZVec<float> realeigen(dim,0.);
    TPZVec<float> imageigen(dim,0.);

    TPZVec<float> beta(dim);
    TPZVec<float> work(lwork);

    sggev_(jobvl, jobvr, &dim, A.fElem, &dim , B.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = (realeigen[i] + I*imageigen[i]) / beta[i];
        }
    }

    return 1;
}

template<>
int
TPZLapackWrapper<double>::SolveGeneralisedEigenProblem(TPZFMatrix<double> &A ,TPZFMatrix<double> &B , TPZVec <SPZAlwaysComplex<double>::type > &eigenvalues, TPZFMatrix < SPZAlwaysComplex<double>::type > &eigenvectors)
{
    if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
    {
        TPZFMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }

    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< double > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);

    TPZVec<double> beta(dim);

    TPZVec<double> work(lwork);

dggev_(jobvl, jobvr, &dim, A.fElem, &dim , B.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = (realeigen[i] + I*imageigen[i]) / beta[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        if(imageigen[i] == 0){
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i);
            }
        }
        else{
            for( int iV = 0 ; iV < dim ; iV++ ){
                eigenvectors(iV,i) = VR(iV,i) + I * VR(iV,i+1) ;
                eigenvectors(iV,i + 1) = VR(iV,i) - I * VR(iV,i+1) ;
            }
            i++;
        }
    }

    return 1;
}


template<>
int
TPZLapackWrapper<double>::SolveGeneralisedEigenProblem(TPZFMatrix<double> &A, TPZFMatrix<double> &B , TPZVec <SPZAlwaysComplex<double>::type > &eigenvalues)
{
    if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
    {
        TPZFMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }

    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< double > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    double testwork;
    int lwork = 10+20*dim;
    int info;
    std::complex<double> I(0,1.);
    TPZVec<double> realeigen(dim,0.);
    TPZVec<double> imageigen(dim,0.);

    TPZVec<double> beta(dim);

    TPZVec<double> work(lwork);

dggev_(jobvl, jobvr, &dim, A.fElem, &dim , B.fElem, &dim , &realeigen[0], &imageigen[0], &beta[0]  , VL.fElem, &dim , VR.fElem, &dim, &work[0], &lwork, &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = (realeigen[i] + I*imageigen[i]) / beta[i];
        }
    }

    return 1;
}

template<>
int
TPZLapackWrapper<std::complex<float> >::SolveGeneralisedEigenProblem(TPZFMatrix<std::complex<float> > &A, TPZFMatrix<std::complex<float> > &B , TPZVec <SPZAlwaysComplex<std::complex<float>>::type> &eigenvalues, TPZFMatrix < SPZAlwaysComplex<std::complex<float>>::type > &eigenvectors)
{
    if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
    {
        TPZFMatrix<std::complex<float>>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }

    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< std::complex<float> > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    TPZVec<std::complex<float> > eigen(dim,0.);

    TPZVec<std::complex<float> > beta(dim);

    TPZVec<std::complex<float> > work(lwork);
    TPZVec<float> rwork( 8 * dim );

cggev_(jobvl, jobvr, &dim, (varfloatcomplex *)A.fElem, &dim , (varfloatcomplex *)B.fElem, &dim , (varfloatcomplex *)&eigen[0], (varfloatcomplex *)&beta[0]  , (varfloatcomplex *)VL.fElem, &dim , (varfloatcomplex *)VR.fElem, &dim, (varfloatcomplex *)&work[0], &lwork, &rwork[0], &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = eigen[i] / beta[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        for( int iV = 0 ; iV < dim ; iV++ ){
            eigenvectors(iV,i) = VR(iV,i);
        }
    }

    return 1;
}


template<>
int
TPZLapackWrapper<std::complex<float> >::SolveGeneralisedEigenProblem(TPZFMatrix<std::complex<float> > &A, TPZFMatrix<std::complex<float> > &B , TPZVec <SPZAlwaysComplex<std::complex<float>>::type> &eigenvalues)
{
    if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
    {
        TPZFMatrix<std::complex<float>>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }

    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< std::complex<float> > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    TPZVec<std::complex<float> > eigen(dim,0.);

    TPZVec<std::complex<float> > beta(dim);

    TPZVec<std::complex<float> > work(lwork);
    TPZVec<float> rwork( 8 * dim );

    cggev_(jobvl, jobvr, &dim, (varfloatcomplex *)A.fElem, &dim , (varfloatcomplex *)B.fElem, &dim , (varfloatcomplex *)&eigen[0], (varfloatcomplex *)&beta[0]  , (varfloatcomplex *)VL.fElem, &dim , (varfloatcomplex *)VR.fElem, &dim, (varfloatcomplex *)&work[0], &lwork, &rwork[0], &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = eigen[i] / beta[i];
        }
    }
    return 1;

}

template<>
int
TPZLapackWrapper<std::complex<double> >::SolveGeneralisedEigenProblem(TPZFMatrix<std::complex<double> > &A, TPZFMatrix<std::complex<double> > &B , TPZVec <SPZAlwaysComplex<std::complex<double>>::type> &eigenvalues, TPZFMatrix < SPZAlwaysComplex<std::complex<double>>::type > &eigenvectors)
{
    if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
    {
        TPZFMatrix<std::complex<double>>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }

    char jobvl[] = "None", jobvr[] = "Vectors";
    TPZFMatrix< std::complex<double> > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    TPZVec<std::complex<double> > eigen(dim,0.);

    TPZVec<std::complex<double> > beta(dim);

    TPZVec<std::complex<double> > work(lwork);
    TPZVec<double> rwork( 8 * dim );

    zggev_(jobvl, jobvr, &dim, (vardoublecomplex *)A.fElem, &dim , (vardoublecomplex *)B.fElem, &dim , (vardoublecomplex *)&eigen[0], (vardoublecomplex *)&beta[0]  , (vardoublecomplex *)VL.fElem, &dim , (vardoublecomplex *)VR.fElem, &dim, (vardoublecomplex *)&work[0], &lwork, &rwork[0], &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvectors.Redim(dim,dim);
    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = eigen[i] / beta[i];
        }
    }
    for(int i = 0 ; i < dim ; i ++){
        for( int iV = 0 ; iV < dim ; iV++ ){
            eigenvectors(iV, i) = VR(iV,i);
        }
    }

    return 1;
}


template<>
int
TPZLapackWrapper<std::complex<double> >::SolveGeneralisedEigenProblem(TPZFMatrix<std::complex<double> > &A, TPZFMatrix<std::complex<double> > &B , TPZVec <SPZAlwaysComplex<std::complex<double>>::type> &eigenvalues)
{
    if (  A.Rows() != B.Rows() || A.Cols() != B.Cols() || A.Cols() != A.Cols() )
    {
        TPZFMatrix<std::complex<double>>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }

    char jobvl[] = "None", jobvr[] = "None";
    TPZFMatrix< std::complex<double> > VL(A.Rows(),A.Cols()),VR(A.Rows(),A.Cols());
    int dim = A.Rows();
    float testwork;
    int lwork = 10+20*dim;
    int info;
    TPZVec<std::complex<double> > eigen(dim,0.);

    TPZVec<std::complex<double> > beta(dim);
    TPZVec<std::complex<double> > work(lwork);
    TPZVec<double> rwork( 8 * dim );

    zggev_(jobvl, jobvr, &dim, (vardoublecomplex *)A.fElem, &dim , (vardoublecomplex *)B.fElem, &dim , (vardoublecomplex *)&eigen[0], (vardoublecomplex *)&beta[0]  , (vardoublecomplex *)VL.fElem, &dim , (vardoublecomplex *)VR.fElem, &dim, (vardoublecomplex *)&work[0], &lwork, &rwork[0], &info);

    if (info != 0) {
        DebugStop();
    }
//    VR.Print("VR = ",std::cout,EMathematicaInput);

    eigenvalues.Resize(dim,0.);
    for(int i = 0 ; i < dim ; i ++){
        if( IsZero(beta[i])){
            DebugStop(); //fran: i really dont know what to do with this result
        }
        else{
            eigenvalues[i] = eigen[i] / beta[i];
        }
    }

    return 1;

}


template class TPZLapackWrapper< std::complex<float> >;
template class TPZLapackWrapper< std::complex<double> >;
template class TPZLapackWrapper< std::complex<long double> >;

template class TPZLapackWrapper<long >;

template class TPZLapackWrapper<float >;
template class TPZLapackWrapper<double >;
template class TPZLapackWrapper<long double>;

template class TPZLapackWrapper<int >;
template class TPZLapackWrapper<TPZFlopCounter>;
