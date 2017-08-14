/**
 * @file
 */


#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"

#include "TPZSandlerExtendedPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticResponse.h"

#include "pzlog.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.material.TPZPlasticStepPV"));
//#endif
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic2"));
#endif

#ifdef LOG4CXX
static LoggerPtr logger2(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::ApplyStrain(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma) {
    TPZTensor<REAL>::TPZDecomposed sigmaDecomposition; 
    TPZTensor<REAL> sigTrial;

    TPZTensor<REAL> epsTrial = epsTotal - fPlasticState.fEpsP;

    // Compute and Decomposition of SigTrial
    fElasticResponse.Compute(epsTrial, sigTrial); // sigma = lambda Tr(E)I + 2 mu E
    sigTrial.EigenSystemJacobi(sigmaDecomposition);
    TPZManVector<REAL, 3> sigmaTrialEigenValues(sigmaDecomposition.fEigenvalues), sigmaProjectedEigenValues(3, 0.);

    // ReturnMap in the principal values
    STATE nextAlpha = -6378.;
    fYieldCriterion.ProjectSigma(sigmaTrialEigenValues, fPlasticState.fAlpha, sigmaProjectedEigenValues, nextAlpha);
    fPlasticState.fAlpha = nextAlpha;
#ifdef LOG4CXX_KEEP
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Sig Trial " << sigmaTrialEigenValues << "\nSig Project " << sigmaProjectedEigenValues << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    
    
    
    TPZTensor<REAL> newElasticDeformation;
    // Reconstruction of sigmaprTensor
    sigmaDecomposition.fEigenvalues = sigmaProjectedEigenValues; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
    sigma = TPZTensor<REAL>(sigmaDecomposition);

    fElasticResponse.ComputeDeformation(sigma, newElasticDeformation);
    fPlasticState.fEpsT = epsTotal;
    fPlasticState.fEpsP = epsTotal - newElasticDeformation;
}
template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma) {
    TPZTensor<REAL>::TPZDecomposed DecompSig; // It may be SigTr or SigPr Decomposition, it depends on the part of this method
    TPZTensor<REAL> sigtr;

    TPZTensor<REAL> epsTr(epsTotal - fPlasticState.fEpsP);

    // Compute and Decomposition of SigTrial
    fElasticResponse.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
    sigtr.EigenSystemJacobi(DecompSig);
    TPZManVector<REAL, 3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3, 0.);

    // ReturnMap in the principal values
    STATE nextalpha = -6378.;
    fYieldCriterion.ProjectSigma(sigtrvec, fPlasticState.fAlpha, sigprvec, nextalpha);
    fPlasticState.fAlpha = nextalpha;
#ifdef LOG4CXX_KEEP
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // Reconstruction of sigmaprTensor
    DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
    sigma = TPZTensor<REAL>(DecompSig);

    TPZTensor<REAL> epsElaNp1;
    fElasticResponse.ComputeDeformation(sigma, epsElaNp1);
    fPlasticState.fEpsT = epsTotal;
    fPlasticState.fEpsP = epsTotal - epsElaNp1; // Transforma epsPN em epsPNp1
}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) {
    TPZTensor<REAL>::TPZDecomposed DecompSig, DecompEps; // It may be SigTr or SigPr Decomposition, it depends on the part of this method
    TPZTensor<REAL> sigtr;

    TPZTensor<REAL> epsTr, epsPN, epsElaNp1;
    epsPN = fPlasticState.fEpsP;
    epsTr = epsTotal;
    epsTr -= epsPN; // Porque soh tem implementado o operator -=

    // Compute and Decomposition of SigTrial
    fElasticResponse.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
    epsTr.EigenSystemJacobi(DecompEps);
    sigtr.EigenSystemJacobi(DecompSig);

    TPZManVector<REAL, 3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3, 0.);

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        DecompSig.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // ReturnMap in the principal values
    STATE nextalpha = -6378.;
    STATE printPlastic = fPlasticState.Alpha();
    TPZFNMatrix<9> GradSigma(3, 3, 0.);
    fYieldCriterion.ProjectSigmaDep(sigtrvec, fPlasticState.fAlpha, sigprvec, nextalpha, GradSigma);
    //GradSigma.Print("Grad");
    fPlasticState.fAlpha = nextalpha;

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        GradSigma.Print("GradSigma", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // Aqui calculo minha matriz tangente ------------------------------------
    // Criando matriz tangente
    TPZFNMatrix<36> dSigDe(6, 6, 0.);

    //Montando a matriz tangente
    int kival[] = {0, 0, 0, 1, 1, 2};
    int kjval[] = {0, 1, 2, 1, 2, 2};
    REAL G = fElasticResponse.G();
    REAL lambda = fElasticResponse.Lambda();
    // Coluna da matriz tangente
    for (unsigned int k = 0; k < 6; ++k) {
        const unsigned int ki = kival[k];
        const unsigned int kj = kjval[k];
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                REAL temp = 2 * G * DecompSig.fEigenvectors[j][kj] * DecompSig.fEigenvectors[j][ki];

                if (ki == kj) { // diagonal
                    temp += lambda;
                } else { // off-diagonal
                    temp *= 2.;
                }
                for (int l = 0; l < 6; ++l) {
                    const unsigned int li = kival[l];
                    const unsigned int lj = kjval[l];
                    dSigDe(l, k) += temp * GradSigma(i, j) * DecompSig.fEigenvectors[i][li] * DecompSig.fEigenvectors[i][lj];
                }/// l
            }///j
        }///i
    }///k

    REAL deigensig, deigeneps;
    TPZFNMatrix<36> RotCorrection(6, 6, 0.);
    // Correcao do giro rigido
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = i + 1; j < 3; ++j) {
            deigeneps = DecompEps.fEigenvalues[i] - DecompEps.fEigenvalues[j];
            deigensig = sigprvec[i] - sigprvec[j];
            TPZFMatrix<REAL> tempMat(3, 3, 0.);
            REAL factor = 0.;
            if (!IsZero(deigeneps)) {
                factor = deigensig / deigeneps;
            } else {
                factor = fElasticResponse.G() * (GradSigma(i, i) - GradSigma(i, j) - GradSigma(j, i) + GradSigma(j, j));
            }
            tempMat = TPZFMatrix<REAL>(DecompEps.fEigenvectors[i], DecompEps.fEigenvectors[j]) + TPZFMatrix<REAL>(DecompEps.fEigenvectors[j], DecompEps.fEigenvectors[i]);
            for (unsigned int k = 0; k < 6; ++k) {
                const unsigned int ki = kival[k];
                const unsigned int kj = kjval[k];
                TPZFMatrix<REAL> ColCorr(3, 3, 0.);
                TPZFMatrix<REAL> ColCorrV(6, 1, 0.);
                if (ki == kj) { // diagonal
                    ColCorr = (DecompEps.fEigenvectors[j][ki] * DecompEps.fEigenvectors[i][kj]) * factor * tempMat;
                } else { // off-diagonal
                    ColCorr = (DecompEps.fEigenvectors[j][ki] * DecompEps.fEigenvectors[i][kj] + DecompEps.fEigenvectors[j][kj] * DecompEps.fEigenvectors[i][ki]) * factor * tempMat;
                }
                ColCorrV = FromMatToVoight(ColCorr);
                for (int l = 0; l < 6; l++) {
                    RotCorrection(l, k) += ColCorrV(l, 0);
                }
            }
        } // j
    } // i

    dSigDe += RotCorrection;

#ifdef LOG4CXX
    {
        if (logger->isDebugEnabled()) {
            std::stringstream str;
            str << "\n**********************MATRIZ TANGENTE**********************" << endl;
            dSigDe.Print("Matriz Tangente:", str);
            str << "\n**********************CORRECAO GIRO**********************" << endl;
            RotCorrection.Print("GiroCorrection", str);
            LOGPZ_DEBUG(logger, str.str())
        }
    }
#endif

    // Reconstruction of sigmaprTensor

    DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
    sigma = TPZTensor<REAL>(DecompSig);

    fElasticResponse.ComputeDeformation(sigma, epsElaNp1);
    fPlasticState.fEpsT = epsTotal;
    epsPN = epsTotal;
    epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
    fPlasticState.fEpsP = epsPN;
    Dep = dSigDe;


#ifdef LOG4CXX
    if (logger2->isDebugEnabled()) {
        if (fabs(printPlastic - fPlasticState.fAlpha) > 1.e-4) {
            std::stringstream sout;
            TPZVec<STATE> phi;
            TPZTensor<STATE> epsElastic(fPlasticState.fEpsT);
            epsElastic -= fPlasticState.fEpsP;
            Phi(epsElastic, phi);
            sout << " \n phi = [";
            for (int i = 0; i < phi.size(); i++) {
                sout << phi[i] << " ";
            }

            sout << " ] " << endl;

            sout << " \n eigenvalues Sigma = [";
            for (int i = 0; i < 3; i++) {
                sout << DecompSig.fEigenvalues[i] << " ";
            }

            sout << " ] " << endl;



            LOGPZ_DEBUG(logger2, sout.str())
        }
    }
#endif
}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv) {
    TPZTensor<REAL> eps1, eps2, SigmaTemp, Sigma1, Sigma2;
    TPZFNMatrix <36> dSigDe(6, 6, 0.);
    TPZStack<REAL> coef;

    fPlasticState.fEpsP.Scale(0.);
    fPlasticState.fEpsT.Scale(0.);
    fPlasticState.fAlpha = kprev;
    this->ApplyStrainComputeDep(EpsIni, SigmaTemp, dSigDe);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "EpsIni " << EpsIni << "\nSigmaTemp " << SigmaTemp << "\ndSidDe " << dSigDe << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fPlasticState.fEpsP.Scale(0.);
    fPlasticState.fEpsT.Scale(0.);
    fPlasticState.fAlpha = kprev;

    REAL scale = 1.;
    REAL alphatable[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (int i = 0; i < 6; i++) {
        alphatable[i] *= scale;
    }
    for (int ia = 0; ia < 5; ia++) {
        REAL alpha1 = alphatable[0];
        REAL alpha2 = alphatable[ia + 1];
        eps1.Scale(0.);
        eps2.Scale(0.);
        eps1 = EpsIni;
        eps2 = EpsIni;
        eps1.Add(deps, alpha1);
        eps2.Add(deps, alpha2);

        fPlasticState.fEpsT = EpsIni;
        this->ApplyStrainComputeSigma(eps1, Sigma1);
        fPlasticState.fEpsP.Scale(0.);
        fPlasticState.fEpsT.Scale(0.);
        fPlasticState.fAlpha = kprev;

        fPlasticState.fEpsT = EpsIni;
        this->ApplyStrainComputeSigma(eps2, Sigma2);
        fPlasticState.fEpsP.Scale(0.);
        fPlasticState.fEpsT.Scale(0.);
        fPlasticState.fAlpha = kprev;

        TPZFNMatrix <6> deps1(6, 1, 0.), deps2(6, 1, 0.);
        TPZFNMatrix <9> depsMat(3, 3, 0.);
        depsMat = deps;
        deps1 = FromMatToVoight(depsMat);
        deps2 = FromMatToVoight(depsMat);

        TPZFNMatrix <6> tanmult1(6, 1, 0.), tanmult2(6, 1, 0.);
        dSigDe.Multiply(deps1, tanmult1);
        dSigDe.Multiply(deps2, tanmult2);

        for (int i = 0; i < 6; i++) {
            tanmult1(i, 0) *= alpha1;
            tanmult2(i, 0) *= alpha2;
        }

        TPZFNMatrix <9> SigMatTemp33(3, 3, 0.);
        TPZFNMatrix <6> sigprMat(6, 1, 0.), sigpr1Mat(6, 1, 0.), sigpr2Mat(6, 1, 0.);
        SigMatTemp33 = SigmaTemp;
        sigprMat = FromMatToVoight(SigMatTemp33);
        SigMatTemp33 = Sigma1;
        sigpr1Mat = FromMatToVoight(SigMatTemp33);
        SigMatTemp33 = Sigma2;
        sigpr2Mat = FromMatToVoight(SigMatTemp33);

        TPZFNMatrix<6> error1(6, 1, 0.), error2(6, 1, 0.);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sigprMat.Print("sigprMat", sout);
            sigpr1Mat.Print("sigpr1Mat", sout);
            tanmult1.Print("tanmult1", sout);
            sigpr2Mat.Print("sigpr2Mat", sout);
            tanmult2.Print("tanmult2", sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i = 0; i < 6; i++) {
            error1(i, 0) = sigpr1Mat(i, 0) - sigprMat(i, 0) - tanmult1(i, 0);
            error2(i, 0) = sigpr2Mat(i, 0) - sigprMat(i, 0) - tanmult2(i, 0);
        }
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            error1.Print("error1:", sout);
            error2.Print("error2:", sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        REAL n;
        REAL norm1, norm2;
        norm1 = Norm(error1);
        norm2 = Norm(error2);
        n = (log(norm1) - log(norm2)) / (log(alpha1) - log(alpha2));
        coef.push_back(n);
    }
    conv = coef;
    std::cout << "coef = " << coef << std::endl;
}

template <class YieldCriterion_type, class ElasticResponse_type >
REAL TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat) {
    REAL norm1, norm2, n;
    norm1 = Norm(error1Mat);
    norm2 = Norm(error2Mat);
    n = log(norm1 / norm2) / log(alpha1 / alpha2);
    return n;
}

TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat) {
    TPZFNMatrix <6> voi(6, 1, 0.);
    int k = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            voi(k++, 0) = mat(i, j);
        }
    }
    return voi;
}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::ApplyStrain(const TPZTensor<REAL> &epsTotal) {

    std::cout << " \n this method is not implemented in PlasticStepPV. ";
    DebugStop();

}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::ApplyLoad(const TPZTensor<REAL> & GivenStress, TPZTensor<REAL> &epsTotal) {
    TPZPlasticState<STATE> prevstate = GetState();
    epsTotal = prevstate.fEpsP;
    TPZTensor<STATE> GuessStress, Diff, Diff2, deps;
    TPZFNMatrix<36, STATE> Dep(6, 6);
    TPZFNMatrix<6, STATE> DiffFN(6, 1);

    ApplyStrainComputeDep(epsTotal, GuessStress, Dep);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        Dep.Print("Dep = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Diff = GivenStress;
    Diff -= GuessStress;

    STATE norm = Norm(Diff), normprev;
    STATE tol = 1.e-7;
    int counter = 0;

    while (norm > tol && counter < 30) {
        CopyFromTensorToFMatrix(Diff, DiffFN);
        std::list<long> singular;
        Dep.Solve_LU(&DiffFN, singular);
        CopyFromFMatrixToTensor(DiffFN, Diff);
        TPZTensor<STATE> epsprev(epsTotal);
        normprev = norm;
        STATE scale = 1.;
        int counter2 = 0;
        do {
            for (unsigned int i = 0; i < 6; ++i)epsTotal.fData[i] = epsprev.fData[i] + scale * Diff.fData[i];

            ApplyStrainComputeDep(epsTotal, GuessStress, Dep);
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                Dep.Print("Dep = ", sout, EMathematicaInput);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif

            fPlasticState = prevstate;
            Diff2 = GivenStress;
            Diff2 -= GuessStress;
            CopyFromTensorToFMatrix(Diff2, DiffFN);
            Dep.Solve_LU(&DiffFN, singular);
            norm = Norm(Diff2);
            //scale*=0.5;
            counter2++;
        } while (norm >= normprev && counter2 < 30);
        Diff = Diff2;
        counter++;
    }
    ApplyStrainComputeDep(epsTotal, GuessStress, Dep);

}

template <class YieldCriterion_type, class ElasticResponse_type >
TPZPlasticState<STATE> TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::GetState() const {
    return fPlasticState;
}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const {
    TPZTensor<STATE> sigma;
    fElasticResponse.Compute(eps, sigma);
    TPZVec<STATE> sigvec(3);
    sigma.EigenValuesJacobi(sigvec);
    fYieldCriterion.Phi(sigvec, fPlasticState.Alpha(), phi);
}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::SetState(const TPZPlasticState<REAL> &state) {
    fPlasticState = state;
}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::Read(TPZStream &buf) {
    fYieldCriterion.Read(buf);
    fElasticResponse.Read(buf);
    buf.Read(&fResidualTol);
    buf.Read(&fMaxNewtonIterations);
    fPlasticState.Read(buf);
}

/** @brief Object which represents the yield criterion */
//YieldCriterion_type fYC;

/** @brief Object representing the elastic response */
//ElasticResponse_type fER;

/** @brief Residual tolerance accepted in the plastic loop processes */
//REAL fResTol;

/** @brief Maximum number of Newton interations allowed in the nonlinear solvers */
//int fMaxNewton;	// COLOCAR = 30 (sugestao do erick!)




/** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */
//TPZPlasticState<STATE> fN;

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::Write(TPZStream &buf) const {
    fYieldCriterion.Write(buf);
    fElasticResponse.Write(buf);
    buf.Write(&fResidualTol);
    buf.Write(&fMaxNewtonIterations);
    fPlasticState.Write(buf);

}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM, TPZTensor<STATE> &copy) const {
    FNM.Resize(6, 1);
    copy.XX() = FNM(0, 0);
    copy.XY() = FNM(1, 0);
    copy.XZ() = FNM(2, 0);
    copy.YY() = FNM(3, 0);
    copy.YZ() = FNM(4, 0);
    copy.ZZ() = FNM(5, 0);
}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::CopyFromTensorToFMatrix(TPZTensor<STATE> tensor, TPZFMatrix<STATE> &copy) const {
    copy(0, 0) = tensor.XX();
    copy(1, 0) = tensor.XY();
    copy(2, 0) = tensor.XZ();
    copy(3, 0) = tensor.YY();
    copy(4, 0) = tensor.YZ();
    copy(5, 0) = tensor.ZZ();
}

template <class YieldCriterion_type, class ElasticResponse_type>
void TPZPlasticStepPV<YieldCriterion_type, ElasticResponse_type>::SetElasticResponse(TPZElasticResponse &elasticResponse) {
    fElasticResponse = elasticResponse;
    fYieldCriterion.SetElasticResponse(elasticResponse);
}

template class TPZPlasticStepPV<TPZSandlerExtendedPV>;
template class TPZPlasticStepPV<TPZYCMohrCoulombPV>;

/*
 // Correcao do giro rigido
 for (int i = 0; i < 2; i++) {
 for (int j = i+1; j<3 ; j++) {
 deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
 deigensig = sigprvec[i] - sigprvec[j];
 TPZFNMatrix<9,REAL> tempMat(3,1,0.);
 depsMat.Multiply(epsegveFromProj[i], tempMat);
 REAL deij = InnerVecOfMat(tempMat,epsegveFromProj[j]);
 REAL factor = 0.;
 if (!IsZero(deigeneps)) {
 factor = deigensig * deij / deigeneps;
 }
 else {
 factor = fER.G() * ( GradSigma(i,i) - GradSigma(i,j) - GradSigma(j,i) + GradSigma(j,j) ) * deij;
 }
 std::cout << "factor = " << factor << std::endl;
 std::cout << "G = " << fER.G() << std::endl;
 GradSigma.Print("GradSigma");
 tempMat.Redim(3, 3);
 tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
 factorMat += tempMat * factor;
 
 }
 }
 */