#include "TPZMatHelmholtz2D.h"
#include "TPZVecL2.h"

#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.fran"));
#endif

TPZMatHelmholtz2D::TPZMatHelmholtz2D(int id,
                                     STATE (&cFunc)(const TPZVec<REAL> &))
    : TPZVecL2(id), fC(cFunc) {
    fDim = 2;
}

TPZMatHelmholtz2D::TPZMatHelmholtz2D(int id) : TPZVecL2(id), fC(urDefault) {
    fDim = 2;
}

/** @brief Default constructor */
TPZMatHelmholtz2D::TPZMatHelmholtz2D() : TPZVecL2(), fC(urDefault) { fDim = 2; }

TPZMatHelmholtz2D::TPZMatHelmholtz2D(const TPZMatHelmholtz2D &mat)
    : TPZVecL2(mat), fC(mat.fC) {
    fDim = mat.fDim;
}

TPZMatHelmholtz2D::~TPZMatHelmholtz2D() {}

void TPZMatHelmholtz2D::Contribute(TPZMaterialData &data, REAL weight,
                                   TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) {
    /*********************CREATE HCURL FUNCTIONS*****************************/
    TPZFNMatrix<36, REAL> phiHCurlAxes = data.phi;
    TPZFNMatrix<40, REAL> curlPhiDAxes = data.dphix;

    TPZFNMatrix<40, REAL> curlPhi, phiHCurl;

    TPZAxesTools<REAL>::Axes2XYZ(phiHCurlAxes, phiHCurl, data.axes, false);

    TPZManVector<REAL, 3> ax1(3), ax2(3), elNormal(3);
    for (int i = 0; i < 3; i++) {
        ax1[i] = data.axes(0, i);
        ax2[i] = data.axes(1, i);
    }
    Cross(ax1, ax2, elNormal);
    TPZFNMatrix<3, REAL> normalVec(1, 3);
    normalVec(0, 0) = elNormal[0];
    normalVec(0, 1) = elNormal[1];
    normalVec(0, 2) = elNormal[2];
    TPZAxesTools<REAL>::Axes2XYZ(curlPhiDAxes, curlPhi, normalVec);

    TPZManVector<REAL, 3> x = data.x;
    const STATE cVal = fC(x);

    //*****************GET FORCING FUNCTION****************//
    TPZManVector<STATE, 3> force(3);
    force.Fill(0.);
    if (fForcingFunction) {
        fForcingFunction->Execute(data.x, force);
    } else {
        DebugStop(); // RHS not set!
    }
    //*****************ACTUAL COMPUTATION OF CONTRIBUTION****************//

    const int nHCurlFunctions = phiHCurl.Rows();

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        STATE load = 0.;
        load += phiHCurl(iVec, 0) * force[0];
        load += phiHCurl(iVec, 1) * force[1];
        load += phiHCurl(iVec, 2) * force[2];
        ef(iVec) += load * weight;
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE stiff = 0.;

            STATE phiIdotPhiJ = 0.;
            phiIdotPhiJ += phiHCurl(iVec, 0) * phiHCurl(jVec, 0);
            phiIdotPhiJ += phiHCurl(iVec, 1) * phiHCurl(jVec, 1);
            phiIdotPhiJ += phiHCurl(iVec, 2) * phiHCurl(jVec, 2);

            STATE curlPhiIvecCurlPhiJ = 0.;
            curlPhiIvecCurlPhiJ += curlPhi(0, iVec) * curlPhi(0, jVec);
            curlPhiIvecCurlPhiJ += curlPhi(1, iVec) * curlPhi(1, jVec);
            curlPhiIvecCurlPhiJ += curlPhi(2, iVec) * curlPhi(2, jVec);

            stiff = curlPhiIvecCurlPhiJ + cVal * phiIdotPhiJ;

            ek(iVec, jVec) += stiff * weight;
        }
    }
}

void TPZMatHelmholtz2D::ContributeBC(TPZMaterialData &data, REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    TPZFMatrix<REAL> &phiHCurl = data.phi;

    int nHCurlFunctions = phiHCurl.Rows();
    REAL BIG = TPZMaterial::gBigNumber;

    // const STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de
    // condicao mista
    const STATE v2 = bc.Val2()(0, 0); // sera posto no vetor F

    switch (bc.Type()) {
    case 0:
        for (int i = 0; i < nHCurlFunctions; i++) {
            const STATE rhs = phiHCurl(i, 0) * BIG * v2;
            ef(i, 0) += rhs * weight;
            for (int j = 0; j < nHCurlFunctions; j++) {
                const STATE stiff = phiHCurl(i, 0) * phiHCurl(j, 0) * BIG;
                ek(i, j) += stiff * weight;
            }
        }
        break;
    case 1:
        DebugStop();
        break;
    case 2:
        DebugStop();
        break;
    }
}

int TPZMatHelmholtz2D::IntegrationRuleOrder(int elPMaxOrder) const {
    return 2 + elPMaxOrder * 2;
}

int TPZMatHelmholtz2D::VariableIndex(const std::string &name) {
    if (strcmp(name.c_str(), "E") == 0)
        return 0;
    DebugStop();
    return 1;
}

/**
 * @brief Returns the number of variables associated with the variable indexed
 * by var.
 * @param var Index variable into the solution, is obtained by calling
 * VariableIndex
 */
int TPZMatHelmholtz2D::NSolutionVariables(int var) {
    switch (var) {
    case 0: // Et
        return 2;
        break;
    default:
        DebugStop();
        break;
    }
    return 1;
}

/** @brief Returns the solution associated with the var index based on the
 * finite element approximation */
void TPZMatHelmholtz2D::Solution(TPZMaterialData &data, int var,
                                 TPZVec<STATE> &Solout) {
    switch (var) {
    case 0: // E
    {
        Solout = data.sol[0];
    } break;
    default:
        DebugStop();
        break;
    }
}

void TPZMatHelmholtz2D::Errors(
    TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &curlU,
    TPZFMatrix<REAL> &axes, TPZVec<STATE> & /*flux*/, TPZVec<STATE> &u_exact,
    TPZFMatrix<STATE> &curlU_exact,
    TPZVec<REAL> &values) { // TODO:Implementar direito

    values.Resize(NEvalErrors());
    values.Fill(0.0);

    TPZManVector<STATE> sol(1), dsol(3, 0.);
    int id;
    // values[0] : erro em norma H1 <=> norma Energia
    // values[1] : eror em norma L2
    // values[2] : erro em semi norma H1
    // values[2] : erro em norma HCurl
    values[0] = 0;
    REAL diff = 0.;
    for (id = 0; id < fDim; id++) {
        diff =
            std::real((u[id] - u_exact[id]) * std::conj(u[id] - u_exact[id]));
        values[0] += diff;
    }

    values[1] = 0.;
    for (id = 0; id < fDim; id++) {
        diff = std::real((curlU(id, 0) - curlU_exact(id, 0)) *
                         std::conj(curlU(id, 0) - curlU_exact(id, 0)));
        values[1] += diff;
    }
}
