/**
* @file
* @brief Contains the implementation of the TPZHCurlNedFLinEl::Shape method.
*/

#include "TPZHCurlNedFLinEl.h"
#ifdef HCURL_HIERARCHICAL
#include "pzshapelinear.h"

using namespace pzshape;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZHCurlNedFLinEl"));
#endif

void TPZHCurlNedFLinEl::Shape(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi,
                             TPZFMatrix<REAL> &curlPhiHat) {
 const int nCon = NConnects();
 const int dim = Dimension();
 const int firstSide = TPZShapeLinear::NSides - TPZShapeLinear::NFaces - 1;

#ifdef PZDEBUG
 if (!(nCon == 1 && dim == 1 && firstSide == 0)) {
   DebugStop();
 }
#endif

 const int lastFuncPos = NConnectShapeF(nCon - 1, ConnectOrder(nCon - 1)) - 1;

 phi.Resize(lastFuncPos + 1, dim);
 curlPhiHat.Resize(1, 1);  //The curl wont be calculated in boundary
							//elements for now.
 const int pOrder = ConnectOrder(0);
 int currentFuncPos = lastFuncPos;

 switch (pOrder) {
   case 15:
     phi(currentFuncPos, 0) =
         (1163381400 *
          (1 - 210 * qsi[0] + 10920 * pow(qsi[0], 2) -
           247520 * pow(qsi[0], 3) + 3063060 * pow(qsi[0], 4) -
           23279256 * pow(qsi[0], 5) + 116396280 * pow(qsi[0], 6) -
           399072960 * pow(qsi[0], 7) + 960269310 * pow(qsi[0], 8) -
           1636014380 * pow(qsi[0], 9) + 1963217256 * pow(qsi[0], 10) -
           1622493600 * pow(qsi[0], 11) + 878850700 * pow(qsi[0], 12) -
           280816200 * pow(qsi[0], 13) + 40116600 * pow(qsi[0], 14)));
     currentFuncPos--;
   case 14:
     phi(currentFuncPos, 0) =
          (280816200 *
                   (-1 + 182 * qsi[0] - 8190 * pow(qsi[0], 2) +
                    160160 * pow(qsi[0], 3) - 1701700 * pow(qsi[0], 4) +
                    11027016 * pow(qsi[0], 5) - 46558512 * pow(qsi[0], 6) +
                    133024320 * pow(qsi[0], 7) - 261891630 * pow(qsi[0], 8) +
                    355655300 * pow(qsi[0], 9) - 327202876 * pow(qsi[0], 10) +
                    194699232 * pow(qsi[0], 11) - 67603900 * pow(qsi[0], 12) +
                    10400600 * pow(qsi[0], 13)));
     currentFuncPos--;
   case 13:
     phi(currentFuncPos, 0) =
         
         (67603900 * (1 - 156 * qsi[0] + 6006 * pow(qsi[0], 2) -
                      100100 * pow(qsi[0], 3) + 900900 * pow(qsi[0], 4) -
                      4900896 * pow(qsi[0], 5) + 17153136 * pow(qsi[0], 6) -
                      39907296 * pow(qsi[0], 7) + 62355150 * pow(qsi[0], 8) -
                      64664600 * pow(qsi[0], 9) + 42678636 * pow(qsi[0], 10) -
                      16224936 * pow(qsi[0], 11) + 2704156 * pow(qsi[0], 12)));
     currentFuncPos--;
   case 12:
     phi(currentFuncPos, 0) =
         
         (16224936 *
          (-1 + 132 * qsi[0] - 4290 * pow(qsi[0], 2) + 60060 * pow(qsi[0], 3) -
           450450 * pow(qsi[0], 4) + 2018016 * pow(qsi[0], 5) -
           5717712 * pow(qsi[0], 6) + 10501920 * pow(qsi[0], 7) -
           12471030 * pow(qsi[0], 8) + 9237800 * pow(qsi[0], 9) -
           3879876 * pow(qsi[0], 10) + 705432 * pow(qsi[0], 11)));
     currentFuncPos--;
   case 11:
     phi(currentFuncPos, 0) =
         
         (3879876 * (1 - 110 * qsi[0] + 2970 * pow(qsi[0], 2) -
                     34320 * pow(qsi[0], 3) + 210210 * pow(qsi[0], 4) -
                     756756 * pow(qsi[0], 5) + 1681680 * pow(qsi[0], 6) -
                     2333760 * pow(qsi[0], 7) + 1969110 * pow(qsi[0], 8) -
                     923780 * pow(qsi[0], 9) + 184756 * pow(qsi[0], 10)));
     currentFuncPos--;
   case 10:
     phi(currentFuncPos, 0) =
         
         (923780 *
          (-1 + 90 * qsi[0] - 1980 * pow(qsi[0], 2) + 18480 * pow(qsi[0], 3) -
           90090 * pow(qsi[0], 4) + 252252 * pow(qsi[0], 5) -
           420420 * pow(qsi[0], 6) + 411840 * pow(qsi[0], 7) -
           218790 * pow(qsi[0], 8) + 48620 * pow(qsi[0], 9)));
     currentFuncPos--;
   case 9:
     phi(currentFuncPos, 0) =
          (218790 * (1 - 72 * qsi[0] + 1260 * pow(qsi[0], 2) -
                             9240 * pow(qsi[0], 3) + 34650 * pow(qsi[0], 4) -
                             72072 * pow(qsi[0], 5) + 84084 * pow(qsi[0], 6) -
                             51480 * pow(qsi[0], 7) + 12870 * pow(qsi[0], 8)));
     currentFuncPos--;
   case 8:
     phi(currentFuncPos, 0) =
          (51480 * (-1 + 56 * qsi[0] - 756 * pow(qsi[0], 2) +
                            4200 * pow(qsi[0], 3) - 11550 * pow(qsi[0], 4) +
                            16632 * pow(qsi[0], 5) - 12012 * pow(qsi[0], 6) +
                            3432 * pow(qsi[0], 7)));
     currentFuncPos--;
   case 7:
     phi(currentFuncPos, 0) =
          (12012 * (1 - 42 * qsi[0] + 420 * pow(qsi[0], 2) -
                            1680 * pow(qsi[0], 3) + 3150 * pow(qsi[0], 4) -
                            2772 * pow(qsi[0], 5) + 924 * pow(qsi[0], 6)));
     currentFuncPos--;
   case 6:
     phi(currentFuncPos, 0) =
          (2772 * (-1 + 30 * qsi[0] - 210 * pow(qsi[0], 2) +
                           560 * pow(qsi[0], 3) - 630 * pow(qsi[0], 4) +
                           252 * pow(qsi[0], 5)));
     currentFuncPos--;
   case 5:
     phi(currentFuncPos, 0) =
          (630 * (1 - 20 * qsi[0] + 90 * pow(qsi[0], 2) -
                          140 * pow(qsi[0], 3) + 70 * pow(qsi[0], 4)));
     currentFuncPos--;
   case 4:
     phi(currentFuncPos, 0) =
          (140 * (-1 + 12 * qsi[0] - 30 * pow(qsi[0], 2) +
                          20 * pow(qsi[0], 3)));
     currentFuncPos--;
   case 3:
     phi(currentFuncPos, 0) =
          (30 * (1 - 6 * qsi[0] + 6 * pow(qsi[0], 2)));
     currentFuncPos--;
   case 2:
     phi(currentFuncPos, 0) =  (-6 + 12 * qsi[0]);
     currentFuncPos--;
   case 1:
     phi(currentFuncPos, 0) =  (1);
     break;
   default:
     DebugStop();  // polynomial order not implemented!
 }
}
#endif
