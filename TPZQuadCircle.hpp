/**
 * @file
 * @brief Contains the TPZQuadCircle class. It is a special map.
 */

#ifndef TPZQuadCircle_hpp
#define TPZQuadCircle_hpp

#include <stdio.h>
#include "pzgeoquad.h"


namespace pzgeom {

    template<class TGeo = pzgeom::TPZGeoQuad>
    class TPZQuadCircle : public TGeo {
        
    private:
        
        /** @brief radius of the circle */
        REAL qR;
        
        /** @brief coordinates of the center of the circle*/
        TPZVec<REAL> qxc;
        
    protected:
        
    public:
        
        //contrutores da classe Quad
        //write e read tpzgeoquad
        
        /** @brief Constructor with list of nodes */
        
        TPZQuadCircle(TPZVec<long> &nodeindexes) : TGeo(nodeindexes)
        {
        }
        
        /** @brief Empty constructor */
        
        TPZQuadCircle() : TGeo()
        {
        }
        
        /** @brief Constructor with node map */
        
        TPZQuadCircle(const TPZQuadCircle &cp, std::map<long,long> & gl2lcNdMap) : TGeo(cp,gl2lcNdMap)
        {
        }
        
        /** @brief Copy constructor */
        
        TPZQuadCircle(const TPZQuadCircle &cp) : TGeo(cp)
        {
        }
        
        /** @brief Copy constructor */
        
        TPZQuadCircle(const TPZQuadCircle &cp, TPZGeoMesh &) : TGeo(cp)
        {
        }
        
        
        virtual void Write(TPZStream &buf, int withclassid) {
            
            // TPZGeoQuad não possui metodo write
            // TGeo::Write(buf,withclassid);
            
        }
        
        virtual void Read(TPZStream &buf,void *context) {
            
            // TPZGeoQuad não possui metodo read
            // TGeo::Read(buf,context);
            
        }
        
        virtual void SetFather(long fatherindex){
            
        }
        
        virtual void Father() {
            
        }
        
        virtual void GradX(TPZGeoEl &el,TPZVec<REAL> &qsi, TPZFMatrix<REAL> &gradx) const
        {
            TGeo::GradX(el,qsi,gradx);
            
            
        }
        
        virtual void X(TPZGeoEl &el,TPZVec<REAL> &ksi,TPZVec<REAL> &result) const
        {
            TGeo::X(el,ksi,result);
            
        }
        
        virtual void Print (std::ostream &out = std::cout) {
            
            
            
            
        }
        
    };
}



#endif 
