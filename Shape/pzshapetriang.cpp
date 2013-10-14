/**
 * @file
 * @brief Contains the implementation of the TPZShapeTriang methods.
 */
// $Id: pzshapetriang.cpp,v 1.9 2008-04-08 20:10:41 fortiago Exp $
#include "pzshapetriang.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"

using namespace std;

namespace pzshape {
	
	/**Transformation of the point within a triangular face */
	REAL TPZShapeTriang::gTrans2dT[6][2][2] = {//s* , t*
		{ { 1., 0.},{ 0., 1.} },
		{ { 0., 1.},{ 1., 0.} },
		{ { 0., 1.},{-1.,-1.} },//s* = t   t* = -s-t-1 ,  etc
		{ {-1.,-1.},{ 0., 1.} },
		{ {-1.,-1.},{ 1., 0.} },
		{ { 1., 0.},{-1.,-1.} }
	};

	REAL TPZShapeTriang::gVet2dT[6][2] = {  {0.,0.},{0.,0.},{0.,1.},{1.,0.},{1.,0.},{0.,1.} };

	REAL TPZShapeTriang::gRibTrans2dT1d[3][2] = { {2.,1.},{-1.,1.},{-1.,-2.} };//Cedric : 06/03/99

	REAL TPZShapeTriang::gVet1dT[3] = {-1.,0.,1.};

	void TPZShapeTriang::ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi)
    {
		phi(0,0) =  1.-pt[0]-pt[1];
		phi(1,0) =  pt[0];
		phi(2,0) =  pt[1];

        // first edge function
		dphi(0,0) = -1.;    // dx
		dphi(1,0) = -1.;    // dy

        // second edge function
		dphi(0,1) =  1.;    // dx
		dphi(1,1) =  0.;    // dy

        // third edge function
		dphi(0,2) =  0.;    // dx
		dphi(1,2) =  1.;    // dy

        if( dphi.Rows() > 2 )
        {
            dphi(2,0) =  0.;    // dxdx
            dphi(3,0) =  0.;    // dydx
            dphi(4,0) =  0.;    // dxdy
            dphi(5,0) =  0.;    // dydy

            dphi(2,1) =  0.;    // dxdx
            dphi(3,1) =  0.;    // dydx
            dphi(4,1) =  0.;    // dxdy
            dphi(5,1) =  0.;    // dydy

            dphi(2,2) =  0.;    // dxdx
            dphi(3,2) =  0.;    // dydx
            dphi(4,2) =  0.;    // dxdy
            dphi(5,2) =  0.;    // dydy
        }
	}

	/**
	 * Computes the generating shape functions for a quadrilateral element
	 * @param pt (input) point where the shape function is computed
	 * @param phi (input/output) value of the (4) shape functions
	 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
	 */
	void TPZShapeTriang::ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi)
	{
		int is;
		for(is=3; is<6; is++)
		{
			int is1 = is%3;
			int is2 = (is+1)%3;

			phi(is,0) = phi(is1,0)*phi(is2,0);

			dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2); // dx(phi[is1]*phi[is2])
			dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2); // dy(phi[is1]*phi[is2])

            if (dphi.Rows() > 2)
            {
                dphi(2,is) = dphi(0,is1)*dphi(0,is2) + dphi(0,is1)*dphi(0,is2);   // dxdx(phi[is1]*phi[is2])
                dphi(3,is) = dphi(0,is1)*dphi(1,is2) + dphi(1,is1)*dphi(0,is2);   // dydx(phi[is1]*phi[is2])
                dphi(4,is) = dphi(1,is1)*dphi(0,is2) + dphi(0,is1)*dphi(1,is2);   // dxdy(phi[is1]*phi[is2])
                dphi(5,is) = dphi(1,is1)*dphi(1,is2) + dphi(1,is1)*dphi(1,is2);   // dydy(phi[is1]*phi[is2])
            }
		}
		int is1 = 0;
		int is2 = 1;
		int is3 = 2;

		phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);   // phi[1]*phi[2]*phi[3]

		dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3); // dx(phi[1]*phi[2]*phi[3])

		dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3); // dy(phi[1]*phi[2]*phi[3])

        if (dphi.Rows() > 2)
        {
            dphi(2,is) = dphi(0,is1)*dphi(0,is2)*phi(is3,0) + dphi(0,is1)*phi(is2,0)*dphi(0,is3) +
                    dphi(0,is1)*dphi(0,is2)*phi(is3,0) + phi(is1,0)*dphi(0,is2)*dphi(0,is3) +
                    dphi(0,is1)*phi(is2,0)*dphi(0,is3) + phi(is1,0)*dphi(0,is2)*dphi(0,is3);  // dxdx(phi[1]*phi[2]*phi[3])

            dphi(3,is) = dphi(0,is1)*dphi(1,is2)*phi(is3,0) + dphi(0,is1)*phi(is2,0)*dphi(1,is3) +
                    dphi(1,is1)*dphi(0,is2)*phi(is3,0) + phi(is1,0)*dphi(0,is2)*dphi(1,is3) +
                    dphi(1,is1)*phi(is2,0)*dphi(0,is3) + phi(is1,0)*dphi(1,is2)*dphi(0,is3);  // dydx(phi[1]*phi[2]*phi[3])

            dphi(4,is) = dphi(1,is1)*dphi(0,is2)*phi(is3,0) + dphi(1,is1)*phi(is2,0)*dphi(0,is3) +
                    dphi(0,is1)*dphi(1,is2)*phi(is3,0) + phi(is1,0)*dphi(1,is2)*dphi(0,is3) +
                    dphi(0,is1)*phi(is2,0)*dphi(1,is3) + phi(is1,0)*dphi(0,is2)*dphi(1,is3);  // dxdy(phi[1]*phi[2]*phi[3])

            dphi(5,is) = dphi(1,is1)*dphi(1,is2)*phi(is3,0) + dphi(1,is1)*phi(is2,0)*dphi(1,is3) +
                    dphi(1,is1)*dphi(1,is2)*phi(is3,0) + phi(is1,0)*dphi(1,is2)*dphi(1,is3) +
                    dphi(1,is1)*phi(is2,0)*dphi(1,is3) + phi(is1,0)*dphi(1,is2)*dphi(1,is3);  // dydy(phi[1]*phi[2]*phi[3])
        }

		// Make the generating shape functions linear and unitary
		/* MODIFIED BY FREDERICO. SEP 19TH, 12:30. 
        REAL mult[] = {1.,1.,1.,4.,4.,4.,27.};
		for(is=3;is<NSides; is++)
		{
			phi(is,0) *= mult[is];
			dphi(0,is) *= mult[is];
			dphi(1,is) *= mult[is];
            if( dphi.Rows() > 2)
            {
                dphi(2,is) *= mult[is];
                dphi(3,is) *= mult[is];
                dphi(4,is) *= mult[is];
                dphi(5,is) *= mult[is];
            }
		}
         */
	}

	void TPZShapeTriang::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
							   TPZFMatrix &phi,TPZFMatrix &dphi)
    {
        //std::cout << "id: " << id[0] << ", " << id[1] << ", " << id[2] << std::endl;

		ShapeCorner(pt,phi,dphi);
		if (order[0] < 2 && order[1] < 2 && order[2] < 2 && order[3] < 3) return;
		int is,d;
        if (dphi.Rows() > 2)
        {
            TPZFNMatrix<100> phiblend(NSides,1),dphiblend(Dimension+Dimension*Dimension,NSides);
            for(is=0; is<NCornerNodes; is++)
            {
                phiblend(is,0) = phi(is,0);
                for(d=0; d<(Dimension+Dimension*Dimension); d++)
                {
                    dphiblend(d,is) = dphi(d,is);
                }
            }
            ShapeGenerating(pt,phiblend,dphiblend);

            REAL out;
            int shape = 3;
            for (int rib = 0; rib < 3; rib++)
            {
                ProjectPoint2dTriangToRib(rib,pt,out);
                TPZManVector<REAL,1> outvec(1,out);
                TPZVec<int> ids(2);
                ids[0] = id[rib%3];
                ids[1] = id[(rib+1)%3];
                REAL store1[20],store2[40];
                int ord2 = order[rib]-1;//numero de shapes por lado rib
                TPZFMatrix phin(ord2,1,store1,20),dphin(2,ord2,store2,40);
                TPZShapeLinear *shplin=0;
                REAL y;
                // this operation below is performed inside "TPZShapeLinear::ShapeInternal"
                int transf = shplin->GetTransformId1d(ids);
                shplin->TransformPoint1d( transf, outvec[0], y);
                outvec[0] = y;
                shplin->ShapeInternal( outvec, order[rib], phin, dphin, transf);
                shplin->TransformDerivative1d( transf, order[rib]-1, dphin);
                // THIS TRANSFORMATION WILL BE PERFORMED INSIDE THE LOOP BELOW.
                // TransformDerivativeFromRibToTriang(rib,ord2,dphin);
                for (int i = 0; i < ord2; i++) // loop sobre as "ord2" funções de Chebyshev.
                {
                    phi(shape,0) = phiblend(rib+3,0)*phin(i,0);

                    dphi(0,shape) = dphiblend(0,rib+3)*phin(i,0) + phiblend(rib+3,0)*dphin(0,i)*gRibTrans2dT1d[rib][0]; // dx.

                    dphi(1,shape) = dphiblend(1,rib+3)*phin(i,0) + phiblend(rib+3,0)*dphin(0,i)*gRibTrans2dT1d[rib][1]; // dy.

                    dphi(2, shape) = dphiblend(2,rib+3)*phin(i,0) + dphiblend(0,rib+3)*dphin(0,i)*gRibTrans2dT1d[rib][0] + 
                                dphiblend(0,rib+3)*dphin(0,i)*gRibTrans2dT1d[rib][0]
                                + phiblend(rib+3,0)*dphin(1,i)*gRibTrans2dT1d[rib][0]*gRibTrans2dT1d[rib][0]; // dxdx

                    dphi(3, shape) = dphiblend(3,rib+3)*phin(i,0) + dphiblend(0,rib+3)*dphin(0,i)*gRibTrans2dT1d[rib][1] + 
                                dphiblend(1,rib+3)*dphin(0,i)*gRibTrans2dT1d[rib][0]
                                + phiblend(rib+3,0)*dphin(1,i)*gRibTrans2dT1d[rib][0]*gRibTrans2dT1d[rib][1]; // dydx

                    dphi(4, shape) = dphiblend(4,rib+3)*phin(i,0) + dphiblend(1,rib+3)*dphin(0,i)*gRibTrans2dT1d[rib][0] + 
                                dphiblend(0,rib+3)*dphin(0,i)*gRibTrans2dT1d[rib][1]
                                + phiblend(rib+3,0)*dphin(1,i)*gRibTrans2dT1d[rib][1]*gRibTrans2dT1d[rib][0]; // dxdy

                    dphi(5, shape) = dphiblend(5,rib+3)*phin(i,0) + dphiblend(1,rib+3)*dphin(0,i)*gRibTrans2dT1d[rib][1] + 
                                dphiblend(1,rib+3)*dphin(0,i)*gRibTrans2dT1d[rib][1]
                                + phiblend(rib+3,0)*dphin(1,i)*gRibTrans2dT1d[rib][1]*gRibTrans2dT1d[rib][1]; // dydy

                    shape++;
                }
            }

            if (order[3] < 3) return;//ordem na face
            int ord =  order[3]-2;//num de shapes da face
            int nsh = (ord*(ord+1))/2;

            REAL store1[20],store2[20],store3[20],store4[20];
            TPZFMatrix phi0(nsh,1,store1,20), phi1(nsh,1,store2,20), dphi0(2,nsh,store3,20),dphi1(2,nsh,store4,20);

            // esta função será retirada e reimplementada abaixo.
            // ShapeInternal(pt,order[3]-2,phin,dphin,GetTransformId2dT(id));

            int transid = GetTransformId2dT(id);
            TPZManVector<REAL> outPt(2);

            //TransformPoint2dT( transid, pt, outPt);

            TPZShapeLinear::fOrthogonal(2.*pt[0]-1., nsh, phi0, dphi0);
            TPZShapeLinear::fOrthogonal(2.*pt[1]-1., nsh, phi1, dphi1);

            TPZFMatrix temp_dphi( 6, nsh );
            TPZFMatrix temp_phi(nsh, 1);
            int index = 0;
            int i;
            for (int iplusj = 0; iplusj < ord; iplusj++)
            {
                for (int j=0;j<=iplusj;j++)
                {
                    i = iplusj-j;
                    temp_phi(index, 0) = phi0(i,0) * phi1(j,0);

                    temp_dphi(0,index) = 2. * dphi0(0,i)*phi1(j,0); // derivada do temp_phi(index,0) em relação à qsi.
                    temp_dphi(1,index) = 2. * phi0(i,0)*dphi1(0,j); // derivada do temp_phi(index,0) em relação à nu.

                    temp_dphi(2,index) = 4. * dphi0(1,i)*phi1(j,0); // derivada do temp_dphi(0,index) em relação à qsi.
                    temp_dphi(3,index) = 4. * dphi0(0,i)*dphi1(0,j); // derivada do temp_dphi(0,index) em relação à nu.

                    temp_dphi(4,index) = 4. * dphi0(0,i)*dphi1(0,j); // derivada do temp_dphi(1,index) em relação à qsi.
                    temp_dphi(5,index) = 4. * phi0(i,0)*dphi1(1,j); // derivada do temp_dphi(1,index) em relação à nu.

                    index++;
                }
            }

            for(index = 0; index < nsh; index++) //number of internal shape equal maximal order
            {
                phi(shape,0) = phiblend(6,0)*temp_phi(index,0);

                dphi(0,shape) = dphiblend(0,6)*temp_phi(index,0) + phiblend(6,0)*temp_dphi(0,index);    // dx

                dphi(1,shape) = dphiblend(1,6)*temp_phi(index,0) + phiblend(6,0)*temp_dphi(1,index);    // dy

                dphi(2,shape) = dphiblend(2,6)*temp_phi(index,0) + dphiblend(0,6)*temp_dphi(0,index)
                                + dphiblend(0,6)*temp_dphi(0,index)  + phiblend(6,0)*temp_dphi(2,index);    // dxdx

                dphi(3,shape) = dphiblend(3,6)*temp_phi(index,0) + dphiblend(0,6)*temp_dphi(1,index)
                                + dphiblend(1,6)*temp_dphi(0,index)  + phiblend(6,0)*temp_dphi(3,index);    // dydx

                dphi(4,shape) = dphiblend(4,6)*temp_phi(index,0) + dphiblend(1,6)*temp_dphi(0,index)
                                + dphiblend(0,6)*temp_dphi(1,index)  + phiblend(6,0)*temp_dphi(4,index);    // dxdy

                dphi(5,shape) = dphiblend(5,6)*temp_phi(index,0) + dphiblend(1,6)*temp_dphi(1,index)
                                + dphiblend(1,6)*temp_dphi(1,index)  + phiblend(6,0)*temp_dphi(5,index);    // dydy

                shape++;
            }

/*
            std::cout << " phi and dphi at integration node: " << pt[0] << ", " << pt[1] << "\n";
            phi.Print();
            dphi.Print();
*/
        }
        else
        {
            TPZFNMatrix<100> phiblend(NSides,1),dphiblend(Dimension,NSides);
            for(is=0; is<NCornerNodes; is++)
            {
                phiblend(is,0) = phi(is,0);
                for(d=0; d<Dimension; d++)
                {
                    dphiblend(d,is) = dphi(d,is);
                }
            }
            ShapeGenerating(pt,phiblend,dphiblend);
            
            REAL out;
            int shape = 3;
            for (int rib = 0; rib < 3; rib++)
            {
                ProjectPoint2dTriangToRib(rib,pt,out);
                TPZManVector<REAL,1> outvec(1,out);
                TPZVec<int> ids(2);
                ids[0] = id[rib%3];
                ids[1] = id[(rib+1)%3];
                REAL store1[20],store2[40];
                int ord2 = order[rib]-1;//numero de shapes por lado rib
                TPZFMatrix phin(ord2,1,store1,20),dphin(2,ord2,store2,40);
                TPZShapeLinear *shplin=0;
                shplin->ShapeInternal(outvec,order[rib],phin,dphin,shplin->GetTransformId1d(ids));
                TransformDerivativeFromRibToTriang(rib,ord2,dphin);
                for (int i = 0; i < ord2; i++) {
                    phi(shape,0) = phiblend(rib+3,0)*phin(i,0);
                    for(int xj=0;xj<2;xj++) {
                        dphi(xj,shape) = dphiblend(xj,rib+3)* phin( i, 0)+ phiblend(rib+3, 0 )* dphin(xj,i);
                    }
                    shape++;
                }
            }
            if (order[3] < 3) return;//ordem na face
            REAL store1[20],store2[40];
            int ord =  order[3]-2;//num de shapes da face
            int nsh = (ord*(ord+1))/2;
            TPZFMatrix phin(nsh,1,store1,20),dphin(2,nsh,store2,40);
            ShapeInternal(pt,order[3]-2,phin,dphin,GetTransformId2dT(id));
            for(int i=0;i<nsh;i++)	{//number of internal shape equal maximal order
                phi(shape,0) = phiblend(6,0)*phin(i,0);
                for(int d=0;d<2;d++) {
                    dphi(d,shape) = dphiblend(d,6)* phin(i,0) + phiblend(6,0)*dphin(d,i);
                }
                shape++;
            }

        }
	}

	void TPZShapeTriang::SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
								   TPZFMatrix &phi,TPZFMatrix &dphi) {
		if(side<0 || side>6) PZError << "TPZShapeTriang::SideShape. Bad paramenter side.\n";
		else if(side==6) Shape(pt,id,order,phi,dphi);
		else if(side<3) {
			TPZShapePoint::Shape(pt,id,order,phi,dphi);
		} else {
			TPZShapeLinear::Shape(pt,id,order,phi,dphi);
		}
	}

	void TPZShapeTriang::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
									   TPZFMatrix &dphi,int triangle_transformation_index)
    {
		if(order < 0) return;
		int ord1 = order;
		int numshape = (ord1*(ord1+1))/2;
		TPZManVector<REAL> out(2);
		TransformPoint2dT(triangle_transformation_index,x,out);

		if (phi.Rows() < numshape || dphi.Cols() < numshape)
        {
			PZError << "\nTPZCompEl::Shape2dTriangleInternal phi or dphi resized\n";
			phi.Resize(numshape,1);
			dphi.Resize(dphi.Rows(),numshape);
		}

		REAL store1[20],store2[20],store3[20],store4[20];
		TPZFMatrix phi0(numshape,1,store1,20),phi1(numshape,1,store2,20),dphi0(2,numshape,store3,20),dphi1(2,numshape,store4,20);

		TPZShapeLinear::fOrthogonal(2.*out[0]-1.,numshape,phi0,dphi0);
		TPZShapeLinear::fOrthogonal(2.*out[1]-1.,numshape,phi1,dphi1);
		int index = 0;
		int i;
		for (int iplusj=0;iplusj<ord1;iplusj++)
        {
			for (int j=0;j<=iplusj;j++)
            {
				i = iplusj-j;
				phi(index,0) = phi0(i,0)*phi1(j,0);
				dphi(0,index) = 2.*dphi0(0,i)*phi1(j,0);
				dphi(1,index) = 2.*phi0(i,0)*dphi1(0,j);
				index++;
			}
		}

		TransformDerivative2dT(triangle_transformation_index,numshape,dphi);
	}
	
	void TPZShapeTriang::ProjectPoint2dTriangToRib(int rib, TPZVec<REAL> &in, REAL &out) {
		
		out = gRibTrans2dT1d[rib][0]*in[0]+gRibTrans2dT1d[rib][1]*in[1]+gVet1dT[rib];
	}
	
	void TPZShapeTriang::TransformDerivativeFromRibToTriang(int rib,int num,TPZFMatrix &dphi) {
		
		for (int j = 0;j<num;j++) {
			
			dphi(1,j) = gRibTrans2dT1d[rib][1]*dphi(0,j);
			dphi(0,j) = gRibTrans2dT1d[rib][0]*dphi(0,j);
		}
	}
	
	int TPZShapeTriang::GetTransformId2dT(TPZVec<int> &id) {
		
		int id0,id1,minid;
		id0 = (id[0] < id[1]) ? 0 : 1;
		minid = (id[2] < id[id0]) ? 2 : id0;
		id0 = (minid+1)%3;
		id1 = (minid+2)%3;
		
		if (id[id0] < id[id1]) {//antihorario
			
			if (minid == 0) return 0;
			if (minid == 1) return 2;
			if (minid == 2) return 4;
			
		} else {//horario
			
			if (minid == 0) return 1;
			if (minid == 1) return 3;
			if (minid == 2) return 5;
		}
		return 0;
	}
	
	//transf. o ponto dentro da face triangular
	void TPZShapeTriang::TransformPoint2dT(int transid, TPZVec<REAL> &in, TPZVec<REAL> &out) {
		
		out[0] = gTrans2dT[transid][0][0]*in[0]+gTrans2dT[transid][0][1]*in[1]+gVet2dT[transid][0];
		out[1] = gTrans2dT[transid][1][0]*in[0]+gTrans2dT[transid][1][1]*in[1]+gVet2dT[transid][1];
	}
	
	void TPZShapeTriang::TransformDerivative2dT(int transid, int num, TPZFMatrix &in) {
		
		int i;
		for(i=0;i<num;i++) { //ds/dcsi
			REAL aux[2];
			aux[0] = in(0,i);
			aux[1] = in(1,i);
			in(0,i) = gTrans2dT[transid][0][0]*aux[0]+gTrans2dT[transid][1][0]*aux[1];
			in(1,i) = gTrans2dT[transid][0][1]*aux[0]+gTrans2dT[transid][1][1]*aux[1];
		}
	}
	
	int TPZShapeTriang::NConnectShapeF(int side, int order) {
		switch(side) {
			case 0:
			case 1:
			case 2:
				return 1;
			case 3:
			case 4:
			case 5:
				return order-1;
			case 6:
				return (order-2) < 0 ? 0 : ((order-2)*(order-1))/2;
			default:
				PZError << "TPZShapeTriang::NConnectShapeF, bad parameter iconnect " << side << endl;
				return 0;
		}
	}
	
	int TPZShapeTriang::NShapeF(TPZVec<int> &order) {
		int in,res=NCornerNodes;
		for(in=NCornerNodes;in<NSides;in++) res += NConnectShapeF(in,order[in-NCornerNodes]);
		return res;
	}

};
