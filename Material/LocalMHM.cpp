#include "LocalMHM.h"

float eps = 1.0;
float rhs = 1.0;

LocalMHM::LocalMHM( int i, int d ) : TPZDiscontinuousGalerkin( i ), dim(d), fEpsilon(dim, dim, 0.0), LambdaProblem(false)
{
    for(int diag = 0; diag < dim; diag++)
        fEpsilon(diag, diag) = eps;
}

void LocalMHM::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef)
{
	TPZFMatrix temp_vec( dim, 1 );
    int k, l;
    double temp_sum;

    int r = data.phi.Rows();
    TPZVec<REAL> res(1, 0.0);

    if(LambdaProblem)
        res[0] = 0;
    else
        res[0] = rhs;   // inserir função...
    
	for(int j = 0; j < r; j++)
	{
		temp_vec.Zero();

   		for(k = 0; k < dim; k++)
        {
   			for(l = 0; l < dim; l++)
         		temp_vec(k,0) += fEpsilon(k,l) * data.dphix(l,j);
        }

        for(int i = 0; i < r; i++)
        {
            ef(i,0) += weight * res[0] * data.phi(i,0);

            temp_sum = 0.0;
			for(k = 0; k < dim; k++)
				temp_sum += temp_vec(k,0) * data.dphix(k,i);

			ek(i,j) += weight * temp_sum;
        }
	}
}

void LocalMHM::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef)
{
    TPZMaterialData& dataDiff = datavec[0];

    TPZFMatrix  &phiu =  datavec[0].phi;    
    TPZFMatrix &phip = datavec[1].phi;

    int phi_u = phiu.Rows();
    int phi_p = phip.Rows();

    TPZFMatrix cur_ek(phi_u, phi_u, 0.0);
    TPZFMatrix cur_ef(phi_u, 1, 0.0);

    // =========
    // Diffusion
    // =========

    Contribute(dataDiff, weight, cur_ek, cur_ef);

    for(int i = 0; i < phi_u; i++)
    {
        ef(i,0) += cur_ef(i,0);
        for(int j = 0; j < phi_u; j++)
            ek(i,j) += cur_ek(i,j);
    }

    // ===================
    // Lagrange Multiplier
    // ===================

    for(int i = 0; i < phi_p; i++)
    {
        ef(phi_u+i, 0) += 0.0;

        for(int j = 0; j < phi_u; j++)
        {
            ek(j, phi_u+i) += weight * phiu(j,0) * phip(i,0);  // B
            ek(phi_u+i, j) += weight * phiu(j,0) * phip(i,0);  // C
        }

        for(int j = 0; j < phi_p; j++)
            ek(phi_u+i, phi_u+j) += 0.0;    // 0
    }
}

void LocalMHM::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc)
{
    if( bc.Type() == 1 ) // COMPULSORY: It indicates to the NeoPZ that this brick is a Neumann BC.
    {
        int r = data.phi.Rows();

        TPZVec<REAL> res;
        res.Resize( 1 );

        if(LambdaProblem)
            bcFunc(data, res);
        else
            res[0] = 0.0;

        // ============
        // A formulação dos problemas locais exibe um produto interno entre dois vetores normais: um relativo à malha e outro relativo ao elemento.
        // ============

        for (int i = 0; i < r; i++)
            ef(i,0) +=  (-1) * weight * res[0] * data.phi(i,0);
    }
}

void LocalMHM::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc)
{
    if( bc.Type() == 1 )
    {
        TPZMaterialData& data = datavec[0];
        int sizephi = data.phi.Rows();
        TPZFMatrix cur_ek(sizephi, sizephi, 0.0);
        TPZFMatrix cur_ef(sizephi, 1, 0.0);

        ContributeBC(data, weight, cur_ek, cur_ef, bc);

        for(int i = 0; i < sizephi; i++)
        {
            ef(i,0) += cur_ef(i,0);

            for(int j = 0; j < sizephi; j++)
                ek(i,j) += cur_ek(i,j);
        }
    }
}

void LocalMHM::bcFunc(TPZMaterialData& data, TPZVec<REAL>& res)
{
    // check if x belongs to the curFace.
    TPZVec<REAL> qsi(3, 0.0);
    //std::cout << "face[" << curFace << "]: " <<  data.x[0] << ", " << data.x[1] << ", " << data.x[2] << std::endl;

    if ( IsXOnTheFace(geoFaces[curFace], data.x) ) // if true: we evaluate the basis functions which relates to the curDof.
    {
        TPZInterpolationSpace* cur = dynamic_cast<TPZInterpolationSpace*>( compFaces[curFace] );

        int nshape = cur->NShapeF();
        TPZFMatrix phi(nshape, 1, 0.0);
        TPZFMatrix dphi(nshape, nshape, 0.0);

        qsi.Resize( geoFaces[curFace]->Dimension() );
        TransformX(data.x, qsi); // it transforms the data.x in to a node in geoFaces[curFace]'s reference element.

        //std::cout << "qsi: " << qsi[0];

        cur->Shape(qsi, phi, dphi);
        res[0] = phi(curDof, 0);
    }
    else // if false: x belongs to the other face.
        res[0] = 0.0;

    //std::cout << ". Dof[" << curDof << "]: " << res[0] << std::endl;
}

bool LocalMHM::IsXOnTheFace( TPZGeoEl* geo, TPZVec<REAL>& X )
{
    REAL tol = 10E-6;
    TPZGeoNode *np;
    int nnodes = geo->NNodes();
    TPZFMatrix coord(3, nnodes, 0.0);

    for(int i = 0; i < nnodes; i++)
    {
        np = geo->NodePtr(i);
        for(int j = 0; j < 3; j++)
            coord(j,i) = np->Coord(j);
    }

    int dim = geo->Dimension();
    REAL line = sqrt( pow(coord(0,1) - coord(0,0),2) + pow(coord(1,1) - coord(1,0),2) );
    REAL dis = 0;
    switch (dim)
    {
        case (1) :
            // distance between point and line
            dis = fabs( (coord(0,1) - coord(0,0))*(coord(1,0) - X[1]) - (coord(0,0) - X[0])*(coord(1,1) - coord(1,0)) ) / line;
            break;
        case (2) :
            // distance between point and plane.
            dis = 1.0;
            break;
        default:
            DebugStop();
            break;
    }

    if( fabs(dis) < tol )
        return true;
    else 
        return false;
}

void LocalMHM::TransformX( const TPZManVector<REAL,3>& dataX, TPZVec<REAL>& qsi )
{
    TPZGeoNode *np;
    int nnodes = geoFaces[curFace]->NNodes();
    TPZFMatrix coord(3, nnodes, 0.0);

    for(int i = 0; i < nnodes; i++)
    {
        np = geoFaces[curFace]->NodePtr(i);
        for(int j = 0; j < 3; j++)
            coord(j,i) = np->Coord(j);
    }

    int dim = geoFaces[curFace]->Dimension();
    qsi.Resize(dim);

    switch (dim) 
    {
        case 1:
        {
            REAL line = sqrt( pow(coord(0,1) - coord(0,0),2) + pow(coord(1,1) - coord(1,0),2) );
            TPZFMatrix temp(3, 1, 0.0);
            for(int i = 0; i < 3; i++)
                temp(i,0) = dataX[i] - coord(i,0);
            qsi[0] = Norm(temp) / line;
            break;
        }

        case 2:
            // Under construction...
            break;

        default:
            break;
    }
}

int LocalMHM::VariableIndex(const std::string &name)
{
    if(!strcmp("Solution1", name.c_str() ) )
        return 19;

    if(!strcmp("Solution2", name.c_str() ) )
        return 20;

	return TPZMaterial::VariableIndex(name);
}

int LocalMHM::NSolutionVariables(int var)
{
    if( var == 19) return 1;
    if( var == 20) return 1;
    
	return TPZMaterial::NSolutionVariables(var);
}

void LocalMHM::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout)
{
	Solout.Resize( this->NSolutionVariables(var) );
	TPZVec<REAL> SolU, SolP;
	SolU = datavec[0].sol[0];
	SolP = datavec[1].sol[0];

    // function (state variable 1)
	if(var == 19)
    {
		Solout[0] = SolU[0];
		return;
	}

    // function (state variable 2)
	if(var == 20)
    {
		Solout[0] = SolP[0];
		return;
	}
}
