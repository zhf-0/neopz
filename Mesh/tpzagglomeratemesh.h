//
// C++ Interface: tpzagglomeratemesh
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@corona>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZAGGLOMERATEMESH_H
#define TPZAGGLOMERATEMESH_H

#include <pzflowcmesh.h>

/**
This class contains both discontinuous, continuous and agglomerated elements
Its distinction from other meshes is that it points to a reference fine mesh

@author Philippe R. B. Devloo
*/
class TPZAgglomerateMesh : public TPZFlowCompMesh
{
public:
    TPZAgglomerateMesh() : TPZFlowCompMesh(0)
    {
      fFineMesh = 0;
    }
    
    /**
    An agglomeratemesh needs a fine mesh to relate to, because its elements may
    point to elements of the finemesh
    */
    TPZAgglomerateMesh(TPZCompMesh *finemesh) : TPZFlowCompMesh(finemesh->Reference()),
        fFineMesh(finemesh)
    {
    }

virtual ~TPZAgglomerateMesh()
{
}
    
  /**
  Return a pointer to the associated fine mesh
  */
  TPZCompMesh *FineMesh()
  {
    return fFineMesh;
  }



private:
  
   /**
   Reference Mesh for agglomerated elements
   */
    TPZCompMesh *fFineMesh;

};

#endif
