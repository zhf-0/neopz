
#ifndef GRAFGRIDH
#define GRAFGRIDH


#include <iostream>

#include <string.h>
#include "pzcompel.h"
#include "pzadmchunk.h"
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzmaterial.h"

class TPZGraphNode;
class TPZCompMesh;
class TPZGraphEl;
class TPZFMatrix;
class TPZBlock;

enum TPZDrawStyle {EDXStyle,EMVStyle,EV3DStyle};

/// This class represents a graphical mesh used for post processing purposes
/**
The  graphical mesh represents a logically refined version of the computational mesh
This logical refinement means that the refined element object are not actually created
They only exist in the output file. 
*/
class TPZGraphMesh{
public:
  TPZGraphMesh(TPZCompMesh *cm, int dimension, TPZAutoPointer<TPZMaterial> mat);
  virtual ~TPZGraphMesh(void);
	
  TPZGraphNode &FindNode(int side);
  TPZGraphEl *FindElement(int sid);
  TPZAdmChunkVector<TPZGraphEl *> &ElementList();//MElementType type
  TPZAdmChunkVector<TPZGraphNode> &NodeMap();
  long NPoints();
  long NElements(MElementType type);
  int Res() {return fResolution;}
  void SetMaterial(TPZAutoPointer<TPZMaterial> mat) {fMaterial = mat;}
  virtual void SetCompMesh(TPZCompMesh *mesh, TPZAutoPointer<TPZMaterial> mat);
  std::ostream *Out();
  virtual void DrawNodes();
  virtual void DrawMesh(int numcases);
  virtual void DrawConnectivity(MElementType type);
  virtual void DrawSolution(TPZBlock &Sol);
  virtual void DrawSolution(char * var = 0);
  virtual void DrawSolution(int step, REAL time);
  TPZDrawStyle Style();
  void SetOutFile(std::ostream &out);
  void SetResolution(int res){ fResolution = res; SequenceNodes();}
	
  void Print(std::ostream &out);
  void SetNames(TPZVec<char *>&scalarnames, TPZVec<char *>&vecnames);
	
protected:
  TPZCompMesh *fCompMesh;
  TPZAutoPointer<TPZMaterial> fMaterial;
  int fDimension;
  TPZAdmChunkVector<TPZGraphEl *> fElementList;
  TPZAdmChunkVector<TPZGraphNode> fNodeMap;
  int fResolution;
  TPZDrawStyle fStyle;
  std::ostream *fOutFile;
  TPZVec<char *> fScalarNames, fVecNames;
  virtual void SequenceNodes();

  TPZCompEl *FindFirstInterpolatedElement(TPZCompMesh *mesh,int dimension);
	
public:
  virtual TPZAutoPointer<TPZMaterial> Material();
  virtual TPZCompMesh *Mesh() { return fCompMesh;}
};

inline TPZDrawStyle TPZGraphMesh::Style(){
  return fStyle;
}


#endif

