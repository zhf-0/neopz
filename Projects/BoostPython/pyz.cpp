#include "PZPythonFunc.h"
#include <boost/python.hpp>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzanalysis.h"
#include "TPZMatModelProblem.h"


BOOST_PYTHON_MODULE(pyz)
{
  using namespace boost::python;
  
  def( "readFile", readFile );
  
  class_<TPZGeoMesh>("TPZGeoMesh")
    .def("NElements", &TPZGeoMesh::NElements)
    .def("NNodes", &TPZGeoMesh::NNodes)
  ;
  
  class_<TPZCompMesh>("TPZCompMesh")
    .def("NEquations", &TPZCompMesh::NEquations)
    .def("NMaterials", &TPZCompMesh::NMaterials)
    .def("printMat", &TPZCompMesh::printMat)
  ;
  
  class_<TPZAnalysis>("TPZAnalysis",init<TPZCompMesh*,bool>())
    .def("Assemble", &TPZAnalysis::Assemble)
    .def("Run", &TPZAnalysis::Run)
    .def("Print", &TPZAnalysis::Print)
    .def("GetTime", &TPZAnalysis::GetTime)
  
  ;
  
  def("CreateGMesh", &CreateGMesh);
  
  def("CMesh", &CMesh);
  
  
//  def( "postProcess", postProcess );
 
  def("returnAnInt", returnAnInt );
  
}
