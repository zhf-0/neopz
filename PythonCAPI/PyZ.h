#ifndef PYZH
#define PYZH

#include <python2.7/Python.h>

//#include <boost/shared_ptr.hpp>
//#include "config.h"
//
//#include "pyanalysis.h"
//#include "pygraphmesh.h"
//#include "pymesh.h"
//#include "pyqa.h"

// Singleton class used to define a Python interface to ISET

template <class T, Subscript n0> class ConcreteRigidArray1d;
template <class T>               class ConcreteFormedArray1d;

class Analysis;
class ExtractSIF;
class DataFileO;

namespace {

#ifdef __ICC
// JAP: apparently we need 'static inline' for Intel...
PyMODINIT_FUNC static inline
#else
PyMODINIT_FUNC inline
#endif
registerPythonModule(const std::string& moduleName,
                     PyMethodDef* pythonMethods) {
  Py_InitModule(moduleName.c_str(), pythonMethods);
}

}

class PyZ {

 protected:
  // Clients must use method instance (see below)
  PyZ() {}

 public:
  // The only way to get the singleton object of this class
  static PyZ& instance();
  ~PyZ() {}

  // Python interface registration: Create interface for ISET.
  void registerPythonModules() const;

  // ***** Python/C++ data interaction utility methods *****

  // utility to get a global or local analysis
  static bool getAnalysis(int myLocalID,
                          boost::shared_ptr<Analysis>& shp_anal);

  // utility to parse a set of coordinates given a Python list
  static bool parseCoordinateList(PyObject* listObj,
                                  ConcreteRigidArray1d<realType,3>& coord);

  // utility to parse a set of p-orders given a Python list
  static bool parsePolyOrderList(PyObject* listObj,
                                 ConcreteRigidArray1d<unsigned char,3>& order,
                                 int& numDir);

  // utility to parse a set of strings given a Python list
  static bool parseStringList(PyObject* listObj,
                              ConcreteFormedArray1d<std::string>& strings);

 private:

  PyMesh myMesh;
  PyAnalysis myAnalysis;
  PyGraphMesh myGraphMesh;
  PyQA myQA;

  // ***** Functions that implement Python interface *****

  // Implements "create_local_problem"
  static PyObject* cmdCreateLocalProblem(PyObject* self, PyObject* args,
                                         PyObject* keys);

  // Implements "set_parameters"
  static PyObject* cmdSetParameters(PyObject* self, PyObject* args,
                                    PyObject* keys);

  static PyMethodDef pythonMethods[3];

  PyZ(PyZ const&); // no copy constructor
  void operator=(PyZ const&); // no copy assign

};

#endif
