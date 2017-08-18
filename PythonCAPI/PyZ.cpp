// ---- pyiset.C ---- Fri Oct 14 21:02:43 MDT 2016
// (c) Copyright: C. Armando Duarte, 2002-2016. All rights reserved
// See README.CopyRight file for further details.
//

#define DBUG_ON 1

#include <limits>
#include <string>
#include <sstream>
#include <typeinfo>

#include "PyZ.h"


PyMethodDef PyISET::pythonMethods[] = {
  { "create_local_problem", (PyCFunction) cmdCreateLocalProblem, METH_VARARGS | METH_KEYWORDS,
    "Create a local problem." },
  { "set_parameters", (PyCFunction) cmdSetParameters, METH_VARARGS | METH_KEYWORDS,
    "Set custom parameters for the problem." },
  { nullptr, nullptr, 0, nullptr }
};

// *****************************************************
// *****************************************************

PyISET& PyISET::instance() {
  static PyISET theInstance;
  return theInstance;
}

// *****************************************************
// *****************************************************

void PyISET::registerPythonModules() const {
  ::registerPythonModule("pyiset", pythonMethods);

  // register mesh, analysis, graphmesh, and QA modules
  myMesh.registerPythonModule();
  myAnalysis.registerPythonModule();
  myGraphMesh.registerPythonModule();
  myQA.registerPythonModule();
}

// *****************************************************
// *****************************************************

bool PyISET::getAnalysis(int myLocalID,
                         boost::shared_ptr<Analysis>& shp_anal) {

  // Get global or local analysis object

  if ( myLocalID == SpecialBasis::invalidSpBas ) {
    if ( !GlobalISET::instance().analysis().get() ) {
      PyErr_SetString(PyExc_RuntimeError, "A global Analysis must be created");
      return false;
    }

    shp_anal = GlobalISET::instance().analysis();
    return true;
  }

  SpecialBasis* p_spBasis =
      GlobalISET::instance().compMesh()->specialBasisWithID(myLocalID);
  if ( !p_spBasis ) {
    PyErr_Format(PyExc_ValueError, "Invalid local problem ID: %i", myLocalID);
    return false;
  }

  LocalProblem* p_localProb = dynamic_cast<LocalProblem*>( p_spBasis );
  if ( !p_localProb ) {
    PyErr_Format(PyExc_ValueError, "SpecialBasis with ID %i is not a "
                 "LocalProblem", myLocalID);
    return false;
  }

  DBGPRNT2(cout << "\nPyISET::getLocalGlobalAnalysis LocalProblem:"
           << *p_localProb << endl);

  shp_anal = p_localProb->analysis();
  return true;
}

// *****************************************************
// *****************************************************

bool PyISET::parseCoordinateList(PyObject* listObj,
                                 ConcreteRigidArray1d<realType,3>& coord) {

  // Parse a set of coordinates in a Python list, e.g., if a bounding box or
  // similar is specified.

  // TBD JAP: should we use NumPy PyArrayObj instead of lists?

  const int listSize = PyList_Size(listObj);
  if ( listSize != 3 ) {
    PyErr_SetString(PyExc_ValueError, "Coordinate list size must be 3 "
                    "(x,y,z)");
    return false;
  }

  for (int idir=0; idir<3; ++idir) {
    PyObject* item = PyList_GetItem(listObj, idir);
    coord(idir) = PyFloat_AsDouble(item);
  }

  return true;
}

// *****************************************************
// *****************************************************

bool PyISET::parseStringList(PyObject* listObj,
                             ConcreteFormedArray1d<std::string>& strings) {

  // Parse a set of strings in a Python list

  const int listSize = PyList_Size(listObj);
  strings.reshape(listSize);
  for (int iname=0; iname<listSize; ++iname) {
    PyObject* item = PyList_GetItem(listObj, iname);
    strings(iname) = PyString_AsString(item);
  }

  return true;
}

// *****************************************************
// *****************************************************

bool PyISET::parsePolyOrderList(PyObject* listObj,
                                ConcreteRigidArray1d<unsigned char,3>& order,
                                int &numDir) {

  // Parse a list of polynomial orders (x,y,z) in a Python list

  // TBD JAP: should we use NumPy PyArrayObj instead of lists?

  numDir = PyList_Size(listObj);
  if ( numDir != 3 && numDir != 2 ) {
    PyErr_SetString(PyExc_ValueError, "List size must be 2 or 3 (x,y[,z]) "
                    "depending on dimensionality of the problem");
    return false;
  }

  for (int idir=0; idir<numDir; ++idir) {
    PyObject* item = PyList_GetItem(listObj, idir);
    order(idir) = static_cast<unsigned char>( _PyInt_AsInt(item) );
  }

  DBGPRNT(std::cout << "\nparsePolyOrderList: porder: ["
                    << int(order(0)) << ", "
                    << int(order(1)) << ", "
                    << int(order(2)) << "]\n");

  return true;
}

// *****************************************************
// *****************************************************

PyObject* PyISET::cmdCreateLocalProblem(PyObject* self, PyObject* args,
                                        PyObject* keys) {

  // Python:
  // -------
  // myid = create_local_problem(bcname, localid=0, bbox_min=None, bbox_max=None,
  //                             file=None, physics='solids', nlayers=1)
  //
  // Arguments:
  // ----------
  // bcname:   <string> Name of the local boundary condition to use
  // localid:  <int> ID of local problem to create (default=next available)
  // bbox_min: <list> Specifies lower limit for bounding box containing local
  //           problem nodes (default=None)
  // bbox_max: <list> Specifies upper limit for bounding box containing local
  //           problem nodes (default=None)
  // filename: <string> Name of file containing predefined local problem mesh
  //           (default=None)
  // physics:  <string> Physics type associated with the mesh in file 'filename'
  //           (default='solids')
  // nlayers:  <int> Number of layers outside of selected region to create the
  //           local problem (default=1)
  //
  // Return value:
  // -------------
  // myid:     <int> handle to/ID of the local problem created

  if ( GlobalISET::instance().geoMesh().get() == nullptr ||
       GlobalISET::instance().compMesh().get() == nullptr ) {
    PyErr_SetString(PyExc_RuntimeError,
                     "Global GeoMesh and CompMesh must be created");
    return NULL;
  }

  static char* kwlist[] = {
    (char*) "bcname",
    (char*) "localid",
    (char*) "bbox_min",
    (char*) "bbox_max",
    (char*) "filename",
    (char*) "physics",
    (char*) "nlayers",
    (char*) NULL
  };

  // default values
  int my_localid = SpecialBasis::invalidSpBas;
  int nlayers = 1;
  PyObject* bbox_min_list = nullptr;
  PyObject* bbox_max_list = nullptr;
  char* filename = (char*) NULL;
  char* bcname   = (char*) NULL;
  char* physics  = (char*) "solids";

  if ( !PyArg_ParseTupleAndKeywords(args, keys, "s|iOOss", kwlist,
                                    &bcname,
                                    &my_localid,
                                    &bbox_min_list,
                                    &bbox_max_list,
                                    &filename,
                                    &physics) ) {
    PyErr_SetString(PyExc_TypeError, "Invalid arguments to create_local_problem()");
    return NULL;
  }

  // // If we have a local id, get a local analysis. Otherwise, get the global
  // // analysis
  // boost::shared_ptr<Analysis> shp_anal;
  // if ( !PyISET::instance().getAnalysis(my_localid, shp_anal) )
  //   return NULL;

  // name for local boundary condition
  const std::string locBCName(bcname);

  // if an ID was not passed in, generate one automatically
  if ( my_localid == SpecialBasis::invalidSpBas ) {
    my_localid = SpecialBasis::createUniqueUserID();
  }
  // negative IDs are not allowed as user IDs (only internally generated IDs
  // are negative)
  else if ( my_localid < 0 ) {
    PyErr_SetString(PyExc_ValueError, "A valid user special basis ID must be "
                    "an integer > 0");
    return NULL;
  }
  // if an ID was specified, make sure it isn't already in use
  else if ( GlobalISET::instance().compMesh()->specialBasisWithID(my_localid) ) {
    int next_available_ID = SpecialBasis::createUniqueUserID();
    PyErr_Format(PyExc_ValueError, "Special basis ID %i is already in use. "
                 "Next available ID is %i", my_localid, next_available_ID);
    return NULL;
  }

  LocalProblem* p_localProb = nullptr;

  //
  // Create local problem from bounding box limits
  //
  if ( bbox_min_list != nullptr || bbox_max_list != nullptr ) {
    // require both upper and lower limits for bounding box
    if ( bbox_min_list == nullptr || bbox_max_list == nullptr ) {
      PyErr_SetString(PyExc_RuntimeError, "bbox_min and bbox_max are both "
                      "required arguments for creating local problem in "
                      "bounding box");
      return NULL;
    }

    // ensure there are no conflicting arguments
    if ( filename != (char*) NULL ) {
      PyErr_SetString(PyExc_RuntimeError, "Conflicting arguments to "
                      "create_local_problem()");
      return NULL;
    }

    ConcreteRigidArray1d<realType,3> my_bbox_min, my_bbox_max;
    my_bbox_min = std::numeric_limits<realType>::max();
    my_bbox_max = std::numeric_limits<realType>::max();

    if ( !PyISET::parseCoordinateList(bbox_min_list, my_bbox_min) ||
         !PyISET::parseCoordinateList(bbox_max_list, my_bbox_max) )
      return NULL;

    // FIXME (JAP): linear hard-coded for now!
    p_localProb = new LocalProblem(my_localid, GlobalISET::instance().geoMesh(),
                                   GlobalISET::instance().multiPhysDir(),
                                   my_bbox_min, my_bbox_max, nlayers, locBCName,
                                   true /*isLinearAnal*/);
  }

  //
  // Create local problem from file
  //
  else if ( filename != (char*) NULL ) {
    const std::string localFileName(filename);
    const std::string physName(physics);
    // FIXME (JAP): linear hard-coded for now!
    p_localProb = new LocalProblem(my_localid, GlobalISET::instance().geoMesh(),
                                   GlobalISET::instance().multiPhysDir(),
                                   localFileName, physName,
                                   nlayers, locBCName, true /*isLinearAnal*/);
  }

  else {
    PyErr_SetString(PyExc_RuntimeError, "Local problem creation requires either a "
                    "bounding box or a mesh file to be specified");
    return NULL;
  }

  if( p_localProb == nullptr ) {
    PyErr_SetString(PyExc_RuntimeError, "Invalid local problem");
    return NULL;
  }

  // set up the default local problem job name for extract file
  std::stringstream localJobName;
  localJobName << GlobalISET::instance().jobName() << "_local_" << my_localid;
  p_localProb->analysis()->setJobName( localJobName.str() );

  // store this SpecialBasis in CompMesh
  GlobalISET::instance().compMesh()->store(p_localProb);
  GlobalISET::instance().compMesh()->increaseNumLocProbs();

  // return the ID of the local problem created.
  return Py_BuildValue("i", my_localid);
}

// **********************************************************
// **********************************************************

PyObject* PyISET::cmdSetParameters(PyObject* self, PyObject* args,
                                   PyObject* keys) {

  // Python:
  // -------
  // set_parameters(params, localid=None, physics='solids')
  //
  // Arguments:
  // ----------
  // params:   <dict> A dictionary of string parameter names and values of the
  //           corresponding parameters
  // localid:  <int> ID of local problem to set parameters for (default=None)
  // physics:  <string> Physics type to set parameters for (default='solids')

  static char* kwlist[] = {
    (char*) "params",
    (char*) "localid",
    (char*) "physics",
    (char*) NULL
  };

  // default values
  PyObject* my_params = nullptr;
  int my_localid      = SpecialBasis::invalidSpBas;
  char* my_physics    = (char*) "solids";

  if ( !PyArg_ParseTupleAndKeywords(args, keys, "O|iz", kwlist,
                                    &my_params,
                                    &my_localid,
                                    &my_physics) ) {
    PyErr_SetString(PyExc_TypeError, "Invalid arguments to set_parameters()");
    return NULL;
  }

  const std::string physName(my_physics);

  // If we have a local id, get a local analysis. Otherwise, get the global
  // analysis
  boost::shared_ptr<Analysis> shp_anal;
  if ( !PyISET::instance().getAnalysis(my_localid, shp_anal) )
    return NULL;

  // TODO JAP: add in multiphysics support
  //
  // EPhysType phys_type = EPhysNone;
  // if ( !PyISET::instance().getPhysType(physics, phys_type) )
  //   return NULL;

  //
  // begin parsing the dictionary of parameters
  //
  const int numParams = PyDict_Size(my_params);
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  for (int i=0; i<numParams; ++i) {

    // get the next key in the dict
    PyDict_Next(my_params, &pos, &key, &value);
    const std::string parName = PyString_AsString(key);

    // get the corresponding next bool value in the dict
    //
    // (JAP) only supporting bool-type parameters for now!
    if ( !PyBool_Check(value) ) {
      PyErr_SetString(PyExc_ValueError, "set_parameters() currently only "
                      "supports boolean parameter values");
      return NULL;
    }

    // (JAP) is there a better solution than this cast from int to bool?
    const bool parValue = _PyInt_AsInt(value);

    if (parName == "exactSolAsBC") {
      shp_anal->setExactSolAsBC(parValue);
    } else if (parName == "exactSolAsPtBC") {
      shp_anal->setExactSolAsPtBC(parValue);
    } else if (parName == "exactSolAsBF") {
      shp_anal->setExactSolAsBF(parValue);
    } else if (parName == "ignoreBCs") {
      shp_anal->setIgnoreBCs(parValue);
    } else if (parName == "directMethodForPtDirichletBC") {
      shp_anal->compMesh()->useDirectMethodPtDirichletBC() = parValue;
    } else if (parName == "scaleAndStabAlgo") {
      shp_anal->setScaleAndStabB4Fact(parValue);
    } else if (parName == "scaleGlobalMatrixB4Fact") {
      shp_anal->setScaleAndStabB4Fact(parValue);
    } else if (parName == "SGFEMSpecialBasisEnrichments") {
      shp_anal->compMesh()->useSGFEMSpecialBasisEnrichments(parValue);
    } else if (parName == "SGFEMBranchFnEnrichments" ) {
      shp_anal->compMesh()->useSGFEMBranchFnEnrichments(parValue);
    } else if (parName == "SGFEMPolyEnrichments") {
      shp_anal->compMesh()->useSGFEMPolyEnrichments(parValue);
    } else {
      PyErr_Format(PyExc_ValueError, "Invalid parameter name '%s'",
                   parName.c_str());
      return NULL;
    }
  }

  return Py_None;
}

// **********************************************************
// **********************************************************
