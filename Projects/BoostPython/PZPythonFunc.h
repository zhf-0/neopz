#ifndef PZPythonFunc_H
#define PZPythonFunc_H

#include <boost/python.hpp>
#include "pzgmesh.h"

char const* readFile();
boost::python::object CreateGMesh(long nel, REAL elsize);
boost::python::object CMesh(TPZGeoMesh *gmesh, int pOrder);
//char const* createAnalysis();
//
//char const* createExactSol();
//char const* parSetExactSolAsBF();
//char const* parSetExactSolAsBC();
//
//char const* createGraphMesh();
//char const* assemble();
//char const* solve();
//char const* postProcess(const int &resolution);

const int returnAnInt();

#endif // PZPythonFunc_H
