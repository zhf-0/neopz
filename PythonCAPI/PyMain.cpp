// ---- pyfunct.C ---- Sat Oct 15 12:25:17 MDT 2016
// (c) Copyright: C. Armando Duarte, 2002-2016. All rights reserved
// See README.CopyRight file for further details.
//

#define DBUG_ON 1

#include <iostream>

#include <python2.7/Python.h>

#include "PyZ.h"

// **********************************************************
// **********************************************************

int mainPython(int argc, char** argv) {


  // set the program name
  Py_SetProgramName( argv[0] );

  // initialize the Python interpreter and create the __main__ module
  Py_Initialize();

  // add modules
  PyISET::instance().registerPythonModules();

  // set job name based on the Python file name
  if ( argc > 1 ) {
    const std::string fileName( argv[1] );
    const std::string::size_type idx = fileName.find_last_of('.');

    GlobalISET::instance().jobName() = ( idx == std::string::npos ) ?
        fileName : fileName.substr(0,idx);
  } else {
    GlobalISET::instance().jobName() = "ISET Python Analysis";
  }

  int rc2 = Py_Main(argc, argv);

  Py_Finalize();

  return rc2;
}



