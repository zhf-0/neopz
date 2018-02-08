/* ---------------------------------------------------------------------
 * The Parameter Handler class was adapted from deal.II library.
 * Therefore, the original copyright notice is preserved here.
 * ---------------------------------------------------------------------
 * Author: Francisco Orlandini, 2018
 */

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2005 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Luca Heltai, Wolfgang Bangerth, 2005
 */

#include "parameter_handler.h"
#include <list>
#include <iostream>
#include <fstream>

void
print_usage_message (const ParameterHandler &prm)
{
  static const char *message
    =
      "\n"
      "Example of reading parameters from a text file.\n"
      "\n"
      "Usage:\n"
      "    ./readParamProj file \n"
      "    ./readParamProj -p | --print <defaultfile> \n"
      "    ./readParamProj -h | --help\n"
      "\n"
      "Options:\n"
      "    -p --print Print default values to the file defaultfile\n"
      "    -h --help  Shows this usage message. \n"
      "\n"
      "The parameter file has the following format and allows the following\n"
      "values:\n"
      "\n";
  std::cout << message;
  prm.print_parameters (std::cout, ParameterHandler::Text);
}

void declare_parameters (ParameterHandler *prm)
{
  prm->declare_entry ("Output file", "",
                     Patterns::Anything(),
                     "The name of the output file to be generated");
  prm->declare_entry ("Dummy iterations", "42",
                     Patterns::Integer (1,1000),
                     "A dummy parameter asking for an integer");
  prm->enter_subsection ("Dummy subsection");
  {
    prm->declare_entry ("Dummy generate output", "true",
                       Patterns::Bool(),
                       "A dummy parameter that can be fed with either "
                       "'true' or 'false'");
    prm->declare_entry ("Dummy color of output", "red",
                       Patterns::Selection("red|black|blue"),
                       "A dummy parameter that shows how one can define a "
                       "parameter that can be assigned values from a finite "
                       "set of values");
  }
  prm->leave_subsection ();
}

void
parse_command_line (ParameterHandler &prm, const int     argc,
                    char *const *argv)
{
  if (argc < 2 || argc > 3)
  {
    print_usage_message (prm);
    exit (1);
  }
  //testing for -p or --print
  if(std::string(argv[1]) == std::string("-p")
    ||
    std::string(argv[1]) == std::string("--print")){
    if(argc != 3){
      exit(1);
    }
    std::ofstream outFile(argv[2]);
    prm.print_parameters (outFile, ParameterHandler::Text);
    outFile.close();
    return;
  }

  //testing for -h or --help
  if(std::string(argv[1]) == std::string("-h")
    ||
    std::string(argv[1]) == std::string("--help")){
    print_usage_message (prm);
    if(argc != 2){
      exit(1);
    }
    return;
  }

  const std::string parameter_file = argv[1];
  prm.parse_input (parameter_file);

  std::cout<<"Output format: "<<prm.get ("Output format")<<std::endl;
  std::cout<<"Output file: "<<prm.get ("Output file")<<std::endl;
  prm.enter_subsection ("Dummy subsection");
  {
    std::cout<<"Dummy generate output: ";
    std::cout<<prm.get_bool ("Dummy generate output")<<std::endl;
    std::cout<<"Dummy color of output: ";
    std::cout<<prm.get ("Dummy color of output")<<std::endl;
  }
  prm.leave_subsection ();
}  

int main (int argc, char **argv)
{
  ParameterHandler prm;
  try
    {
      declare_parameters (&prm);
      parse_command_line (prm, argc, argv);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
  return 0;
}