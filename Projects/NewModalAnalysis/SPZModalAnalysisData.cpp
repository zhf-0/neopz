
#include "SPZModalAnalysisData.h"

SPZModalAnalysisData::SPZModalAnalysisData(ParameterHandler &paramHandler) :
  prm(paramHandler)
{

}

void SPZModalAnalysisData::DeclareParameters() {
  prm.enter_subsection ("Physical options");
  {
    prm.declare_entry("Mesh file", "",
                      Patterns::Anything(),
                      "Path to .msh gmsh generated mesh");
    prm.declare_entry("Cut-off analysis", "false",
                      Patterns::Bool(),
                      "Whether to perform cut-off analysis on the waveguide or"
                          "to calculate its modes for a given "
                          "operational frequency");
    prm.declare_entry("Operational frequency", "25e+9",
                      Patterns::Double(0.),
                      "Operational frequency of the waveguide (it will "
                          "be ignored if Cut-off analysis is true");
    prm.declare_entry("Number of materials", "1",
                      Patterns::Integer(1.),
                      "How many dielectric materials are present "
                          "in the waveguide");
    prm.declare_entry("Electric permittivity vector", "1.",
                      Patterns::List(Patterns::Double(0.),1),
                      "The REAL part of the electric permittivity of "
                          "the dielectric materials, separated by commas");
    prm.declare_entry("Is lossless(epsilon)","true", Patterns:: Bool(),
                      "Whether the dielectric material has complex "
                          "permittivity or not.");
    prm.declare_entry("Dielectric losses vector", " ",
                      Patterns::List(Patterns::Double(0.),1),
                      "The IMAGINARY part of the electric permittivity"
                          " of the dielectric materials, separated by "
                          "commas (it will be ignored"
                          " if Is lossless(epsilon) is false");
    prm.declare_entry("Magnetic permeability vector", "1.",
                      Patterns::List(Patterns::Double(0.),1),
                      "The REAL part of the magnetic permeability of the "
                          "dielectric materials, separated by commas");
    prm.declare_entry("Is lossless(mu)","true", Patterns:: Bool(),
                      "Whether the dielectric material has complex permeability or not.");
    prm.declare_entry("Dielectric losses vector", " ",
                      Patterns::List(Patterns::Double(0.),1),
                      "The IMAGINARY part of the magnetic permeability"
                          " of the dielectric materials, "
                          "separated by commas (it will be ignored"
                          " if Is lossless(mu) is false");
  }
  prm.leave_subsection ();

  prm.enter_subsection("NeoPZ options");
  {
    prm.declare_entry("Polynomial order","1", Patterns::Integer(1),
                      "Default polynomial order of the Pk space used to build"
                          "the Nédélec elements(The H1 elements"
                          " will be built accordingly");
    prm.declare_entry("Generate VTK","false", Patterns::Bool(),
                      "If set to true, a .vtk plot of the electric field will be generated"
                          "for the calculated modes");
    prm.declare_entry("VTK file name","electricFieldPlot.vtk", Patterns::Anything(),
                      "File name in which the .vtk plots will be saved "
                          "if Generate VTK is true");
    prm.declare_entry("L2 error","false", Patterns::Bool(),
                      "If set to true, the error (L2 norm) will be calculated"
                          " for the electric field");
    prm.declare_entry("L2 error file name","", Patterns::Anything(),
                      "File name in which the L2 error will be saved "
                          "if L2 error is true. If not set, it will "
                          "be shown in std::cout only");
    prm.declare_entry("Export eigenvalues","false",Patterns::Bool(),
                      "If set to true, eigenvalues will"
                          "be exported to a text file.");
    prm.declare_entry("Eigenvalues file","eigenvalues.txt",Patterns::Anything(),
    "File in which the eigenvalues will be saved.");
    prm.declare_entry("Number of threads","0",Patterns::Integer(0),
                      "Number of threads to use in NeoPZ assembly.");
  }
  prm.leave_subsection ();
}

void SPZModalAnalysisData::ReadParameters() {
  DeclareParameters();
}
