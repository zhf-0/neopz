
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
    prm.declare_entry("L2 error","false", Patterns::Bool(),
                      "If set to true, the error (L2 norm) will be calculated"
                          " for the electric field");
    prm.declare_entry("Export L2 error","false", Patterns::Bool(),
                      "File name in which the L2 error will be saved "
                          "if L2 error is true. If not set, it will "
                          "be shown in std::cout only");
    prm.declare_entry("Export eigenvalues","false",Patterns::Bool(),
                      "If set to true, eigenvalues will"
                          "be exported to a text file.");
    prm.declare_entry("Number of threads","4",Patterns::Integer(0),
                      "Number of threads to use in NeoPZ assembly.");
    prm.declare_entry("Suffix","",Patterns::Anything(),
                      "Suffix to be added to exported files");
  }
  prm.leave_subsection ();
  //IN THE FOLLOWING OPTIONS, -1 =  PETSC_DECIDE and -2 = PETSC_DEFAULT
  prm.enter_subsection("SLEPc solver options");
  {
    prm.declare_entry("Problem type", "EPS_NHEP",
                      Patterns::Selection("EPS_HEP|EPS_GHEP|EPS_NHEP"
                                              "|EPS_GNHEP|EPS_PGNHEP|EPS_GHIEP"),
                      "Sets the type of the eigenvalue problem type.\n"
                          "G stands for Generalized, H for Hermitian,\n"
                          "P for positive semi-definite B matrix,\n"
                          "and I for generalized Hermitian-indefinite.");
    prm.declare_entry("Eigensolver","EPSKRYLOVSCHUR",Patterns::Selection("EPSPOWER|EPSSUBSPACE|EPSARNOLDI|EPSLANCZOS|"
                                                                             "EPSKRYLOVSCHUR|EPSGD|EPSJD|EPSRQCG|"
                                                                             "EPSLOBPCG|EPSCISS|EPSLAPACK|EPSARPACK|"
                                                                             "EPSBLZPACK|EPSTRLAN|EPSBLOPEX|"
                                                                             "EPSPRIMME|EPSFEAST"),
    "Sets the particular eigensolver to be used");
    prm.declare_entry("Krylov locking","false",Patterns::Bool(),"If Eigensolver is Krylov-Schur,"
        " whether to use the locking variant or not.");
    prm.declare_entry("Krylov restart","0.5",Patterns::Double(0.,1.),"Sets the restart parameter for the "
        "Krylov-Schur method, in particular the proportion of basis vectors that must be kept after restart.");
    prm.declare_entry("Convergence test", "EPS_CONV_REL",
                      Patterns::Selection("EPS_CONV_ABS|EPS_CONV_REL"
                                              "|EPS_CONV_NORM"),
                      "Specifies how to compute the error"
                      " estimate used in the convergence test.");
    prm.declare_entry("True residual","false",Patterns::Bool(),"Whether to calculate the true residual"
        " if using shift-and-invert as a spectral transformation.");
    prm.declare_entry("Which eigenvalues", "EPS_LARGEST_MAGNITUDE",
                      Patterns::Selection("EPS_LARGEST_MAGNITUDE|EPS_SMALLEST_MAGNITUDE"
                                              "|EPS_LARGEST_REAL|EPS_SMALLEST_REAL|"
                                              "EPS_LARGEST_IMAGINARY|EPS_SMALLEST_IMAGINARY|"
                                              "EPS_TARGET_MAGNITUDE|EPS_TARGET_REAL|"
                                              "EPS_TARGET_IMAGINARY|EPS_ALL"),
                      "Specifies which portion of the spectrum is to be sought.");
    prm.declare_entry("Target eigenvalue", "-1.", Patterns::Double(),
                      "Target eigenvalue to be sought if Which eigenvalues is equal"
                          "to EPS_TARGET_(MAGNITUDE|IMAGINARY|REAL).");
    prm.declare_entry("Eigensolver tolerance", "-2.", Patterns::Double(),
                      "Eigensolver convergence tolerance(Default value is PETSC_DEFAULT)");
    prm.declare_entry("Eigensolver maximum iterations", "-2.", Patterns::Integer(),
                      "Maximum number of iterations of the eigensolver.");
    prm.declare_entry("How many eigenvalues(nev)", "5", Patterns::Integer(),
                      "Number of eigenvalues that must satisfy the convergence test");
    prm.declare_entry("Solver subspace dimension(ncv)", "-1", Patterns::Integer(0),
                      "The maximum dimension of the subspace to be "
                          "used by the eigensolver(Default value is PETSC_DECIDE)");
    prm.declare_entry("Maximum projected dimension(mpd)", "-1", Patterns::Integer(0),
                      "The maximum projected dimension of the "
                          "eigensolver(Default value is PETSC_DECIDE)");
    prm.declare_entry("Maximum projected dimension(mpd)", "-1", Patterns::Integer(0),
                      "The maximum projected dimension of the "
                          "eigensolver(Default value is PETSC_DECIDE)");
    prm.declare_entry("Solver verbosity", "false", Patterns::Bool(),
                      "Whether the solver should print data during the simulation.");
    prm.declare_entry("Preconditioner","PCNONE",Patterns::Selection("PCNONE|PCJACOBI|PCSOR|PCLU|"
                                                                      "PCSHELL|PCBJACOBI|PCMG|PCEISENSTAT|PCILU|"
                                                                      "PCICC|PCASM|PCGASM|PCKSP|"
                                                                      "PCCOMPOSITE|PCREDUNDANT"),
                      "Sets the preconditioner to be used.");
    prm.declare_entry("Linear solver","KSPPREONLY",Patterns::Selection("KSPRICHARDSON|KSPCHEBYSHEV|KSPCG|KSPGROPPCG|"
                                                                           "KSPPIPECG|KSPPIPECGRR|KSPCGNE|KSPCGNASH|"
                                                                           "KSPCGSTCG|KSPCGGLTR|KSPFCG|KSPPIPEFCG|"
                                                                           "KSPGMRES|KSPPIPEFGMRES|KSPFGMRES|"
                                                                           "KSPLGMRES|KSPDGMRES|KSPPGMRES|KSPTCQMR|"
                                                                           "KSPBCGS|KSPIBCGS|KSPFBCGS|KSPFBCGSR|"
                                                                           "KSPBCGSL|KSPPIPEBCGS|KSPCGS|KSPTFQMR|"
                                                                           "KSPCR|KSPPIPECR|KSPLSQR|KSPPREONLY"),
                      "Sets the preconditioner to be used.");
    prm.declare_entry("Linear solver relative tolerance", "-2", Patterns::Double(0.),
                      "Sets the relative convergence tolerance, relative decrease in"
                          " the residual norm");
    prm.declare_entry("Linear solver absolute tolerance", "-2", Patterns::Double(0.),
                      "Sets the absolute convergence tolerance absolute size of the "
                          "residual norm");
    prm.declare_entry("Linear solver divergence tolerance", "-2", Patterns::Double(0.),
                      "Sets the divergence tolerance, amount residual norm can increase before "
                          "the convergence test concludes that the method is diverging");
    prm.declare_entry("Linear solver maximum iterations", "-2", Patterns::Integer(0),
                      "Maximum number of iterations to use");
    prm.declare_entry("Spectral transformation","STSINVERT",Patterns::Selection("STSHELL|STSHIFT|STSINVERT|"
                                                                                    "STCAYLEY|STPRECOND|STFILTER"),
                      "Sets the spectral transformation to be used.");

  }
  prm.leave_subsection ();
}

void SPZModalAnalysisData::ReadParameters() {
  DeclareParameters();
}
