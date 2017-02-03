// ---- Cylinder wellbore Region Gmsh scritp ----
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured tri region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

// Settings
ExpertMode = 1;


Include "CadReservoir.geo";
Include "BoxReservoir.geo";
Include "BoxSideBurden.geo";
Include "BuildOmegas.geo";
Include "drill_well.geo";
Include "PhysicalEntities.geo";

well_index = 0;
well_lids = {};

well_p_bores = {};
well_p_regions = {};
well_p_v_regions = {};

well_i_bores = {};
well_i_regions = {};
well_i_v_regions = {};


geomechanicQ = 1;
dimension = 3;
nolinearQ = 0;
CADReservoirQ = 0;

xzQ = 1;
hexahedronsWQ = 0;
hexahedronsRQ = 0;
hexahedronsSBQ = 0;

If (nolinearQ == 1)
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
EndIf

If (hexahedronsWQ == 1 || hexahedronsRQ == 1 || hexahedronsSBQ == 1)
Mesh.Algorithm3D = 6 ;
Else
Mesh.Algorithm3D = 1 ;
EndIf


// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 1;
cl2 = 0.1;
cl3 = 10.0;
cl4 = 150.0;
cl5 = 1000.0;

////////////////////////////////////////////////////////////////////////////
// reservoir region geometry
////////////////////////////////////////////////////////////////////////////

// reservoir box dimensions
x_length = 1000.0;
y_length = 1000.0;
z_length = 200.0;

////////////////////////////////////////////////////////////////////////////
// side-burden region geometry
////////////////////////////////////////////////////////////////////////////

// side-burden box dimensions
sb_x_length = 10000.0;
sb_y_length = 10000.0;
sb_z_length = 4000.0;

////////////////////////////////////////////////////////////////////////////
// reservoir rock
////////////////////////////////////////////////////////////////////////////
If(CADReservoirQ == 1)
Call ReservoirCAD;
Else
Call ReservoirBox;
EndIf

////////////////////////////////////////////////////////////////////////////
// side-burden rock
////////////////////////////////////////////////////////////////////////////
If (geomechanicQ == 1)
Call SideBurdenBox;
EndIf


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Drilling wells 
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

rw = 0.1;
wl = 50.0;

wbr = 20.0;
ela = 50.0;
rw_cell= 1.0;
wr_cell= 20.0;

// Orientation and length
alfa = Pi/2.0;
beta = 0.0;

////////////////////////////////////////////////////////////////////////////
// Drill producer 1 
////////////////////////////////////////////////////////////////////////////

// well location
wcx = 50.0;
wcy = 100.0;
wcz = 0.0;
IsInjectorQ = 0;
Call DrillWell;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Computational domain, boundaries:: Tagging physical entities
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

Call DefineOmegas;
Call DrawBoundaries;



