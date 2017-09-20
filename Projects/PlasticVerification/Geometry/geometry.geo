
/*********************************************************************
 *
 *  Gmsh control scritp
 *	Two dimensional examples with three embeded fractures
 *
 *********************************************************************/

Include "rock_sample.geo";
Include "wellbore_region.geo";

Mesh.Algorithm = 1;
Mesh.CharacteristicLengthExtendFromBoundary = 1;
General.ExpertMode = 1;

// Controls 
WellboreGeometryQ = 1;
QuadrilateralMeshQ = 1;

If (WellboreGeometryQ == 1)
Call WellboreDomain;
Else
Call RockSampleDomain;
EndIf

