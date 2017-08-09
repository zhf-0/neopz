/*********************************************************************
 *
 *  Gmsh circle
 *
 *********************************************************************/
Mesh.Algorithm = 3;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
General.ExpertMode = 1;
//h_cell_overall = 1.0;
//h_cell_azimuthal = 20;
r=0.01;
lc=r;

Point(1) = {0,0,0,lc};

Point(2) = {r,0,0,lc};
Point(3) = {0,r,0,lc};
Point(4) = {-r,0,0,lc};
Point(5) = {0,-r,0,lc};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Field[1] = Attractor;
Field[1].NodesList = {1};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc;
Field[2].LcMax = lc/3;
Field[2].DistMin = 0.25*r;
Field[2].DistMax = 0.9*r;
Field[2].Sigmoid = 0;

Background Field = 2;

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

//Transfinite Line {l1,l2,l3,l4} = h_cell_azimuthal;

// Tagging volume and boundary domains
Physical Surface("Omega") = {6};
Physical Line("Gamma") = {1,2,3,4};

Color Blue{ Line{ 1,2,3,4 }; }
Color Red{ Surface{ 6 }; }
