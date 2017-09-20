/*********************************************************************
 *
 *  Gmsh circle
 *	Define a macro that draw a rectangle with three embeded fractures
 *
 *********************************************************************/



Macro RockSampleDomain

lx= 0.5; // a/2
ly= 1.0; // b/2
n_cartesian = 20;
r=Sqrt(lx*lx+ly*ly);
lc=0.2*r;


p1 = newp; Point(p1) = {0,0,0,lc};
p2 = newp; Point(p2) = {lx,0,0,lc};
p3 = newp; Point(p3) = {lx,ly,0,lc};
p4 = newp; Point(p4) = {0,ly,0,lc};

l0 = newl; Line(l0) = {p1,p2};
l1 = newl; Line(l1) = {p2,p3};
l2 = newl; Line(l2) = {p3,p4};
l3 = newl; Line(l3) = {p4,p1};


ll0 = newll; Line Loop(ll0) = {l0,l1,l2,l3};
s0 = news; Plane Surface(s0) = {ll0};



///////////////////////////////////////////////////////////
// Refinment towards stress boundaries
///////////////////////////////////////////////////////////

Transfinite Line {l1,l2} = n_cartesian Using Progression 1;

///////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////
// Tagging surface and boundary domains
Physical Surface("Omega") = {s0};
Physical Line("Gamma_D_uy") = {l0};
Physical Line("Gamma_D_ux") = {l3};
Physical Line("Gamma_N_S_geo") = {l2};
Physical Line("Gamma_N_P_inner") = {l3};

///////////////////////////////////////////////////////////
// Coloring surface and boundary domains
Color Green{ Surface{ s0 }; }
Color Blue{ Line{ l0,l1,l2,l3 }; }


If (QuadrilateralMeshQ == 1)
Recombine Surface {s0};
EndIf

Return
