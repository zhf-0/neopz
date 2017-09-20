/*********************************************************************
 *
 *  Gmsh circle
 *	Define a circle with embeded fractures
 *
 *********************************************************************/

Macro WellboreDomain

r_inner=0.1;
r_outer=10.0;
lc=0.2*r_outer;
n_wellbore = 10;
n_radial = 4*n_wellbore;
ratio = 1.1;

p0 = newp; Point(p0) = {0,0,0,lc};

p1 = newp; Point(p1) = {r_inner,0,0,lc};
p2 = newp; Point(p2) = {r_outer,0,0,lc};
p3 = newp; Point(p3) = {0,r_outer,0,lc};
p4 = newp; Point(p4) = {0,r_inner,0,lc};

l0 = newl; Line(l0) = {p1,p2};
l1 = newl; Circle(l1) = {p2,p0,p3};
l2 = newl; Line(l2) = {p3,p4};
l3 = newl; Circle(l3) = {p4,p0,p1};

ll0 = newll; Line Loop(ll0) = {l0,l1,l2,l3};
s0 = news; Plane Surface(s0) = {ll0};


///////////////////////////////////////////////////////////
// Refinment towards wellbore
///////////////////////////////////////////////////////////

Transfinite Line {l3} = n_wellbore Using Progression 1;
Transfinite Line {l0,-l2} = n_radial Using Progression ratio;

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
