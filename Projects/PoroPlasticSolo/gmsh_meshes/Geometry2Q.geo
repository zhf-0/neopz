

IsquadQ = 1;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

lc = 0.9;
n_bc = 4;
nLayers    = 8;


  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {1, 0, 0, 1e+22};
  Point(3) = {1, 10, 0, 1e+22};
  Point(4) = {0, 10, 0, 1e+22};

// Fronteiras

  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};  

  Transfinite Line{1,3} = 1 ;
  Transfinite Line{2,-4} = nLayers Using Progression lc;



// Definição da superfície 

  Line Loop(93) = {1, 2, 3, 4};
 
  Plane Surface(110) = {93};

  If(IsquadQ)

  Recombine Surface {110};

  EndIf

 // Physical Volume("internal") = {1};
 // Extrude {0, 0, 10} {
 //  //Surface{110};
 //  Layers{1};
 //  Recombine;
 // }

  Physical Surface("Omega") = {110};
  Physical Line("bottom") = {1};
  Physical Line("top") = {2};
  Physical Line("left") = {3};
  Physical Line("right") = {4};
  
  //Physical Line("holes") = {holes[]};  
  //Physical Surface("interface") = {23};

  
  Coherence Mesh;




