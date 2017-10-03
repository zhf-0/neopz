//+
SetFactory("OpenCASCADE");
angle = 2* Pi/9;
cangle = Pi/2-angle;
Rc = 15;
thick = 0.06;
Rd = Rc/Sin(angle);
R1 = Rd-thick/2;
R2 = Rd+thick/2;
Point(0) = {0, 0, 0, 1.0};
Point(1) = {0, R1, 0, 0.5};
Point(2) = {0, R2, 0, 0.5};
msh=0.05;
msh2 = 0.2;
Point(3) = {R1 * Cos(cangle), R1 * Sin(cangle), 0, msh};
Point(4) = {R2 * Cos(cangle), R2 * Sin(cangle), 0, msh};
Point(5) = {R1 * Cos(cangle) + 0.5, R2 * Sin(cangle) , 0, msh2};
Point(6) = {R1 * Cos(cangle), R2 * Sin(cangle) - 0.5, 0, msh2};
Point(7) = {R1 * Cos(cangle) + 0.5, R2 * Sin(cangle) - 0.5,0, msh2};

//+
Circle(1) = {3, 0, 1};
//+
Circle(2) = {4, 0, 2};
//+
Line(3) = {1, 2};
//+
Line(4) = {4, 3};
//+
Line(5) = {3, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 5};
//+
Line(8) = {5, 4};
//+
Line Loop(10) = {4, 5, 6, 7, 8};
//+
Plane Surface(11) = {10};

Line Loop(12) = {1,4,-2,3};

Plane Surface(13) = {12};//+
//+
Physical Line("ZERO") = {7, 5, 8, 2, 1};
//+
Physical Line("SYM") = {3};
//+
Physical Line("MOMENT") = {4};
//+
Physical Surface("ELASTICITY") = {11, 13};
//+
Physical Line("SUPPORT") = {6};
