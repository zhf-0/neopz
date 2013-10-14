#include "pztransientmat.h"

#include "pzpoisson3d.h"
template class TPZTransientMaterial< TPZMatPoisson3d >;

#include "pznonlinearpoisson3d.h"
template class TPZTransientMaterial< TPZNonLinearPoisson3d >;

#include "pzburger.h"
template class TPZTransientMaterial< TPZBurger >;

// teste..
// #include "Users/moniquetelo/Documents/Fred/LuaPZ/Bricks/TPZMatComposite.h"
// template class TPZTransientMaterial< TPZMatComposite >;

/** @brief Instantiations to TPZMatPoisson3d, TPZNonLinearPoisson3d and TPZBurger. */
void TestInstantiations(){
	TPZTransientMaterial< TPZMatPoisson3d > A(1,1,1.);
	TPZTransientMaterial< TPZNonLinearPoisson3d > B(1,1,1.);
	TPZTransientMaterial< TPZBurger > C(1,1,1.);
}