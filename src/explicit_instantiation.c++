/*
 * explicit_instantiation.c++
 *
 *  Created on: May 19, 2009
 *      Author: defeo
 */

#include "faast.hpp"

#include "Types.c++"
#include "Couveignes2000.c++"
#include "FE-Liftup-Pushdown.c++"
#include "FE-Trace-Frob.c++"
#include "Field.c++"
#include "FieldAlgorithms.c++"
#include "FieldElement.c++"
#include "FieldPolynomial.c++"
#include "FieldPrecomputations.c++"
#include "Minpols.c++"
#include "utilities.c++"
#include "NTLhacks.c++"

using namespace FAAST;

/* Classes */
template class Field<zz_p_Algebra>;
template class FieldElement<zz_p_Algebra>;
template class FieldPolynomial<zz_p_Algebra>;

template class Field<ZZ_p_Algebra>;
template class FieldElement<ZZ_p_Algebra>;
template class FieldPolynomial<ZZ_p_Algebra>;

template class Field<GF2_Algebra>;
template class FieldElement<GF2_Algebra>;
template class FieldPolynomial<GF2_Algebra>;
