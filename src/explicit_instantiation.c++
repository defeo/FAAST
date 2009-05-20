/*
 * explicit_instantiation.c++
 *
 *  Created on: May 19, 2009
 *      Author: defeo
 */

#include "Artin-Schreier.hpp"

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

using namespace AS;

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

/* From FieldElement.hpp */
/*
template void pushDown(const FieldElement<zz_p_Algebra>& e, vector<FieldElement<zz_p_Algebra> >& v) throw(NoSubFieldException);
template void liftUp(const vector<FieldElement<zz_p_Algebra> >& v, FieldElement<zz_p_Algebra>& e) throw(NotInSameFieldException, NoOverFieldException);

template void pushDown(const FieldElement<ZZ_p_Algebra>& e, vector<FieldElement<ZZ_p_Algebra> >& v) throw(NoSubFieldException);
template void liftUp(const vector<FieldElement<ZZ_p_Algebra> >& v, FieldElement<ZZ_p_Algebra>& e) throw(NotInSameFieldException, NoOverFieldException);

template void pushDown(const FieldElement<GF2_Algebra>& e, vector<FieldElement<GF2_Algebra> >& v) throw(NoSubFieldException);
template void liftUp(const vector<FieldElement<GF2_Algebra> >& v, FieldElement<GF2_Algebra>& e) throw(NotInSameFieldException, NoOverFieldException);
*/

/* From utilities.hpp */
/*
template void expand(zz_pX& res, const zz_pX& P, const long n);
template void contract(zz_pX& res, const zz_pX& P, const long n);
template void compose(zz_pX& res, const zz_pX& Q, const zz_pX& R, const long& p);
template void cyclotomic(zz_pX& res, const long n, const long& p);

template void expand(ZZ_pX& res, const ZZ_pX& P, const long n);
template void contract(ZZ_pX& res, const ZZ_pX& P, const long n);
template void compose(ZZ_pX& res, const ZZ_pX& Q, const ZZ_pX& R, const ZZ& p);
template void cyclotomic(ZZ_pX& res, const long n, const ZZ& p);

template void expand(GF2X& res, const GF2X& P, const long n);
template void contract(GF2X& res, const GF2X& P, const long n);
template void compose(GF2X& res, const GF2X& Q, const GF2X& R, const int& p);
template void cyclotomic(GF2X& res, const long n, const int& p);
*/
