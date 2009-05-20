/* This header contains small hacks to add routines lacking in NTL.
 *
 * NTLhacks.hpp
 *
 *  Created on: May 19, 2009
 *      Author: defeo
 */

#ifndef NTLHACKS_HPP_
#define NTLHACKS_HPP_

#include <NTL/GF2X.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/GF2.h>
#include <NTL/lzz_p.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/tools.h>

NTL_OPEN_NNS
// this one is obsolete
void RightShiftAdd(GF2X& c, const GF2X& a, long n);
NTL_CLOSE_NNS

namespace NTL {
	/* Transforms a polynomial over GF2 in a monic polynomial.
	 * Of course, there is nothing to do.
	 */
	inline void MakeMonic(GF2X& x) {}
	/* Compute
	 * 		x = a^e (e may be negative)
	 * The routine lacks in NTL since in these cases
	 * it is stupid to pass a ZZ as exponent.
	 *
	 * Warning : we use the fact that p is prime !
	 */
	void power(zz_p& x, const zz_p& a, const ZZ& e);
	void power(GF2& x, const GF2& a, const ZZ& e);

	// this one is obsolete
	// assumes input does not alias output
	void RightShiftAdd(ZZ_pX& U, const ZZ_pX& V, long n);

	// this one is obsolete
	// assumes input does not alias output
	void RightShiftAdd(zz_pX& U, const zz_pX& V, long n);
}

#endif /* NTLHACKS_HPP_ */
