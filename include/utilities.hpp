#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <NTL/ZZ.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Types.hpp"
#include "Exceptions.hpp"

/* Miscellaneous routines */

namespace AS {
	/* Pollard Rho factorisation algorithm */
	void factor(long n, vector<pair<long,int> >& factors,
	const bool trial = true, const int multiplicity = 1);
	/* return the gratest power of p less than or equal to n */
	long NumPits(const long p, const long n);
	/* Store in res the composition Q(R).
	 * Q and R are two polynomials over GF(p).
	 * See also section 3.
	 */
	template <class T> void compose
	(typename T::GFpX& res, const typename T::GFpX& Q,
	const typename T::GFpX& R, const typename T::BigIny p);
	/* compute the n-th cyclotomic polynomial modulo p */
	template <class T> void cyclotomic
	(typename T::GFpX& res, const long n, const typename T::BigIny p);
}

#include "../src/utilities.c++"

#endif /*UTILITIES_H_*/
