#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <NTL/ZZ.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Types.hpp"
#include "Exceptions.hpp"

namespace AS {
	/**
	 * \cond DEV
	 * \defgroup Utilities Utility routines
	 *  Miscellaneous utility routines
	 *  @{
	 */

	/**
	 * \brief Pollard Rho factorisation algorithm.
	 *
	 * Uses Pollard Rho factorisation + NTL's primality test
	 * (Miller-Rabin test) to completely factor an integer.
	 *
	 * \param [in] n The integer to factor
	 * \param [in,out] factors An ordered vector to hold pairs of
	 * 		<factor,multiplicity> (the order is relative to the factors.
	 * 		If it is not empty, assuming it is ordered,
	 * 		factors of \a n are
	 * 		smoothly added to the vector without smashing the pre-extant
	 * 		<values,multiplicities>.
	 * \param [in] trial Do some trial division before Pollard Rho. Default
	 * 		is yes.
	 * \param [in] multiplicity All the multiplicity of the factors of \a n
	 * 		shall by scaled by this integer. Default is one.
	 */
	void factor(long n, vector<pair<long,int> >& factors,
	const bool trial = true, const int multiplicity = 1);
	/**
	 * \brief The least power of \a p greater than \a n.
	 *
	 * When \a p equals 2, this is equivalent to NTL's \c NumBits(n).
	 */
	long NumPits(const long p, const long n);
	/**
	 *  \brief Compute \a P(X<sup>n</sup>).
	 *
	 *  Compute the polynomial obtained by composing \a P with X<sup>n</sup>.
	 *  It implements the algorithm
	 *  referred as \c Evaluate in Section 4 of [\ref ISSAC "DFS '09"].
	 *
	 *  \param [out] res A polynomial to hold the result.
	 *  \param [in] P A polynomial.
	 *  \param [in] n A strictly positive integer.
	 */
	template <class T> void expand(typename T::GFpX& res,
	const typename T::GFpX& P, const long n);
	/**
	 *  \brief Transposition of expand().
	 *
	 *  From polynomial \a P, deduce a polynomial whose k-th coefficient
	 *  is the nk-th coefficient of \a P. It implements the algorithm
	 *  referred as \c Evaluate* in Section 4 of [\ref ISSAC "DFS '09"].
	 *
	 *  \param [out] res A polynomial to hold the result.
	 *  \param [in] P A polynomial.
	 *  \param [in] n A strictly positive integer.
	 */
	template <class T> void contract(typename T::GFpX& res,
	const typename T::GFpX& P, const long n);
	/**
	 * \brief Compute the composition \a Q(R).
	 *
	 * Given \a Q and \a R, two polynomials over GF(p), compute
	 * their composition. This implements algorithm \c Compose of
	 * [\ref ISSAC "DFS '09"].
	 *
	 *  \param [out] res A polynomial to hold the result.
	 *  \param [in] Q A polynomial.
	 *  \param [in] R A polynomial.
	 *  \param [in] p The characteristic of the base field of \a Q and \a R. Observe that NTL's
	 *  			modulus has to be set accordingly.
	 */
	template <class T> void compose
	(typename T::GFpX& res, const typename T::GFpX& Q,
	const typename T::GFpX& R, const typename T::BigInt& p);
	/**
	 *  \brief Compute the <i>n</i>-th cyclotomic polynomial modulo \a p.
	 *
	 *  It implements the algorithm described in [\ref Brent "Brent '93"].
	 *
	 *  \param [out] res A polynomial to hold the result.
	 *  \param [in] n A strictly positive integer.
	 *  \param [in] p The characteristic to work with. Observe that NTL's
	 *  			modulus has to be set accordingly.
	 */
	template <class T> void cyclotomic
	(typename T::GFpX& res, const long n, const typename T::BigInt& p);

	/**
	 * @}
	 * \endcond
	 */
}

#endif /*UTILITIES_H_*/
