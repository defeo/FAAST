/* This header contains small hacks adding routines and features to NTL.
 *
 * NTLhacks.hpp
 *
 *  Created on: May 19, 2009
 *      Author: defeo
 */

#ifndef NTLHACKS_HPP_
#define NTLHACKS_HPP_

#include <NTL/GF2X.h>
#include <NTL/GF2.h>
#include <NTL/lzz_p.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/tools.h>

namespace NTL_NAMESPACE {
	/**
	 * \cond DEV
	 * \defgroup NTLhacks Hacks to NTL
	 * This section contains hacks inserted in the NTL namespace
	 * in order to improve the genericity of the code.
	 * @{
	 */

	/**
	 * \brief Transforms a polynomial over GF2 in a monic polynomial.
	 *
	 * Of course, there is nothing to do.
	 */
	inline void MakeMonic(GF2X& x) {}
	/**
	 * \brief Compute \a x = \a a<sup>\a e</sup> (\a e may be negative).
	 *
	 * \warning For the moment, we use the fact that \a p is prime!
	 */
	void power(zz_p& x, const zz_p& a, const ZZ& e);
	/** \brief \copybrief power() */
	void power(GF2& x, const GF2& a, const ZZ& e);

	/** @}
	 *  \endcond
	 */
}

#endif /* NTLHACKS_HPP_ */
