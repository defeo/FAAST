/*
	This file is part of the FAAST library.

	Copyright (c) 2009 Luca De Feo and Éric Schost.

	The most recent version of FAAST is available at http://www.lix.polytechnique.fr/~defeo/FAAST

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; see file COPYING. If not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
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
