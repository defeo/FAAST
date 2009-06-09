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

/* Courtesy of E. Schost */

#include <NTL/GF2X.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>

namespace NTL_NAMESPACE {
	/** \ingroup NTLhacks
	 *  \brief Transposed multiplication.
	 *
	 *  This declaration is defined in NTL sources.
	 */
	void TransMulMod(GF2X& x, const GF2X& a, const GF2XTransMultiplier& B, const GF2XModulus& F);

	/** \ingroup NTLhacks
	 *  \brief Transposed multiplication.
	 *  \author É. Schost
	 */
	inline void TransMulMod(ZZ_pX& x, const ZZ_pX& a,
	const ZZ_pXMultiplier& B, const ZZ_pXModulus& F){
	  x = 0;
	  vec_ZZ_p X;
	  UpdateMap(X, a.rep, B, F);
	  for (long i = 0; i < X.length(); i++)
	    SetCoeff(x,i,X[i]);
	  x.normalize();
	}

	/** \ingroup NTLhacks
	 *  \brief \copybrief TransMulMod()
	 *
	 *  \copydetails TransMulMod()
	 */
	inline void TransMulMod(zz_pX& x, const zz_pX& a,
	const zz_pXMultiplier& B, const zz_pXModulus& F){
	  x = 0;
	  vec_zz_p X;
	  UpdateMap(X, a.rep, B, F);
	  for (long i = 0; i < X.length(); i++)
	    SetCoeff(x,i,X[i]);
	  x.normalize();
	}
}
