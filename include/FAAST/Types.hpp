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
#ifndef TYPES_H_
#define TYPES_H_

#include <NTL/ZZ.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2.h>
#include "NTLhacks.hpp"

NTL_CLIENT

namespace FAAST {
	/** \brief A structure to hold predefined constants depending on the type. */
	struct Constants {
		const long HalfGCD_CROSSOVER;
		Constants(const long h) : HalfGCD_CROSSOVER(h) {}
	};

	/**
	 * \if DEV
	 * \ingroup NTLhacks
	 * \endif
	 *
	 * \brief A class providing automatic conversion between
	 * 22\c ZZ and \c long types.
	 *
	 * You can pass a \c ZZ as well as a \c long to
	 * any function requiring a \c ZZ_auto as parameter.
	 *
	 * This class is not meant to be used by the user. Use \c ZZ, \c long
	 * or ZZ_p_Algebra::BigInt instead.
	 */
	class ZZ_auto : public ZZ {
	public:
		ZZ_auto() {}
		ZZ_auto(const long i) : ZZ(to_ZZ(i)) {}
		ZZ_auto(const ZZ& i) : ZZ(i) {}
		ZZ_auto& operator=(const ZZ_auto& i) { ZZ::operator=(i); return *this; }
		ZZ_auto& operator=(const ZZ& i) { ZZ::operator=(i); return *this; }
		ZZ_auto& operator=(const long i) { ZZ::operator=(i); return *this; }
		operator long() const {
			return to_long(*this);
		}
	};

	/**
	 * \defgroup Infrastructures Infrastructures
	 * \NTL provides three different ways of representing modular integers:
	 *  - \c zz_p is a class representing modular integers with word-sized
	 *  modulus,
	 *  - \c ZZ_p is a class representing modular integers with arbitrary sized
	 *  modulus,
	 *  - \c GF2 is a class representing integers modulo 2.
	 *
	 * To each of these types corresponds a whole family of types (\c ZZ_p_X,
	 * \c GF2_E, etc.) representing
	 * polynomials with modular coefficients, modular polynomials with modular
	 * coefficients, etc. The actual algorithms implementing arithmetics
	 * for such types vary and affect remarkably the performances of FAAST.
	 * See the \NTL manual for more details.
	 *
	 * \ref Infrastructures are collections (\c struct 's) of types providing
	 * genericity over
	 * \NTL types. They provide an unique set of names to abstract from the
	 * implementation details of the three NTL families \c ZZ_p*, \c zz_p*
	 * and \c GF2*. You must provide an Infrastructure as template parameter to
	 * the types and most of the functions of FAAST. This parameter tells FAAST
	 * which of the three \NTL families it should use to perform modular
	 * arithmetics.
	 *
	 * Here's an example of how to use the infrastructure FAAST::zz_p_Algebra.
	 * First some \c typedef 's to save typing:
	 * \dontinclude test.c++
	 * \skipline typedef
	 * \line typedef
	 * \line typedef
	 * Then define some parameters (characteristic, degree, etc.),
	 * notice the use of FAAST::Field::Infrastructure as an alias for
	 * FAAST::zz_p_Algebra (useful if you later change your mind about
	 * the Infrastructure):
	 * \skipline gfp
	 * \line long
	 * Finally create a finite field:
	 * \skipline K
	 *
	 * If you are wondering which Infrastructure you should use, then use
	 * FAAST::zz_p_Algebra, as it is quite flexible and way faster than
	 * FAAST::ZZ_p_Algebra.
	 *
	 * If you plan to construct fields with huge characteristics (larger than the
	 * largest \c long),
	 * then you should opt for FAAST::ZZ_p_Algebra; notice however that FAAST does
	 * not let you build Artin-Schreier extensions in characteristics greater than
	 * the greatest \c long, so you will probably miss its most exciting features.
	 *
	 * Finally, if you only work in characteristic 2 and care about
	 * performance, you should consider compiling \NTL with the \gf2x library
	 * and using FAAST::GF2_Algebra as Infrastructure. If you don't use the
	 * \gf2x library, then FAAST::GF2_Algebra will only be interesting for moderate
	 * field cardinalities, but it will give the pace up to FAAST::zz_p_Algebra
	 * for huge fields.
	 * @{
	 */

	/**
	 * \brief The infrastructure corresponding to \c zz_p* types.
	 *
	 * NTL's \c zz_p* types provide arithmetics modulo \a p where \a p
	 * is a word-size integer greater than one.
	 */
	struct zz_p_Algebra {
		/** \brief Elements of the field F<sub>p</sub>. */
		typedef zz_p            GFp;
		/** \brief Matrices over the field F<sub>p</sub>. */
		typedef mat_zz_p        MatGFp;
		/** \brief Vectors over the field F<sub>p</sub>. */
		typedef vec_zz_p        VecGFp;
		/** \brief Polynomials over the field F<sub>p</sub>. */
		typedef zz_pX           GFpX;
		/** \brief Elements of an extension field of F<sub>p</sub>. */
		typedef zz_pE           GFpE;
		/** \brief Polynomials over an extension field of F<sub>p</sub>. */
		typedef zz_pEX          GFpEX;
		/** \brief The type of the characteristic \a p. */
		typedef long            BigInt;
		/**
		 * \brief Modulus switching data.
		 *
		 * NTL's modulus switching mechanism permits to save the data defining a finite
		 * field. This type groups togheter all these data.
		 */
		typedef struct {
			/** \brief The characteristic of a base field. */
			zz_pContext  p;
			/** \brief The defining polynomial of an extension field. */
			zz_pEContext P;
		}                       Context;
		/** \brief Pre-condtioned polynomials over the field F<sub>p</sub>. */
		typedef zz_pXModulus    GFpXModulus;
		/** \brief Pre-condtioned polynomials over the field F<sub>p</sub>. */
		typedef zz_pXMultiplier GFpXMultiplier;

		/** \brief predefined constants */
		static const Constants consts;
		/** \brief The name of the Infrastructure */
		static char const * const name;
	};

	/**
	 * \brief The infrastructure corresponding to \c ZZ_p* types.
	 *
	 * NTL's \c ZZ_p* types provide arithmetics modulo \a p where \a p
	 * is a multiprecision integer greater than one.
	 */
	struct ZZ_p_Algebra {
		/** \brief \copybrief zz_p_Algebra::GFp */
		typedef ZZ_p            GFp;
		/** \brief \copybrief zz_p_Algebra::MatGFp */
		typedef mat_ZZ_p        MatGFp;
		/** \brief \copybrief zz_p_Algebra::VecGFp */
		typedef vec_ZZ_p        VecGFp;
		/** \brief \copybrief zz_p_Algebra::GFpX */
		typedef ZZ_pX           GFpX;
		/** \brief \copybrief zz_p_Algebra::GFpE */
		typedef ZZ_pE           GFpE;
		/** \brief \copybrief zz_p_Algebra::GFpEX */
		typedef ZZ_pEX          GFpEX;
		/** \brief \copybrief zz_p_Algebra::BigInt */
		typedef ZZ_auto         BigInt;
		/**
		 * \brief \copybrief zz_p_Algebra::Context
		 *
		 * \copydetails zz_p_Algebra::Context
		 */
		typedef struct {
			/** \brief \copybrief zz_p_Algebra::Context::p */
			ZZ_pContext  p;
			/** \brief \copybrief zz_p_Algebra::Context::P */
			ZZ_pEContext P;
		}                       Context;
		/** \brief \copybrief zz_p_Algebra::GFpXModulus */
		typedef ZZ_pXModulus    GFpXModulus;
		/** \brief \copybrief zz_p_Algebra::GFpXMultiplier */
		typedef ZZ_pXMultiplier GFpXMultiplier;

		/** \brief \copybrief zz_p_Algebra::consts */
		static const Constants consts;
		/** \brief \copybrief zz_p_Algebra::name */
		static char const * const name;
	};

	/**
	 * \brief The infrastructure corresponding to \c GF2* types.
	 *
	 * NTL's \c GF2* types provide highly optimized arithmetics modulo 2.
	 */
	struct GF2_Algebra {
		/** \brief \copybrief zz_p_Algebra::GFp */
		typedef GF2                 GFp;
		/** \brief \copybrief zz_p_Algebra::MatGFp */
		typedef mat_GF2             MatGFp;
		/** \brief \copybrief zz_p_Algebra::VecGFp */
		typedef vec_GF2             VecGFp;
		/** \brief \copybrief zz_p_Algebra::GFpX */
		typedef GF2X                GFpX;
		/** \brief \copybrief zz_p_Algebra::GFpE */
		typedef GF2E                GFpE;
		/** \brief \copybrief zz_p_Algebra::GFpEX */
		typedef GF2EX               GFpEX;
		/** \brief \copybrief zz_p_Algebra::BigInt */
		typedef int                 BigInt;
		/**
		 * \brief \copybrief zz_p_Algebra::Context
		 *
		 * \copydetails zz_p_Algebra::Context
		 */
		typedef struct {
			/** \brief \copybrief zz_p_Algebra::Context::P */
			GF2EContext P;
		}                           Context;
		/** \brief \copybrief zz_p_Algebra::GFpXModulus */
		typedef GF2XModulus         GFpXModulus;
		/** \brief \copybrief zz_p_Algebra::GFpXMultiplier */
		typedef GF2XTransMultiplier GFpXMultiplier;

		/** \brief \copybrief zz_p_Algebra::consts */
		static const Constants consts;
		/** \brief \copybrief zz_p_Algebra::name */
		static char const * const name;
	};

	/** @} */
}

#endif /*TYPES_H_*/
