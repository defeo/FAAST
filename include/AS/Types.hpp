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

namespace AS {
	/**
	 * \if DEV
	 * \ingroup NTLhacks
	 * \endif
	 * 
	 * \brief A class providing automatic conversion between
	 * \c ZZ and \c long types.
	 * 
	 * You can pass a \c ZZ as well as a \c long to
	 * any function requiring a \c ZZ_auto as parameter.
	 * 
	 * This class is not meant to be used by the user. Use \c ZZ, \c long
	 * or ZZ_p_Algebra::BigInt instead. 
	 */
	class ZZ_auto : public ZZ {
		ZZ_auto(const long i) : ZZ(to_ZZ(i)) {}
	public:
		ZZ_auto() {};
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
	 * Infrastructures are collections of types providing genericity over
	 * NTL types. They provide an unique set of names to abstract from the
	 * implementation details of the three NTL families \c ZZ_p*, \c zz_p* and \c GF2*.
	 * 
	 * Most routines in the library take an Infrastructure as template parameter.
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

		/** \brief \copybrief zz_p_Algebra::name */
		static char const * const name;
	};
	
	/** @} */
}

#endif /*TYPES_H_*/
