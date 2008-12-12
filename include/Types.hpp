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

NTL_CLIENT

namespace AS {
	struct zz_p_Algebra {
		typedef zz_p         GFp;
		typedef mat_zz_p     MatGFp;
		typedef zz_pX        GFpX;
		typedef zz_pE        GFpE;
		typedef zz_pEX       GFpEX;
		typedef long         BigInt;
		typedef struct {
			zz_pContext  p;
			zz_pEContext P;
		}                    Context;
		typedef zz_pXModulus GFpXModulus;
	};

	struct ZZ_p_Algebra {
		typedef ZZ_p         GFp;
		typedef mat_ZZ_p     MatGFp;
		typedef ZZ_pX        GFpX;
		typedef ZZ_pE        GFpE;
		typedef ZZ_pEX       GFpEX;
		typedef ZZ           BigInt;
		typedef struct {
			ZZ_pContext  p;
			ZZ_pEContext P;
		}                    Context;
		typedef ZZ_pXModulus GFpXModulus;
	};

	struct GF2_Algebra {
		typedef GF2         GFp;
		typedef mat_GF2     MatGFp;
		typedef GF2X        GFpX;
		typedef GF2E        GFpE;
		typedef GF2EX       GFpEX;
		typedef int          BigInt;
		typedef struct {
			GF2EContext P;
		}                   Context;
		typedef GF2XModulus GFpXModulus;
	};

	template <class T> long to_long(const typename T::BigInt& B) throw()
	{ return B; }
	template <class T> void to_BigInt(typename T::BigInt& B, const long l)
	throw()	{ B = l; }

	template<> inline long to_long<ZZ_p_Algebra>(const ZZ& B) throw()
	{ return to_long(B); }
}

#endif /*TYPES_H_*/
