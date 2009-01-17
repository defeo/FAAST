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
	class ZZ_auto : public ZZ {
	public:
		ZZ_auto() {};
		ZZ_auto(const long i) : ZZ(to_ZZ(i)) {}
		ZZ_auto(const ZZ& i) : ZZ(i) {}
		ZZ_auto& operator=(const ZZ_auto& i) { ZZ::operator=(i); return *this; }
		ZZ_auto& operator=(const ZZ& i) { ZZ::operator=(i); return *this; }
		ZZ_auto& operator=(const long i) { ZZ::operator=(i); return *this; }
		operator long() const {
			return to_long(*this);
#if AS_DEBUG >= 2
			cout << "Warning : possible loss of precision." << endl;
#endif
		}
	};

	struct zz_p_Algebra {
		typedef zz_p            GFp;
		typedef mat_zz_p        MatGFp;
		typedef zz_pX           GFpX;
		typedef zz_pE           GFpE;
		typedef zz_pEX          GFpEX;
		typedef long            BigInt;
		typedef struct {
			zz_pContext  p;
			zz_pEContext P;
		}                       Context;
		typedef zz_pXModulus    GFpXModulus;
		typedef zz_pXMultiplier GFpXMultiplier;

		static char const * const name;
	};
	char const * const zz_p_Algebra::name = "zz_p";

	struct ZZ_p_Algebra {
		typedef ZZ_p            GFp;
		typedef mat_ZZ_p        MatGFp;
		typedef ZZ_pX           GFpX;
		typedef ZZ_pE           GFpE;
		typedef ZZ_pEX          GFpEX;
		typedef ZZ_auto         BigInt;
		typedef struct {
			ZZ_pContext  p;
			ZZ_pEContext P;
		}                       Context;
		typedef ZZ_pXModulus    GFpXModulus;
		typedef ZZ_pXMultiplier GFpXMultiplier;

		static char const * const name;
	};
	char const * const ZZ_p_Algebra::name = "ZZ_p";

	struct GF2_Algebra {
		typedef GF2                 GFp;
		typedef mat_GF2             MatGFp;
		typedef GF2X                GFpX;
		typedef GF2E                GFpE;
		typedef GF2EX               GFpEX;
		typedef int                 BigInt;
		typedef struct {
			GF2EContext P;
		}                           Context;
		typedef GF2XModulus         GFpXModulus;
		typedef GF2XTransMultiplier GFpXMultiplier;

		static char const * const name;
	};
	char const * const GF2_Algebra::name = "GF2";
}

#endif /*TYPES_H_*/
