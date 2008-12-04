#ifndef TYPES_H_
#define TYPES_H_

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
	struct ZZ_p_Algebra {
		typedef ZZ_p        GFp;
		typedef mat_ZZ_p    MatGFp;
		typedef ZZ_pX       Poly;
		typedef ZZ_pE       ModPoly;
		typedef ZZ_pEX      PolyModPoly;
		typedef ZZ          BigInt;
		typedef struct {
			ZZ_pContext  p;
			ZZ_pEContext P;
		}                   Context;
	};

	struct GF2_Algebra {
		typedef GF2         GFp;
		typedef mat_GF2     MatGFp;
		typedef GF2X        Poly;
		typedef GF2E        ModPoly;
		typedef GF2EX       PolyModPoly;
		typedef ZZ          BigInt;
		typedef GF2EContext Context;
	};
}

#endif /*TYPES_H_*/
