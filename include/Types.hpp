#ifndef TYPES_H_
#define TYPES_H_

#include <NTL/ZZ.h>
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
		typedef ZZ          BigInt;
		typedef struct {
			GF2EContext P;
		}                   Context;
		typedef GF2XModulus GFpXModulus;
	};
}

#endif /*TYPES_H_*/
