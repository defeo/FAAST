#ifndef TYPES_H_
#define TYPES_H_

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2EX.h>

NTL_CLIENT

namespace AS {
	struct ZZ_p_Algebra {
		typedef ZZ_p   GFp;
		typedef ZZ_pX  Poly;
		typedef ZZ_pE  ModPoly;
		typedef ZZ_pEX PolyModPoly;
		typedef ZZ     BigInt;
	};

	struct GF2_Algebra {
		typedef GF2   GFp;
		typedef GF2X  Poly;
		typedef GF2E  ModPoly;
		typedef GF2EX PolyModPoly;
		typedef ZZ    BigInt;
	};
}

#endif /*TYPES_H_*/
