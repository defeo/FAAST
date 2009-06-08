/* Courtesy of E. Schost */

#include <NTL/GF2X.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>

namespace NTL_NAMESPACE {
	/** \ingroup NTLhacks
	 *  \brief Transposed multiplication.
	 *
	 *  This is taken directly from NTL sources.
	 *
	 *  \author V. Shoup
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
