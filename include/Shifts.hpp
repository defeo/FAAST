#include <NTL/GF2X.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>

NTL_OPEN_NNS
void RightShiftAdd(GF2X& c, const GF2X& a, long n);
NTL_CLOSE_NNS

namespace NTL {
	void RightShiftAdd(ZZ_pX& U, const ZZ_pX& V, long n)
	// assumes input does not alias output
	{
	   if (IsZero(V))
	      return;
	
	   long du = deg(U);
	   long dv = deg(V);
	
	   long d = max(du, dv-n);
	
	   U.rep.SetLength(d+1);
	   long i;
	
	   for (i = du+1; i <= d; i++)
	      clear(U.rep[i]);
	
	   for (i = n; i <= dv; i++)
	      add(U.rep[i-n], U.rep[i-n], V.rep[i]);
	
	   U.normalize();
	}
	
	void RightShiftAdd(zz_pX& U, const zz_pX& V, long n)
	// assumes input does not alias output
	{
	   if (IsZero(V))
	      return;
	
	   long du = deg(U);
	   long dv = deg(V);
	
	   long d = max(du, n-dv);
	
	   U.rep.SetLength(d+1);
	   long i;
	
	   for (i = du+1; i <= d; i++)
	      clear(U.rep[i]);
	
	   for (i = n; i <= dv; i++)
	      add(U.rep[i-n], U.rep[i-n], V.rep[i]);
	
	   U.normalize();
	}
}
