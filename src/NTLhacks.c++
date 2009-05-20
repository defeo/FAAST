namespace NTL {
	void power(zz_p& x, const zz_p& a, const ZZ& e) {
		if (e == 0) x = 1;
		else if (a == 0 && e < 0)
			Error("zz_p: division by zero");
		else {
			long p = zz_p::modulus();
			long ee = e % (p-1);
			power(x, a, ee);
		}
	}
	void power(GF2& x, const GF2& a, const ZZ& e) {
		if (e == 0) x = 1;
		else if (a == 0 && e < 0)
			Error("GF2: division by zero");
		else x = a;
	}

	// this one is obsolete
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

	// this one is obsolete
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
