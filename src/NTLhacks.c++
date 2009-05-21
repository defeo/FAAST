namespace NTL_NAMESPACE {
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
}
