namespace AS {
	/* return the least power of p greater than n */
	long NumPits(const long p, long n) {
		if (p == 2) return NumBits(n);
		else {
			if (n == 0) return 0;
			if (n < 0) n = -n;
			long k = floor(log(n) / log(p)) + 1;
			long pk = power_long(p, k);
			if (pk <= n) {
				long i = 0;
				while (pk <= n) { i++; pk *= p; }
#if AS_DEBUG>=2
				if (i > 1) cout << "cmath log imprecise by " << i << " p-its." << endl;
#endif
				return k+i;
			} else {
				long i = -1;
				while (pk > n) { i++; pk /= p; }
#if AS_DEBUG>=2
				if (i > 1) cout << "cmath log imprecise by " << i << " p-its." << endl;
#endif
				return k - i;
			}
		}
	}

	/* Store in res the composition Q(R).
	 * Q and R are two polynomials over GF(p)
	 */
	template <class T> void compose
	(typename T::GFpX& res, const typename T::GFpX& Q,
	const typename T::GFpX& R, const typename T::BigInt p) {
		typedef typename T::GFpX GFpX;

		long degree = max(deg(Q),0);
		long degR = max(deg(R),0);
		long k = NumPits(p, degree);
		if (k > 0) {
			res = 0;
			long splitdegree = power_long(p, k-1);
			for (long i = splitdegree * (degree / splitdegree) ; i >= 0 ; i -= splitdegree) {
				GFpX Q1;
				for (long j = 0 ; j < splitdegree && i+j <= degree ; j++) {
					SetCoeff(Q1, j, coeff(Q, i+j));
				}
				GFpX Q1X; compose<T>(Q1X, Q1, R, p);
				// Horner's rule
				GFpX shifted;
				if (res != 0) {
					for (long h = 0 ; h <= degR ; h++) {
						if (coeff(R, h) != 0)
							shifted += coeff(R, h) * LeftShift(res, splitdegree*h);
					}
				}
				res = Q1X + shifted;
			}
		} else {
			res = Q;
		}
	}
	
	/* Brent variant of Pollard Rho. Returns 0 or a
	 * factor of n.
	 */
	long brent80(long n, long x0, long m, long c) {
/*		ZZ x, y, d, zzn, zzc;
		x = y = x0; d = 1; zzn = n; zzc = c;
		
		while (d == 1) {
			PowerMod(x,x,2,zzn); AddMod(x, x, zzc, zzn);
			PowerMod(y,y,2,zzn); AddMod(y, y, zzc, zzn);
			PowerMod(y,y,2,zzn); AddMod(y, y, zzc, zzn);
			if (x == y) return 0;
			d = GCD(abs(x-y), zzn);
		}
		return to_long(d);*/
		
		ZZ y, x, ys, zzn, zzc, q, g;
		y = x0; zzn = n; zzc = c; q = 1; g = 0;
		long r = 1;

		while (g <= 1) {
			x = y;
			for (long i = 1 ; i <= r ; i++)
				PowerMod(y,y,2,zzn); AddMod(y, y, zzc, zzn);
			long k = 0;
			while(k < r && g <= 1) {
				ys = y;
				for (long i = 1 ; i <= min(m, r-k) ; i++) {
					PowerMod(y,y,2,zzn); AddMod(y, y, zzc, zzn);
					MulMod(q, q, abs(x-y), zzn);
				}
				g = GCD(q, zzn);
				k += m;
			}
			r *= 2;
		}
		if (g == zzn) {
			ZZ gg;
			while (gg <= 1) {
				PowerMod(ys,ys,2,zzn); AddMod(ys, ys, zzc, zzn);
				gg = GCD(abs(x-ys), zzn);
			}
			if (gg == n) return 0;
			else return to_long(gg);
		}
		else return to_long(g);
	}
	
	// Deterministic primality test. Only valid if
	// 17 < n < 341,550,071,728,321 
	bool isJaeschkePrime(long n) {
		ZZ witn[7];
		witn[0] = 2; witn[1] = 3; witn[2] = 5;
		witn[3] = 7; witn[4] = 11; witn[5] = 13;
		witn[6] = 17;
		
		ZZ zzn; zzn = n;
		for (int i = 0 ; i < 7 ; i++) {
			if (MillerWitness(zzn, witn[i])) return false;
		}
		return true;
	}
	
	/* Inserts n with multiplicity m in the vector */
	void insertInPlace(vector<pair<long,int> >& factors, 
	long n, int m) {
		bool inserted = false;
		vector<pair<long,int> >::iterator it;
		for (it = factors.begin() ;
			it != factors.end() && !inserted ; it++) {
			if (it->first == n) {
				it->second += m;
				inserted = true;
			}
			if (it->first > n) {
				factors.insert(it, pair<long,int>(n,m));
				inserted = true;
			}
		}
		if (!inserted) factors.insert(it, pair<long,int>(n,m));
	}

	/* Pollard Rho factorisation algorithm */
	void factor(long n, vector<pair<long,int> >& factors,
	const bool trial, const int multiplicity) {
		int mul;
		if (trial) {
			// trial division up to 17
			long trials[] = {2, 3, 5, 7, 11, 13, 17};
			for (int i = 0 ; i < 7 && n > 1; i++) {
				for (mul = 0 ; n % trials[i] == 0 ; mul++) n /= trials[i];
				if (mul > 0) {
					insertInPlace(factors, trials[i], 
									mul*multiplicity);
				}
			}
		}
		
		if (n <= 1) return;
		// test primality
		bool isprime = isJaeschkePrime(n);
		if (isprime && NumBits(n) > 47) {
			isprime &= ProbPrime(n, 10);
		}
		
		// if n is prime, simply add it to the factorisation
		if (isprime) {
			insertInPlace(factors, n, multiplicity);
		}
		// else, use pollard rho until a factor is found
		else {
			long fact = 0;
			while (fact == 0) {
				long c = 0;
				while (c == 0) { c = RandomBnd(n); }
				fact = brent80(n,1,1,c);
			}
			for (mul = 0 ; n % fact == 0 ; mul++) n /= fact;
			factor(fact, factors, false, mul*multiplicity);
			factor(n, factors, false, multiplicity);
		}
	}
	
	/* compute the n-th cyclotomic polynomial modulo p */
	template <class T> void cyclotomic
	(typename T::GFpX& res, const long n, const typename T::BigInt p) {
		typedef typename T::GFpX GFpX;

		if (n == 1) { SetX(res); SetCoeff(res, 0, -1); }
		else {
			vector<pair<long,int> > factors; factor(n, factors);
			// first process the squarefree factors of n
			GFpX Phi; cyclotomic<T>(Phi, 1, p);
			long m = 1;  // this variable contains the non-squarefree part
			vector<pair<long,int> >::iterator it;
			for (it = factors.begin() ; it != factors.end() ; it++) {
				m *= power_long(it->first, (it->second - 1));
				long q = it->first;
				GFpX Phip; compose<T>(Phip, Phi, GFpX(q,1), p);
#ifdef AS_DEBUG
				if (Phip % Phi != 0)
					throw ASException("Error : computing the cyclotomic polynomial.");
#endif
				Phi = Phip / Phi;
			}
			
			// now add the other factors
			if (m > 1) {
				compose<T>(res, Phi, GFpX(m,1), p);
			} else res = Phi;
		}
	}

}
