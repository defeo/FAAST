namespace AS {
	/* return the gratest power of p less than n */
	long NumPits(const long p, const long n) {
		if (p == 2) return NumBits(n);
		else {
			long k = floor(log(n) / log(p));
			long pk = power_long(p, k);
			if (pk <= n) {
				long i = -1;
				while (pk <= n) { i++; pk *= p; }
#if AS_DEBUG>=2
				if (i > 1) cout << "cmath log imprecise by " << i << " p-its." << endl;
#endif
				return k+i;
			} else {
				long i = 0;
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
	const typename T::GFpX& R, const long p) {
		typedef typename T::GFpX GFpX;
	
		long degree = deg(Q);
		long degR = deg(R);
		long k = NumPits(p, degree);
		if (k > 0) {
			res = 0;
			long splitdegree = power_long(p, k-1);
			for (long i = 0 ; i <= degree ; i += splitdegree) {
				GFpX Q1;
				for (long j = 0 ; j < splitdegree && i+j <= degree ; j++) {
					SetCoeff(Q1, j, coeff(Q, i+j));
				}
				GFpX Q1X; compose(Q1X, Q1, R, p);
				// Horner's rule
				for (long h = 0 ; h <= degR && res != 0 ; h++) {
					if (coeff(R, h) != 0)
						res += coeff(R, h) * LeftShift(res, splitdegree*h);
				}
				res += Q1X;
			}
		} else {
			res = Q;
		}
	}
	
	long brent80(long n, long x0, long m, long c) {
		long y = x0, r = 1, q = 1, g = 0, x, ys;
		while (g <= 1) {
			x = y;
			for (long i = 1 ; i <= r ; i++)
				y = y*y + c;
			long k = 0;
			while(k < r && g <= 1) {
				ys = y;
				for (long i = 1 ; i <= min(m, r-k) ; i++) {
					y = y*y + c;
					q = (q * abs(x-y)) % n;
				}
				g = GCD(q, n);
				k = k + m;
			}
			r = 2*r;
		}
		if (g == n) {
			long gg = 0;
			while (gg <= 1) {
				ys = ys*ys + c;
				gg = GCD(abs(x-ys), n);
			}
			if (gg == n) return 0;
			else return gg;
		}
		else return g;
	}
	
	/* Pollard Rho factorisation algorithm */
	void factor(const long n, vector<pair<long,int> >& factors) {
		
	}
	
	/* compute the n-th cyclotomic polynomial modulo p */
	template <class T> void cyclotomic
	(typename T::GFpX& res, const long n, const long p) {
		typedef typename T::GFpX GFpX;

		if (n == 1) { SetX(res); setCoeff(res, 0, -1); }
		else {
			vector<pair<long,int> > factors; factor(n, factors);
			// first process the squarefree factors of n
			GFpX Phi; cyclotomic(Phi, 1, p);
			long m = 1;  // this variable contains the non-squarefree part
			for (long i = 0 ; i < factors.size() ; i++) {
				m *= power_long(factors[i].first, (factors[i].second - 1));
				long q = factors[i].first;
				GFpX Phip; compose(Phip, Phi, GFpX(q,1), p);
#ifdef AS_DEBUG
				if (Phip % Phi != 0)
					throw ASException("Error : computing the cyclotomic polynomial.");
#endif
				Phi = Phip / Phi;
			}
			// now add the other factors
			if (m > 1) {
				compose(res, Phi, GFpX(m,1), p);
			} else res = Phi;
		}
	}

}