#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <NTL/ZZ.h>
#include <cmath>

namespace AS {
	/* return the gratest power of p less than or equal to n */
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
}

#endif /*UTILITIES_H_*/
