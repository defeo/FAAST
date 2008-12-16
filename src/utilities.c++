#include <NTL/ZZ.h>
#include <cmath>

namespace AS {
	/* return the gratest power of p less than n */
	long NumPits(const long p, const long n) {
		if (p == 2) return NumBits(n);
		else {
			long k = floor(log(n) / log(p));
		}
	}
}