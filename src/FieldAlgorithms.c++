#include "utilities.hpp"

/* This file contains algorithms from Sections 3 and 6 of the paper */

namespace AS {
	// compute Q(R)
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
	
	/* Pollard Rho */
	void factor(const long n, vector<pair<long,int> >& factors) {
		
	}
	
	// compute the n-th cyclotomic polynomial
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

/****************** Field Extensions ******************/
	/* Build a default extension of degree p over this field. 
	 */
	template <class T> const Field<T>& Field<T>::ArtinSchreierExtension()
	const throw (CharacteristicTooLargeException, NotSupportedException) {
#ifdef AS_DEBUG
		if (!stem) throw ASException("Error : Stem is NULL.");
#endif
		// if the extension already exists, return it
		if (stem->overfield) return *(stem->overfield);

		// Build the extension (Section 3)

		// test if the characteristic stays in one word
		if (p != long(p)) throw CharacteristicTooLargeException();
		//TODO
	}
	
	/* Build the splitting field of the polynomial
	 * 			X^p - X - alpha
	 * over this field. This may or may not be an extension of
	 * degree p depending if the polynomial is irreducible.
	 */
//	const Field<T>& ArtinSchreierExtension(const FieldElement<T>& alpha) const throw (CharacteristicTooLargeExeption, NotSupportedException);


/****************** Level embedding ******************/
	/* Push the element e down to this field and store
	 * the result in v.
	 * 
	 * throw : IllegalCoercionException if the field e belongs to
	 *         is not the immediate overfield of this.
	 */
//	void pushDown(FieldElement<T>& e, vector<FieldElement<T> >& v) const
//		 throw(IllegalCoercionException);

	/* Lift the elements in v up to this field and store the result in e.
	 * If v is too short, it is filled with zeros. If v is too
	 * long, the unnecessary elements are ignored.
	 * 
	 * throw : NotInSameFieldException if the elements of v do not
	 *         belong all to the same field.
	 * throw : IllegalCoercionException if the field e belongs to
	 *         is not the immediate subfield of this.
	 */
//	void liftUp(vector<FieldElement<T> >& v, FieldElement<T>& e) const
//		throw(NotInSameFieldException, IllegalCoercionException);

}
