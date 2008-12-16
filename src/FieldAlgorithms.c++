#include "utilities.hpp"

/* This file contains algorithms from Sections 3 and 6 of the paper */

namespace AS {
	// Return Q(X^p-X)
	template <class T> void compose_Xp_minus_X
	(typename T::GFpX& res, const typename T::GFpX& Q, const typename T::BigInt& p) {
		
	}
	// Return Q(X-1)
	template <class T> void compose_X_minus_1
	(typename T::GFpX& res, const typename T::GFpX& Q, const typename T::BigInt& p) {
		typedef typename T::GFpX GFpX;
	
		long degree = deg(Q);
		long k = NumPits(p, degree);
		if (k > 0) {
			res = 0;
			long splitdegree = power_long(p, k-1);
			for (long i = 0 ; i <= degree ; i += splitdegree) {
				GFpX Q1;
				for (long j = 0 ; j < splitdegree && i+j <= degree ; j++) {
					SetCoeff(Q1, j, coeff(Q, i+j));
				}
				GFpX Q1X; compose_X_minus_1(Q1X, Q1);
				// Horner's rule
				res = Q1X + LeftShift(res, splitdegree) - res;
			}
		} else {
			res = Q;
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
