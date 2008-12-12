/* This file contains algorithms from Sections 3 and 6 of the paper */

namespace AS {
	// Return Q(X^p-X)
	template <class T> void compose_Xp_minus_X(GFpX& res, const GFpX& Q) {
		
	}
	// Return Q(X-1)
	template <class T> void compose_X_minus_1(GFpX& res, const GFpX& Q) {
		GFpX Q0, Q1, Q0X, Q1X;
		long degree = deg(f);
		long k = (p == 2) ? NumBits(degree) : //TODO;
		if (k > 1) {
			long halfdegree = power_long(2, k-1);
			for (long i = 0 ; i < halfdegree ; i++)
				SetCoeff(f0, i, coeff(f, i));
			for (long i = halfdegree ; i <= degree ; i++)
				SetCoeff(f1, i - halfdegree, coeff(f, i));
				
			f0X = composeXminus1(f0);
			f1X = composeXminus1(f1);
			result = f0X + f1X + LeftShift(f1X, halfdegree);
		} else {
			GF2X tmp; SetCoeff(tmp,0); SetCoeff(tmp,1);
			result = coeff(f, 0); 
			if (coeff(f, 1)==1) result += tmp;
		}
		
		return result;
	}

/****************** Field Extensions ******************/
	/* Build a default extension of degree p over this field. 
	 */
	template <class T> Field<T>& Field<T>::ArtinSchreierExtension()
	const throw (CharacteristicTooLargeExeption, NotSupportedException) {
#ifdef AS_DEBUG
		if (!stem) throw ASException("Error : Stem is NULL.");
#endif
		// if the extension already exists, return it
		if (stem->overfield) return *(stem->overfield);

		// Build the extension (Section 3)

		// test if the characteristic stays in one word
		if (p > to_long(p)) throw CharacteristicTooLargeException();
		//TODO
	}
	
	/* Build the splitting field of the polynomial
	 * 			X^p - X - alpha
	 * over this field. This may or may not be an extension of
	 * degree p depending if the polynomial is irreducible.
	 */
//	Field<T>& ArtinSchreierExtension(const FieldElement<T>& alpha) const throw (CharacteristicTooLargeExeption, NotSupportedException);


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
