#include "utilities.hpp"

/* This file contains algorithms from Sections 3 and 6
 * of the paper
 */

namespace AS {
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
