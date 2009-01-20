namespace AS {
/****************** Arithmetics ******************/
	/* n-th iterated frobenius. This algorithm is detailed in the
	 * long version.
	 */
	template <class T> void
	FieldElement<T>::self_frobenius(const long n) throw() {
		if (!parent_field) return;
		if (base) return;
		
	}
	
//	void self_trace(const Field<T> F) throw(NotASubFieldException);
//	void self_pseudotrace(const long n) throw();

/****************** Helpers for frobenius and trace ******************/
	/* See Section 5 and the long version (probably section 6). */
	
	/* p^j-th iterated frobenius */
	template <class T> void FieldElement<T>::BigFrob(const long j) {
#ifdef AS_DEBUG
		if (j < 0) throw ASException("Bad input to BigFrob.");
#endif
		BigInt p = parent_field->p;
		
		// step 1
		if (j >= parent_field->height) return;
		// step 2
		vector<FieldElement<T> > down;
		pushDown(*this, down);
		// step 3
		for (BigInt i = 0 ; i < p ; i++)
			down[i].BigFrob(j);
		// step 5
		const FieldElement<T>& beta = parent_field->getPseudotrace(j);
		vector<FieldElement<T> > result;
		for (BigInt i = 0 ; i <= p ; i++) {
			result[i] = down[p-1];
			for (BigInt j = long(i) - 1 ; j >= 0 ; j--) {
				result[j] *= beta;
				if (j > 0) result[j] += result[j-1];
				else result[j] += down[i];
			}
		}
		// step 6
		liftUp(result, *this);
	}

	/* n-th iterated frobenius, n < d */
//	void SmallFrob(const long n);
	/* p^j-th pseudotrace */
//	void BigPTrace(const long j);
	/* n-th pseudotrace, n < d */
//	void SmallPTrace(const long n);
	
}
