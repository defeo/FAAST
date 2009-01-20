namespace AS {
/****************** Arithmetics ******************/
	/* n-th iterated frobenius. This algorithm is detailed in the
	 * long version.
	 */
	template <class T> void
	FieldElement<T>::self_frobenius(const long n) throw() {
		if (isZero()) return;
		if (base) return;
		
		//TODO
	}
	
	template <class T> void
	FieldElement<T>::self_trace(const Field<T> F)
	throw(NotASubFieldException) {
		// zero stays zero
		if (!parent_field) {
			*this = F.zero();
		}
		// handle errors
		else if (!F.isSubFieldOf(*parent_field)) {
			throw NotASubFieldException();
		}
		// if F is a prime field, this is the absolute trace
		else if (F.isPrimeField()) {
			self_trace();
		}
		// else, the trace is minus the coefficient of x_i^(p-1)
		else {
			// move to the stem
			parent_field = parent_field->stem;
			
			vector<FieldElement<T> > down;
			while (parent_field != F.stem) {
				pushDown(*this, down);
				*this = down[long(parent_field->p) - 1];
				negate();
			}
			// move out of the stem
			parent_field = &F;
		}
	}
	
	template <class T> void
	FieldElement<T>::self_pseudotrace(const long n) throw() {
		if (isZero()) return;
		if (base) {
			*this *= parent_field->scalar(n);
			return;
		}
		
		//TODO
	}

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
		const Field<T>* parent = parent_field;
		liftUp(result, *this);
		parent_field = parent;
	}

	/* n-th iterated frobenius, n < d */
	template <class T> void FieldElement<T>::SmallFrob(const long n) {
#ifdef AS_DEBUG
		if (n >= parent_field->d)
			throw ASException("Bad input to SmallFrob.");
#endif
		for (long i = 0 ; i < n ; i++)
			self_frobenius();
	}
	
	/* p^j-th pseudotrace */
	template <class T> void FieldElement<T>::BigPTrace(const long j) {
#ifdef AS_DEBUG
		if (j < 0) throw ASException("Bad input to BigPTrace.");
#endif
		if (j == 0) SmallPTrace(parent_field->d);
		else {
			FieldElement<T> t(*this);
			t.BigPTrace(j-1);
			*this += t;
			for (BigInt i = 1 ; i < parent_field->p ; i++) {
				t.BigFrob(j-1);
				*this += t;
			}
		}
	}

	/* n-th pseudotrace, n < d */
	template <class T> void FieldElement<T>::SmallPTrace(const long n) {
#ifdef AS_DEBUG
		if (n > parent_field->d)
			throw ASException("Bad input to SmallPTrace.");
#endif
		FieldElement<T> t(*this);
		*this += t;
		for (long i = 1 ; i < n ; i++) {
			t.self_frobenius();
			*this += t;
		}
	}
	
}
