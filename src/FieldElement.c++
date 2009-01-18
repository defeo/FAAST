namespace AS {
/****************** Copy ******************/
	template <class T>
	FieldElement<T>::FieldElement(const FieldElement<T>& e) throw() :
	repBase(), repExt(), base(e.base), parent_field(e.parent_field) {
		if (parent_field) {
			parent_field->switchContext();
			if (base) repBase = e.repBase;
			else repExt = e.repExt;
		}
	}

	template <class T> FieldElement<T>&
	FieldElement<T>::operator=(const FieldElement<T>& e) throw() {
		repBase = 0;
		repExt = 0;
		base = e.base;
		parent_field = e.parent_field;
		if (parent_field) {
			parent_field->switchContext();
			if (base) repBase = e.repBase;
			else repExt = e.repExt;
		}
		return *this;
	}

	template <class T> FieldElement<T>&
	FieldElement<T>::operator=(const BigInt& i)
	throw(UndefinedFieldException) {
		if (!parent_field) {
			if (i == 0) return *this;
			else throw UndefinedFieldException();
		}
		parent_field->switchContext();
		if (base) repBase = i;
		else repExt = i;
		return *this;
	}
	
/****************** Arithmetics ******************/
	/* Self-incrementing binary operations. */
	template <class T> void
	FieldElement<T>::operator+=(const FieldElement<T>& e)
	throw(NotInSameFieldException) {
		if (!e.parent_field) return;
		if (!parent_field) {
			*this = e;
			return;
		}
		sameLevel(e);
		parent_field->switchContext();
		if (base) repBase += e.repBase;
		else repExt += e.repExt;
	}

	template <class T> void
	FieldElement<T>::operator-=(const FieldElement<T>& e)
	throw(NotInSameFieldException) {
		if (!e.parent_field) return;
		if (!parent_field) {
			*this = e;
			negate();
			return;
		}
		sameLevel(e);
		parent_field->switchContext();
		if (base) repBase -= e.repBase;
		else repExt -= e.repExt;
	}

	template <class T> void
	FieldElement<T>::operator*=(const FieldElement<T>& e)
	throw(NotInSameFieldException) {
		if (!e.parent_field && !parent_field) return;
		if (!parent_field) {
			*this = e.parent_field->zero();
			return;
		}
		if (!e.parent_field) {
			*this = parent_field->zero();
			return;
		}
		sameLevel(e);
		parent_field->switchContext();
		if (base) repBase *= e.repBase;
		else repExt *= e.repExt;
	}

	template <class T> void
	FieldElement<T>::operator/=(const FieldElement<T>& e)
	throw(NotInSameFieldException, DivisionByZeroException) {
		if (e.isZero()) throw DivisionByZeroException();
		if (!parent_field) {
			*this = e.parent_field->zero();
			return;
		}
		sameLevel(e);
		parent_field->switchContext();
		if (base) repBase *= e.repBase;
		else repExt *= e.repExt;
	}

	/* Unary operations */
	/* Absolute trace over GF(p) */
	template <class T> FieldElement<T> FieldElement<T>::trace() const
	throw() {
		if (!parent_field) return *this;
		parent_field->switchContext();
		if (base) return *this;
		else {
			FieldElement<T> tmp;
			NTL::trace(tmp.repBase, repExt);
			tmp.base = true;
			tmp.parent_field = &(parent_field->baseField());
			return tmp;
		}
	}

	/* Self-incrementing Unary operations */
 	template <class T> void FieldElement<T>::negate() throw() {
 		if (!parent_field) return;
 		parent_field->switchContext();
 		if (base) negate(repBase, repBase);
 		else negate(repExt, repExt);
 	}
 	
	template <class T> void FieldElement<T>::self_inv()
	throw(DivisionByZeroException) {
		if (isZero()) throw DivisionByZeroException();
		parent_field->switchContext();
		if (base) NTL::inv(repBase, repBase);
		else NTL::inv(repExt, repExt);
	}
	
	template <class T> void FieldElement<T>::operator^=(const ZZ& i)
	throw() {
		if (!parent_field) return;
		parent_field->switchContext();
		if (base) power(repBase, repBase, i);
		else power(repExt, repExt, i);
	}
	
	template <class T> void FieldElement<T>::operator^=(const long i)
	throw() {
		if (!parent_field) return;
		parent_field->switchContext();
		if (base) power(repBase, repBase, i);
		else power(repExt, repExt, i);
	}
	
	template <class T> void FieldElement<T>::self_trace() throw() {
		if (!parent_field) return;
		parent_field->switchContext();
		if (base) return;
		else {
			trace(repBase, repExt);
			base = true;
			repExt = 0;
			parent_field = &(parent_field->baseField());
		}
	}
	
/****************** Coercion of elements ******************/
/*	FieldElement<T> operator>>(const Field<T>&) const 
		throw(IllegalCoercionException);
	void operator>>=(const Field<T>&) throw(IllegalCoercionException);
	bool isCoercible(const Field<T>&) const throw();
*/	
/****************** Comparison ******************/
	template <class T> bool
	FieldElement<T>::operator==(const FieldElement<T>& e)
	const throw(NotInSameFieldException) {
		if (!parent_field) return e.isZero();
		if (!e.parent_field) return isZero();
		sameLevel(e);
		return base ? repBase == e.repBase : repExt == e.repExt;
	}
	template <class T> bool FieldElement<T>::operator==(const BigInt& i)
	const throw() {
		if (!parent_field) return i == long(0);
		else return base ? repBase == i : repExt == i;
	}

/****************** Infrastructure ******************/
	/* Interface with infrastructure. Use this only if you are
	 * sure of what you do !
	 */
	template <class T> void FieldElement<T>::toInfrastructure(GFp& i)
	const throw(IllegalCoercionException) {
		if (!parent_field || !base) throw IllegalCoercionException();
		parent_field->switchContext();
		i = repBase;
	}
	template <class T> void FieldElement<T>::toInfrastructure(GFpE& i)
	const throw(IllegalCoercionException) {
		if (!parent_field || !base) throw IllegalCoercionException();
		parent_field->switchContext();
		i = repExt;
	}

/****************** Printing ******************/
	template <class T> ostream& FieldElement<T>::print(ostream& o) const {
		if (!parent_field) return o << 0;
		if (base) return o << repBase;
		else return o << repExt;
	}
	
	/* Print the element as a polynomial over GF(p) in the
	 * variable var
	 */
//	ostream& print(ostream&, const string& var) const;
	/* Print the element as a multivariate polynomial over
	 * GF(p). The number of variables in vars must match one
	 * plus the Artin-Schreier height of the field the element
	 * belongs to.
	 * 
	 * throws : ASException if there's not enough variables
	 *          in var
	 */
//	ostream& print(ostream&, vector<const string>& vars) const;
}
