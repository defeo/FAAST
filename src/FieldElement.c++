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
			if (i == long(0)) return *this;
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
		if (base && !e.base &&
			&(e.parent_field->primeField()) == parent_field) {
			e.parent_field->switchContext();
			repExt = e.repExt;
			repExt += repBase;
			repBase = 0;
			base = false;
			parent_field = e.parent_field;
			return;
		}
		if (!base && e.base &&
			&(parent_field->primeField()) == e.parent_field) {
			parent_field->switchContext();
			repExt += e.repBase;
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
		if (base && !e.base &&
			&(e.parent_field->primeField()) == parent_field) {
			e.parent_field->switchContext();
			repExt = e.repExt;
			repExt -= repBase;
			repBase = 0;
			base = false;
			parent_field = e.parent_field;
			return;
		}
		if (!base && e.base &&
			&(parent_field->primeField()) == e.parent_field) {
			parent_field->switchContext();
			repExt -= e.repBase;
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
		if (base && !e.base &&
			&(e.parent_field->primeField()) == parent_field) {
			e.parent_field->switchContext();
			repExt = e.repExt;
			repExt *= repBase;
			repBase = 0;
			base = false;
			parent_field = e.parent_field;
			return;
		}
		if (!base && e.base &&
			&(parent_field->primeField()) == e.parent_field) {
			parent_field->switchContext();
			repExt *= e.repBase;
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
		if (base && !e.base &&
			&(e.parent_field->primeField()) == parent_field) {
			e.parent_field->switchContext();
			repExt = e.repExt;
			repExt /= repBase;
			repBase = 0;
			base = false;
			parent_field = e.parent_field;
			return;
		}
		if (!base && e.base &&
			&(parent_field->primeField()) == e.parent_field) {
			parent_field->switchContext();
			repExt /= e.repBase;
			return;
		}
		
		sameLevel(e);
		parent_field->switchContext();
		if (base) repBase /= e.repBase;
		else repExt /= e.repExt;
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
			tmp.parent_field = &(parent_field->primeField());
			return tmp;
		}
	}

	/* Self-incrementing Unary operations */
 	template <class T> void FieldElement<T>::negate() throw() {
 		if (!parent_field) return;
 		parent_field->switchContext();
 		if (base) NTL::negate(repBase, repBase);
 		else NTL::negate(repExt, repExt);
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
	
	template <class T> void FieldElement<T>::self_frobenius()
	throw() {
		if (!parent_field) return;
		if (base) return;
		
		parent_field->switchContext();
		
		power(repExt, repExt, parent_field->p);
	}
	
	template <class T> void FieldElement<T>::self_trace() throw() {
		if (!parent_field) return;
		parent_field->switchContext();
		if (base) return;
		else {
			NTL::trace(repBase, repExt);
			base = true;
			repExt = 0;
			parent_field = &(parent_field->baseField());
		}
	}
	
/****************** Coercion of elements ******************/
	template <class T> FieldElement<T>
	FieldElement<T>::toScalar()
	const throw(IllegalCoercionException) {
		if (!parent_field) return *this;
		if (parent_field->d == 1)
			return *this >> parent_field->primeField();
			
#ifdef AS_DEBUG
		if (base) throw
			ASException("Malformed element in isScalar().");
#endif

		parent_field->switchContext();
		if (deg(rep(repExt)) <= 0) {
			GFp e; conv(e, repExt);
			return FieldElement<T>(
					&(parent_field->primeField()), e);
		} else throw IllegalCoercionException();
	}
	
	template <class T> void
	FieldElement<T>::operator>>=(const Field<T>& F)
	throw(IllegalCoercionException) {
		if (!parent_field) {
			*this = F.zero();
			return;
		}
		
		// go up ...
		if (parent_field->isSubFieldOf(F)) {
			while (parent_field->stem != F.stem) {
				vector<FieldElement<T> > v;
				v.resize(1);
				v[0] = *this;
				AS::liftUp(v, *this);
			}
		}
		// or go down ...
		else if (parent_field->isOverFieldOf(F)) {
			FieldElement<T> bak = *this;
			while (parent_field->stem != F.stem) {
				vector<FieldElement<T> > v;
				AS::pushDown(*this, v);
				// verify coercibility
				typename vector<FieldElement<T> >::iterator it = v.begin();
				if (it != v.end()) {
					*this = *it;
					for (it++ ; it != v.end() ; it++) {
						if (*it != 0) {
							*this = bak;
							throw IllegalCoercionException();
						}
					}
				} else *this = F.zero();
			}
		}
		// or go nowhere
		else throw IllegalCoercionException();
		
		// move out of the stem, if needed
		parent_field = &F;
	}
	 
	template <class T> bool
	FieldElement<T>::isCoercible(const Field<T>& F)
	const throw() {
		if (!parent_field) return true;
		
		// go up ...
		if (parent_field->isSubFieldOf(F)) return true;
		// or go down ...
		else if (parent_field->isOverFieldOf(F)) {
			while (parent_field.stem != F.stem) {
				vector<FieldElement<T> > v;
				AS::pushDown(*this, v);
				// verify coercibility
				typename vector<FieldElement<T> >::iterator it = v.begin();
				if (it == v.end()) return true;
				else {
					for (it++ ; it != v.end() ; it++) {
						if (*it != 0) return false;
					}
				}
			}
			return true;
		}
		// or go nowhere
		else return false;
	}
	
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
	template <class T> bool FieldElement<T>::isScalar() const
	throw() {
		if (!parent_field  || parent_field->d == 1) return true;
#ifdef AS_DEBUG
		if (base) throw
			ASException("Malformed element in isScalar().");
#endif
		parent_field->switchContext();
		return deg(rep(repExt)) <= 0;
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
		if (!parent_field || base) throw IllegalCoercionException();
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
