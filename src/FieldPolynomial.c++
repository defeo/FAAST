namespace AS {
/****************** Properties ******************/
	/* Degree. Returns -1 if the polynomial is 0. */
	template <class T> long FieldPolynomial<T>::degree()
	const throw() {
		if (!parent_field) return -1;
		if (base) return deg(repBase);
		else return deg(repExt);
	}

/****************** Copy ******************/
	template <class T>
	FieldPolynomial<T>::FieldPolynomial(const FieldPolynomial<T>& e)
	throw() : repBase(), repExt(), base(e.base),
	parent_field(e.parent_field) {
		if (parent_field) {
			parent_field->switchContext();
			if (base) repBase = e.repBase;
			else repExt = e.repExt;
		}
	}

	template <class T> FieldPolynomial<T>&
	FieldPolynomial<T>::operator=(const FieldPolynomial<T>& e)
	throw() {
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

	/* Starting only from the constant coefficient */
	template <class T>
	FieldPolynomial<T>::FieldPolynomial(const FieldElement<T>& e)
	throw() : repBase(), repExt(), base(e.base),
	parent_field(e.parent_field) {
		if (parent_field) {
			parent_field->switchContext();
			if (base) SetCoeff(repBase, 0, e.repBase);
			else SetCoeff(repExt, 0, e.repExt);
		}
	}

	template <class T> FieldPolynomial<T>&
	FieldPolynomial<T>::operator=(const FieldElement<T>& e)
	throw() {
		repBase = 0;
		repExt = 0;
		base = e.base;
		parent_field = e.parent_field;
		if (parent_field) {
			parent_field->switchContext();
			if (base) SetCoeff(repBase, 0, e.repBase);
			else SetCoeff(repExt, 0, e.repExt);
		}
		return *this;
	}

	template <class T> FieldPolynomial<T>&
	FieldPolynomial<T>::operator=(const BigInt& i)
	throw(UndefinedFieldException) {
		if (!parent_field) {
			if (i == long(0)) return *this;
			else throw UndefinedFieldException();
		}
		parent_field->switchContext();
		repBase = 0;
		repExt = 0;
		if (base) SetCoeff(repBase, 0, i);
		else SetCoeff(repExt, i);
		return *this;
	}

/****************** Coefficients ******************/
	template <class T> void
	FieldPolynomial<T>::getCoeff(const long i, FieldElement<T>& e)
	const throw(BadParametersException) {
		if (i < 0)
			throw BadParametersException("Negative index for polynomial coefficient.");
		if (!parent_field) {
			e = FieldElement<T>();
			return;
		}
		parent_field->switchContext();
		e.parent_field = parent_field;
		e.base = base;
		e.repBase = 0;
		e.repExt = 0;
		if (base) e.repBase = coeff(repBase, i);
		else e.repExt = coeff(repExt, i);
	}

	template <class T> void
	FieldPolynomial<T>::setCoeff(const long i, const FieldElement<T>& e)
	throw(NotInSameFieldException, BadParametersException) {
		if (i < 0)
			throw BadParametersException("Negative index for polynomial coefficient.");
		if (e.isZero()) return;
		if (!parent_field) {
			parent_field = e.parent_field;
			base = e.base;
		} else sameLevel(e);

		parent_field->switchContext();
		if (base) SetCoeff(repBase, i, e.repBase);
		else SetCoeff(repExt, i, e.repExt);
	}

	template <class T> void
	FieldPolynomial<T>::setCoeff(const long i, const BigInt& c)
	throw(UndefinedFieldException, BadParametersException) {
		if (i < 0)
			throw BadParametersException("Negative index for polynomial coefficient.");
		if (c == 0) return;
		if (!parent_field)
			throw UndefinedFieldException();

		parent_field->switchContext();
		if (base) SetCoeff(repBase, i, c);
		else SetCoeff(repExt, i, c);
	}

	template <class T> void
	FieldPolynomial<T>::setCoeff(const long i)
	throw(UndefinedFieldException, BadParametersException) {
		if (i < 0)
			throw BadParametersException("Negative index for polynomial coefficient.");
		if (!parent_field)
			throw UndefinedFieldException();

		parent_field->switchContext();
		if (base) SetCoeff(repBase, i);
		else SetCoeff(repExt, i);
	}

/****************** Arithmetics ******************/
	/* Self-incrementing binary operations. */
	template <class T> void
	FieldPolynomial<T>::operator+=(const FieldPolynomial<T>& e)
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
	FieldPolynomial<T>::operator-=(const FieldPolynomial<T>& e)
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
	FieldPolynomial<T>::operator*=(const FieldPolynomial<T>& e)
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
	FieldPolynomial<T>::operator/=(const FieldPolynomial<T>& e)
	throw(NotInSameFieldException, DivisionByZeroException) {
		if (e.isZero()) throw DivisionByZeroException();
		if (!parent_field) {
			*this = e.parent_field->zero();
			return;
		}
		sameLevel(e);
		parent_field->switchContext();
		if (base) repBase /= e.repBase;
		else repExt /= e.repExt;
	}

	template <class T> void
	FieldPolynomial<T>::operator%=(const FieldPolynomial<T>& e)
	throw(NotInSameFieldException, DivisionByZeroException) {
		if (e.isZero()) throw DivisionByZeroException();
		if (!parent_field) {
			*this = e.parent_field->zero();
			return;
		}
		sameLevel(e);
		parent_field->switchContext();
		if (base) repBase %= e.repBase;
		else repExt %= e.repExt;
	}

	/* Unary operations */
	template <class T> bool
	FieldPolynomial<T>::divides(const FieldPolynomial<T>& e)
	const throw() {
		if (isZero()) return false;
		if (e.isZero()) return true;
		sameLevel(e);
		parent_field->switchContext();
		if (base) return divide(e.repBase, repBase) == 1;
		else return divide(e.repExt, repExt) == 1;
	}

	/* Self-incrementing Unary operations */
	template <class T> void FieldPolynomial<T>::negate() throw() {
 		if (!parent_field) return;
 		parent_field->switchContext();
 		if (base) NTL::negate(repBase, repBase);
 		else NTL::negate(repExt, repExt);
	}

	template <class T> void
	FieldPolynomial<T>::operator^=(const long i) throw() {
		if (!parent_field) return;
		parent_field->switchContext();
		if (base) power(repBase, repBase, i);
		else power(repExt, repExt, i);
	}

	template <class T> void
	FieldPolynomial<T>::self_derivative() throw() {
		if (!parent_field) return;
		parent_field->switchContext();
		if (base) diff(repBase, repBase);
		else diff(repExt, repExt);
	}

	template <class T> void
	FieldPolynomial<T>::normalize() throw() {
		if (!parent_field) return;
		parent_field->switchContext();
		if (base) MakeMonic(repBase);
		else MakeMonic(repExt);
	}

	template <class T> void
	FieldPolynomial<T>::self_frobenius() throw() {
		if (isZero()) return;
		if (base) return;

		parent_field->switchContext();
		GFpE tmp;
		for (long i = deg(repExt) ; i >= 0 ; i--) {
			power(tmp, coeff(repExt, i), parent_field->p);
			SetCoeff(repExt, i, tmp);
		}
	}

	template <class T> void
	FieldPolynomial<T>::self_frobenius(long n) throw() {
		if (isZero()) return;
		if (base) return;
		n %= parent_field->d;
		if (n < 0) n = parent_field->d - n;
		if (n == 0) return;

		FieldElement<T> tmp;
		for (long i = degree() ; i >= 0 ; i--) {
			getCoeff(i, tmp);
			tmp.self_frobenius(n);
			setCoeff(i, tmp);
		}
	}

/****************** Coercion of elements ******************/
	template <class T> FieldPolynomial<T>
	FieldPolynomial<T>::toScalarPolynomial()
	const throw(IllegalCoercionException) {
		if (!parent_field) return *this;
		if (parent_field->d == 1)
			return *this >> parent_field->primeField();

#ifdef AS_DEBUG
		if (base) throw
			ASException("Malformed element in isScalarPolynomial().");
#endif

		parent_field->switchContext();
		if (isScalarPolynomial()) {
			GFpX e; conv(e, repExt);
			return FieldPolynomial<T>(
					&(parent_field->primeField()), e);
		} else throw IllegalCoercionException();
	}

	template <class T> FieldPolynomial<T>
	FieldPolynomial<T>::operator>>(const Field<T>& F)
	const throw(IllegalCoercionException) {
		if (!parent_field) {
			return F.zero();
		} else if (parent_field->isIsomorphic(F)) {
			FieldPolynomial<T> res = *this;
			res.parent_field = &F;
			return res;
		} else if (parent_field->isSubFieldOf(F) ||
			parent_field->isOverFieldOf(F)) {
			FieldPolynomial<T> res = F.zero();
			FieldElement<T> tmp;
			for (long i = degree() ; i >= 0 ; i--) {
				getCoeff(i, tmp);
				tmp >>= F;
				res.setCoeff(i, tmp);
			}
			return res;
		} else throw IllegalCoercionException();
	}

	template <class T> bool
	FieldPolynomial<T>::isCoercible(const Field<T>& F)
	const throw() {
		if (!parent_field) return true;

		// go up ...
		if (parent_field->isSubFieldOf(F)) return true;
		// or go down ...
		else if (parent_field->isOverFieldOf(F)) {
			FieldElement<T> tmp;
			for (long i = degree() ; i  >0 ; i--) {
				getCoeff(i, tmp);
				if (!tmp.isCoercible(F)) return false;
			}
			return true;
		}
	}

/****************** Comparison ******************/
	template <class T> bool
	FieldPolynomial<T>::operator==(const FieldPolynomial<T>& e)
	const throw(NotInSameFieldException) {
		if (!parent_field) return e.isZero();
		if (!e.parent_field) return isZero();
		sameLevel(e);
		return base ? repBase == e.repBase : repExt == e.repExt;
	}

	template <class T> bool
	FieldPolynomial<T>::operator==(const FieldElement<T>& e)
	const throw(NotInSameFieldException) {
		if (!parent_field) return e.isZero();
		if (!e.parent_field) return isZero();

		parent_field->switchContext();
		// if e is a scalar
		if (parent_field->primeField().isIsomorphic(e.parent()))
			return base ? repBase == e.repBase : repExt == e.repBase;
		// generic case
		sameLevel(e);
		return repExt == e.repExt;
	}

	template <class T> bool
	FieldPolynomial<T>::operator==(const BigInt& i) const throw() {
		if (!parent_field) return i == long(0);
		else {
			parent_field->switchContext();
			return base ? repBase == i : repExt == i;
		}
	}

	template <class T> bool
	FieldPolynomial<T>::isScalarPolynomial() const throw() {
		if (!parent_field || parent_field->d == 1) return true;

#ifdef AS_DEBUG
		if (base) throw
			ASException("Malformed element in isScalarPolynomial().");
#endif

		parent_field->switchContext();
		for (long i = degree() ; i >= 0 ; i--) {
			if (deg(rep(coeff(repExt, i))) > 0) return false;
		}
		return true;
	}

/****************** Infrastructure ******************/
	/* Interface with infrastructure. Use this only if you are
	 * sure of what you do !
	 */
	template <class T> void
	FieldPolynomial<T>::toInfrastructure(GFpX& p)
	const throw(IllegalCoercionException) {
		if (!parent_field || !base) throw IllegalCoercionException();
		parent_field->switchContext();
		p = repBase;
	}

	template <class T> void
	FieldPolynomial<T>::toInfrastructure(GFpEX& p)
	const throw(IllegalCoercionException) {
		if (!parent_field || !base) throw IllegalCoercionException();
		parent_field->switchContext();
		p = repExt;
	}

/****************** Printing ******************/
	template <class T> ostream&
	FieldPolynomial<T>::print(ostream& o) const {
		if (!parent_field) return o << 0;
		if (base) return o << repBase;
		else return o << repExt;
	}

	/* Print the element as a polynomial over the base field in the
	 * variable varPoly. varField is used to print the element of the field.
	 */
	template <class T> ostream&
	FieldPolynomial<T>::print(ostream& o, const string& varPoly,
	const string& varField) const {
		if (isZero()) return o << 0;
		bool first = true;
		FieldElement<T> c;
		for (long i = degree() ; i >= 0 ; i--) {
			getCoeff(i, c);
			if (c != 0) {
				if (first) first = false;
				else o << " + ";
				if (!c.isScalar()) o << "(";
				if (c != 1) c.print(o, varField);
				if (!c.isScalar()) o << ")";
				if (c != 1 && i != 0) o << "*";
				if (i != 0) o << varPoly;
				if (i > 1) o << "^" << i;
			}
		}
		if (first) o << 0;
		return o;
	}

	/* Print the element as a polynomial over the base field in the
	 * variable varPoly.
	 * varsField is used to print the element of the field on the
	 * multivariate basis. The number of variables in vars must match
	 * one plus the Artin-Schreier height of the base field.
	 *
	 * throws : ASException if there's not enough variables
	 *          in var
	 */
	template <class T> ostream&
	FieldPolynomial<T>::print(ostream& o, const string& varPoly,
	const vector<string>& varsField) const {
		if (isZero()) return o << 0;
		bool first = true;
		FieldElement<T> c;
		for (long i = degree() ; i >= 0 ; i--) {
			getCoeff(i, c);
			if (c != 0) {
				if (first) first = false;
				else o << " + ";
				if (!c.isScalar()) o << "(";
				if (c != 1) c.print(o, varsField);
				if (!c.isScalar()) o << ")";
				if (c != 1 && i != 0) o << "*";
				if (i != 0) o << varPoly;
				if (i > 1) o << "^" << i;
			}
		}
		if (first) o << 0;
		return o;
	}
}
