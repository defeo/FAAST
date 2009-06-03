#ifndef FIELDELEMENT_H_
#define FIELDELEMENT_H_

#include "Exceptions.hpp"
#include <string>
#include <vector>

namespace AS {
	template <class T> class Field;
	template <class T> class FieldElement;
	template <class T> class FieldPolynomial;

/****************** Level embedding ******************/
/* Find docs for these functions in the friends section of FieldElement */
	template <class T>
	void pushDown(const FieldElement<T>& e, vector<FieldElement<T> >& v)
	throw(NoSubFieldException);

	template <class T>
	void liftUp(const vector<FieldElement<T> >& v, FieldElement<T>& e)
	throw(NotInSameFieldException, NoOverFieldException);


/****************** Class FieldElement ******************/
	/**
	 * \ingroup Field
	 * \brief An element of a finite field.
	 * 
	 * Objects of this class represent elements of a finite field, as represented by the class Field.
	 * With the exception of the zero element created by the \link FieldElement() default constructor\endlink,
	 * any element belongs to one field and binary operations can combine two elements only in one of the
	 * following two cases:
	 *  - the two elements belong to the same field,
	 *  - one element belongs to a field and the other one to its \link Field::primeField() prime field\endlink.
	 * 
	 * Elements created through the \link FieldElement() default constructor\endlink, as for example
	 * \code
	 * FieldElement<T> elt;
	 * \endcode
	 * have a special status as they don't belong to any specific field: their default value is 0 and they
	 * can be combined with any other element. The result of a binary operation involving such a special
	 * element is one of the following:
	 *  - A DivisionByZeroException if the operation is division and the divisor is the special 0 element.
	 *  - An element of field \b K, if the other element belongs to \b K.
	 *  - The special 0 element if the other element is the special 0 element.
	 *  - An UndefinedFieldException if the other element is not the special 0 element and yet does
	 *    not belong to any field, as in this example
	 *    \code
	 *    FieldElement<T> elt;
	 *    return elt + 1;
	 *    \endcode
	 *    See UndefinedFieldException for more details.
	 * 
	 * Elements are internally represented as univariate polynomials with coefficients in F<sub>p</sub>
	 * modulo an irreducible polynomial as described in [\ref ISSAC "DFS '09", Section 3].
	 * The way the arithmetics of the field are actually implemented is
	 * given by the template parameter \a T that must be one of the \ref Infrastructures.
	 * Note that changing the Infrastructure may sensibly change the speed of your code.
	 *
	 * \tparam T An \ref Infrastructures "Infrastructure". It specfies which \NTL types will carry out
	 * the arithmetic operations.
	 * 
	 * \see Field, UndefinedFieldException
	 */
	template <class T> class FieldElement {

	friend class Field<T>;
	friend class FieldPolynomial<T>;
	/**
	 * \brief Convert \a e from the internal (univariate) representation to the bivariate representation
	 * over the immediate subfield in the primitive tower (the stem).
	 * 
	 * If \a e belongs to the base field of an Artin-Schreier tower, then \a v is filled with the
	 * coefficients in F<sub>p</sub> of its univariate representation. Otherwise let
	 * \code
	 * x = e.parent().primitiveElement();
	 * \endcode
	 * This method fills the vector \a v with \ref p elements of 
	 * \link Field::subField() \c e.parent()\c.subField() \endlink such that
	 * \f{equation}{
	 * \mathtt{e} = \mathtt{v[0]} + \mathtt{v[1]}*\mathtt{x} + ... 
	 * 		+ \mathtt{v[p-1]}*\mathtt{x}^{p-1}
	 * 		\mathrm{.}
	 * \f}
	 * 
	 * Let \b K be the field in the primitive tower (the stem) isomorphic to
	 * \link Field::subField() \c e.parent()\c.subField() \endlink, this corresponds to convert \a e from
	 * its internal (univariate) 
	 * representation to the bivariate representation as an element of \b K[\a x].
	 * This routine implements the algorithm \c PushDown of [\ref ISSAC "DFS '09", Section 4.2].
	 * 
	 * \param [in] e An element of any field.
	 * \param [out] v A vector of elements of \link Field::stemField() \c e.parent()\c.stemField() \endlink
	 * that satisfies condition (1).
	 * \throw NoSubFieldException If \a e belongs to F<sub>p</sub>.
	 * \note The result always lies in the primitive tower (the stem), even if \a e does not.
	 * \invariant If \a e belongs to a field in the primitive tower (the stem), this is equivalent to
	 * \code
	 * e.parent().subField().toBivariate(e, v);
	 * \endcode
	 * and this latter form should be preferred for the sake of clarity.
	 * \see Field::toBivariate().
	 *
	 * \relates FieldElement
	 */
	friend void pushDown<T>(const FieldElement<T>& e, vector<FieldElement<T> >& v) throw(NoSubFieldException);
	/**
	 * \brief Convert \a v from the bivariate representation 
	 * over the immediate subfield in the primitive tower (the stem) to the internal
	 * (univariate) representation.
	 * 
	 * If the elements of \a v belong to the prime field, then \a e is the element whose univariate
	 * representation has \a v as cofficients. Otherwise let
	 * \code
	 * x = v[0].parent().primitiveElement();
	 * \endcode
	 * This method stores in \a e an element of
	 * \link Field::overField() \c e.parent()\c.overField() \endlink such that
	 * \f{equation}{
	 * \mathtt{e} = \mathtt{v[0]} + \mathtt{v[1]}*\mathtt{x} + ... 
	 * 		+ \mathtt{v[p-1]}*\mathtt{x}^{p-1}
	 * 		\mathrm{,}
	 * \f}
	 * If \a v is too short, it is filled with zeros. If \a v is too
	 * long, the unnecessary elements are ignored.
	 * 
	 * Let \b K be the field in the primitive tower (the stem) isomorphic to
	 * the field containing the elements of \a v, this corresponds to convert \a v from the
	 * multivariate 
	 * representation as an element of \b K[\a x] to the internal (univariate) representation.
	 * This routine implements the algorithm \c LiftUp of [\ref ISSAC "DFS '09", Section 4.4].
	 *
	 * \param [in] v A vector of elements all belonging to the same field.
	 * \param [out] e An element satisfying condition (2).
	 * \throw NotInSameFieldException If the elements of \a v do not all
	 * belong to the same field.
	 * \throw NoOverFieldException If the field the elements of \a v belong to has no
	 * \link Field::overField() overfield\endlink.
	 * \note The result always lies in the primitive tower (the stem), even if the elements of \a v do not.
	 * \invariant If the elements of \a v belong to a field in the primitive tower (the stem),
	 * this is equivalent to
	 * \code
	 * v[0].parent().overField().toUnivariate(e, v);
	 * \endcode
	 * and this latter form should be preferred for the sake of clarity.
	 * \see Field::toUnivariate().
	 * 
	 * \relates FieldElement
	 */
	friend void liftUp<T>(const vector<FieldElement<T> >& v, FieldElement<T>& e) throw(NotInSameFieldException, NoOverFieldException);

	/** \name Local types
	 * Local types defined in this class. They are aliases to simplify the access
	 * to the \ref Infrastructures "Infrastructure" \a T and its subtypes.
	 *
	 * \see \ref Infrastructures.
	 * @{
	 */
	public:
		typedef T Infrastructure;

	private:
		typedef typename T::GFp     GFp;
		typedef typename T::MatGFp  MatGFp;
		typedef typename T::GFpX    GFpX;
		typedef typename T::GFpE    GFpE;
		typedef typename T::GFpEX   GFpEX;
		typedef typename T::BigInt  BigInt;
		typedef typename T::Context Context;
	/** @} */

	/****************** Members ******************/
	/** \cond DEV */
		/** \brief The \NTL representation of this element if it belongs to a prime field. */
		GFp repBase;
		/** \brief The \NTL representation of this element if it belongs to an extension field. */
		GFpE repExt;
		/** \brief Whether this element belongs to a prime or an extension field. */
		bool base;
		/** \brief The field this element belongs to. NULL if no parent field. */
		const Field<T>* parent_field;
	/** \endcond */

	public:
	/****************//** \name Constructors ******************/
	/** @{ */
		/**
		 * \brief Construct the special 0 element.
		 * 
		 * The special 0 element does not belong to any field, yet it can be added, multiplied, etc.
		 * to any other FieldElement. See the \link FieldElement introduction \endlink for more details.
		 * If you want to construct the 0 element of a specific field, use Field::zero() instead.
		 * 
		 * \see UndefinedFieldException, Field::zero().
		 */
		FieldElement() throw() : parent_field(NULL) {}
	/** @} */
	/****************//** \name Properties ******************/
	/** @{ */
		/**
		 * \brief The field this element belongs to.
		 * \throw UndefinedFieldException If this element does not belong to any field.
		 * \see FieldElement(), UndefinedFieldException.
		 */
		const Field<T>& parent() const
		throw(UndefinedFieldException) {
			if (!parent_field)
				throw UndefinedFieldException();
			return *parent_field;
		}
	/** @} */

	/****************//** \name Copy
	 * You can copy elements using assignment.
	 * @{
	 */
		FieldElement(const FieldElement<T>& e) throw();
		FieldElement<T>& operator=(const FieldElement<T>& e) throw();
		/**
		 * \brief Assign a scalar value to this element.
		 * 
		 * The result belongs to the same field as the element before the assignment.
		 * 
		 * \return A reference to the result. 
		 * \throw UndefinedFieldException If this is the \link FieldElement() special 0 element \endlink
		 * and \a i is different from 0.
		 */
		FieldElement<T>& operator=(const BigInt& i)
		throw(UndefinedFieldException);
	/** @} */

	/****************** Arithmetics ******************/
		/** \name Binary operators
		 * All binary operations throw a NotInSameFieldException if neither of this two conditions
		 * is satisfied:
		 *  - the two operands belong to the same field,
		 *  - one of the elements belongs to a field and the other to its prime field.
		 * @{
		 */
		FieldElement<T> operator+(const FieldElement<T>& e)
		const throw(NotInSameFieldException) {
			FieldElement<T> tmp = *this;
			tmp += e;
			return tmp;
		}
		FieldElement<T> operator-(const FieldElement<T>& e)
		const throw(NotInSameFieldException) {
			FieldElement<T> tmp = *this;
			tmp -= e;
			return tmp;
		}
		FieldElement<T> operator*(const FieldElement<T>& e)
		const throw(NotInSameFieldException) {
			FieldElement<T> tmp = *this;
			tmp *= e;
			return tmp;
		}
		FieldElement<T> operator/(const FieldElement<T>& e)
		const throw(NotInSameFieldException, DivisionByZeroException) {
			FieldElement<T> tmp = *this;
			tmp /= e;
			return tmp;
		}

		/* Self-incrementing binary operators. */
		void operator+=(const FieldElement<T>&)
			throw(NotInSameFieldException);
		void operator-=(const FieldElement<T>&)
			throw(NotInSameFieldException);
		void operator*=(const FieldElement<T>&)
			throw(NotInSameFieldException);
		void operator/=(const FieldElement<T>&)
			throw(NotInSameFieldException, DivisionByZeroException);

		/** \brief Stores \a a + \a b in this element. */
		void sum(const FieldElement<T>& a, const FieldElement<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator+=(b); }
		/** \brief Stores \a a - \a b in this element. */
		void difference(const FieldElement<T>& a, const FieldElement<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator-=(b); }
		/** \brief Stores \a a * \a b in this element. */
		void product(const FieldElement<T>& a, const FieldElement<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator*=(b); }
		/** \brief Stores \a a / \a b in this element. */
		void division(const FieldElement<T>& a, const FieldElement<T>& b)
			throw(NotInSameFieldException, DivisionByZeroException)
		{ operator=(a); operator/=(b); }
		/** @} */


		/** \name Unary operators */
		/** @{ */
		/** \brief Additive inverse. */
		FieldElement<T> operator-() const throw() {
			FieldElement<T> tmp = *this;
			tmp.negate();
			return tmp;
		}
		/** \brief Multiplicative inverse. */
		FieldElement<T> inv() const throw(DivisionByZeroException) {
			FieldElement<T> tmp = *this;
			tmp.self_inv();
			return tmp;
		}
		/** \brief Power. */
		FieldElement<T> operator^(const ZZ& i) const throw() {
			FieldElement<T> tmp = *this;
			tmp ^= i;
			return tmp;
		}
		/** \brief \copybrief operator^() */
		FieldElement<T> operator^(const long i) const throw() {
			FieldElement<T> tmp = *this;
			tmp ^= i;
			return tmp;
		}
		/* Frobenius and iterated frobenius */
		/** \brief <i>p</i>-th power (frobenius morphism). */
		FieldElement<T> frobenius() const throw() {
			FieldElement<T> tmp = *this;
			tmp.self_frobenius();
			return tmp;
		}
		/** \brief <i>p<sup>n</sup></i>-th power (iterated frobenius morphism).
		 * 
		 * This method relies on the algorithm \c IterFrobenius of [\ref ISSAC "DFS '09", Section 5].
		 */
		FieldElement<T> frobenius(const long n) const throw() {
			FieldElement<T> tmp = *this;
			tmp.self_frobenius(n);
			return tmp;
		}
		/**
		 * \brief Trace over the field \a F.
		 *
		 * \throws NotASubFieldException If this element does not belong to
		 *          an overfield of \a F.
		 */
		FieldElement<T> trace(const Field<T>& F)
		const throw(NotASubFieldException) {
			FieldElement<T> tmp = *this;
			tmp.self_trace(F);
			return tmp;
		}
		/** \brief Trace over F<sub>p</sub> */
		FieldElement<T> trace() const throw();
		/**
		 * \brief <i>n</i>-th pseudotrace
		 * 
		 * Let \a x be this element, the <i>n</i>-th pseudotrace, noted
		 * T<sub>n</sub>(x) is
		 * \f[
		 * 	\mathrm{T}_n(x) = \sum_{\ell=0}^{n-1} x^{p^\ell}
		 * 	\mathrm{.}
		 * \f]
		 * This method relies on the algorithm \c Pseudotrace of [\ref ISSAC "DFS '09", Section 5].
		 */
		FieldElement<T> pseudotrace(unsigned long n) const throw() {
			FieldElement<T> tmp = *this;
			tmp.self_pseudotrace(n);
			return tmp;
		}
		/** @} */

		/** \name Self-incrementing Unary operations
		 * These operators store the result of the operation into the element itself.
		 * @{
		 */
		/** \brief Flip the sign of this element. */
	 	void negate() throw();
	 	/** \brief Invert this element. */
		void self_inv() throw(DivisionByZeroException);
		/** \brief \copybrief operator^() */
		void operator^=(const ZZ&) throw();
		/** \brief \copybrief operator^() */
		void operator^=(const long) throw();
		/** \brief \copybrief frobenius() const */
		void self_frobenius() throw();
		/** \brief \copybrief frobenius(const long) const */
		void self_frobenius(long) throw();
		/** \brief \copybrief trace(const Field<T>&) const */
		void self_trace(const Field<T>& F) throw(NotASubFieldException);
		/** \brief \copybrief trace() const */
		void self_trace() throw();
		/** \brief \copybrief pseudotrace()
		 * 
		 * \copydetails pseudotrace()
		 */
		void self_pseudotrace(unsigned long) throw();
		/** @} */

	/****************//** \name Minimal polynomials and Evaluation******************/
	/** @{ */
		/* The minimal polynomial over the base field. */
		FieldPolynomial<T> minimalPolynomial() const throw() {
			parent_field->switchContext();
			GFpX minpol;
			MinPolyMod(minpol, rep(repExt), GFpE::modulus());
			return parent_field->primeField().fromInfrastructure(minpol);
		}
		/* The minimal polynomial over the field F
		 *
		 * throws: NotASubFieldException if F is not a subfield of this.parent
		 * throws: NotSupportedException if F is not part of an Artin-Schreier tower
		 */
		FieldPolynomial<T> minimalPolynomial(const Field<T>& F)
		const throw(NotASubFieldException, NotSupportedException) {
			vector<FieldPolynomial<T> > minpols;
			minimalPolynomials(F, minpols);
			return minpols[0];
		}
		/* All the minimal polynomials up to the field F. If
		 * called this way :
		 *		minimalPolynomials(F, res);
		 * res[0] will contain minimalPolynomials(F), res[1]
		 * will contain minimalPolynomials(F.overfield) and so on.
		 *
		 * throws: NotASubFieldException if F is not a subfield of this.parent
		 * throws: NotSupportedException if F is not part of an Artin-Schreier tower
		 */
		void minimalPolynomials(const Field<T>& F, vector<FieldPolynomial<T> >&)
		const throw(NotASubFieldException, NotSupportedException);

		/* The a-affine minimal polynomial over the field F,
		 * That is the minimum degree polynomial P of F[X] such that
		 * 		P(this) = a.
		 * The optional parameter minpols must contain either the result of
		 * this.minimalPolynomials(F,v), or must be an empty vector, in which
		 * case it is filled with the result of minimalPolynomials(F,v).
		 *
		 * throws: NotASubFieldException if F is not a subfield of this.parent
		 * throws: NoSuchPolynomialException if a is not in the field
		 * 			generated by this
		 * throws: NotSupportedException if F is not part of an Artin-Schreier tower
		 */
		FieldPolynomial<T> affineMinimalPolynomial(
		const Field<T>& F, const FieldElement<T>& a,
		vector<FieldPolynomial<T> >& minpols)
		const throw(NotASubFieldException, NoSuchPolynomialException,
		NotSupportedException, BadParametersException);

		FieldPolynomial<T> affineMinimalPolynomial(
		const Field<T>& F, const FieldElement<T>& a)
		const throw(NotASubFieldException, NoSuchPolynomialException,
		NotSupportedException, BadParametersException) {
			vector<FieldPolynomial<T> > minpols;
			return affineMinimalPolynomial(F,a,minpols);
		}

		/* Returns P(this).
		 *
		 * The optional parameter minpols must contain either the result of
		 * this.minimalPolynomials(P.parent,v), or must be an empty vector,
		 * in which case it is filled with the result of
		 * minimalPolynomials(P.parent,v).
		 *
		 * throws : IllegalCoercionException if this cannot be coerced to P.parent
		 * 			nor can P be coerced to this.parent
		 * throws : BadParametersException if minpols contains bad data
		 */
		FieldElement<T> evaluate(const FieldPolynomial<T>& P,
		vector<FieldPolynomial<T> >& minpols)
		const throw(IllegalCoercionException, BadParametersException);

		FieldElement<T> evaluate(const FieldPolynomial<T>& P)
		const throw(IllegalCoercionException, BadParametersException) {
			vector<FieldPolynomial<T> > minpols;
			return evaluate(P, minpols);
		}
	/** @} */

	/****************** Coercion of elements ******************/
		FieldElement<T> toScalar() const throw(IllegalCoercionException);
		FieldElement<T> operator>>(const Field<T>& F) const
		throw(IllegalCoercionException) {
			FieldElement<T> tmp = *this;
			tmp >>= F;
			return tmp;
		}
		void operator>>=(const Field<T>&) throw(IllegalCoercionException);
		bool isCoercible(const Field<T>&) const throw();

	/****************** Comparison ******************/
		bool operator==(const FieldElement<T>&) const throw(NotInSameFieldException);
		bool operator==(const BigInt&) const throw();
		bool operator!=(const FieldElement<T>& e) const throw(NotInSameFieldException)
		{ return !(*this==e); }
		bool operator!=(const BigInt& i) const throw() { return !(*this==i); }
		bool isZero() const throw() {
			return !parent_field ||
				(base ? IsZero(repBase) : IsZero(repExt));
		}
		bool isOne() const throw() {
			return parent_field &&
				(base ? IsOne(repBase) : IsOne(repExt));
		}
		bool isScalar() const throw();
	/****************** Infrastructure ******************/
		/* Interface with infrastructure. Use this only if you are
		 * sure of what you do !
		 */
		void toInfrastructure(GFp&) const throw(IllegalCoercionException);
		void toInfrastructure(GFpE&) const throw(IllegalCoercionException);

	/****************** Printing ******************/
		ostream& print(ostream&) const;
		/* Print the element as a polynomial over GF(p) in the
		 * variable var
		 */
		ostream& print(ostream&, const string& var) const;
		/* Print the element as a multivariate polynomial over
		 * GF(p). The number of variables in vars must be at least
		 * one plus the Artin-Schreier height of the field the
		 * element belongs to.
		 *
		 * throws : ASException if there's not enough variables
		 *          in var
		 */
		ostream& print(ostream&, const vector<string>& vars) const;
	/****************** Destructor ******************/
		~FieldElement() throw() {}


	/*****************************************************/
	/****************** Private section ******************/
	/*****************************************************/

	private:
	/****************** Helpers for frobenius and trace ******************/
		/* p^jd-th iterated frobenius */
		void BigFrob(const long j);
		/* n-th iterated frobenius, n < d */
		void SmallFrob(const long n);
		/* p^jd-th pseudotrace */
		void BigPTrace(const long j);
		/* Put in the vector v all the p^id pseudotraces for 0 <= i <= j */
		void BigPTraceVector(vector<FieldElement<T> >& v, const long j) const;
		/* n-th pseudotrace, n < d */
		void SmallPTrace(const long n);

	/****************** Internal Constructors ******************/
		/* Construct an element with given representation and parent.
		 * Reserved for used by Field<T>
		 */
		FieldElement(const Field<T>* p, const GFp& PBase, const GFpE& PExt, const bool b) throw() :
			repBase(PBase), repExt(PExt), base(b), parent_field(p) {}
		FieldElement(const Field<T>* p, const GFpE& P) throw() :
			repExt(P), base(false), parent_field(p) {}
		FieldElement(const Field<T>* p, const GFp& P) throw() :
			repBase(P), base(true), parent_field(p) {}
	/****************** Utility Routines ******************/
		void sameLevel(const FieldElement<T>& e) const
		throw(NotInSameFieldException) {
			if (parent_field != e.parent_field)
			throw NotInSameFieldException();
		}
	};

/****************** Printing ******************/
	template <class T> ostream&
	operator<<(ostream& o, const FieldElement<T>& e) {
		return e.print(o);
	}
}

#endif /*FIELDELEMENT_H_*/
