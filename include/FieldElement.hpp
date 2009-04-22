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
	/* Push the element e down along the stem and store
	 * the result in v.
	 *
	 * throw : NoSubFieldException if e belongs to GF(p)
	 */
	template <class T>
	void pushDown(const FieldElement<T>& e, vector<FieldElement<T> >& v)
	throw(NoSubFieldException);

	/* Lift the elements in v up along the stem and store
	 * the result in e.
	 * If v is too short, it is filled with zeros.
	 * If v is too long, the unnecessary elements are ignored.
	 *
	 * throw : NotInSameFieldException if the elements of v do not
	 *         belong all to the same field.
	 * throw : NoOverFieldException if there's no extension to lift
	 *         up to.
	 */
	template <class T>
	void liftUp(const vector<FieldElement<T> >& v, FieldElement<T>& e)
	throw(NotInSameFieldException, NoOverFieldException);


/****************** Class FieldElement ******************/
	template <class T> class FieldElement {

	friend class Field<T>;
	friend class FieldPolynomial<T>;
	friend void pushDown<T>(const FieldElement<T>& e, vector<FieldElement<T> >& v) throw(NoSubFieldException);
	friend void liftUp<T>(const vector<FieldElement<T> >& v, FieldElement<T>& e)	throw(NotInSameFieldException, NoOverFieldException);

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


	/****************** Members ******************/
		/* The representation of this element */
		GFp repBase;
		GFpE repExt;
		bool base;
		/* The field this element belongs to */
		const Field<T>* parent_field;


	public:
	/****************** Constructors ******************/
		/* Constructor by default.
		 * Constructs the 0 element (of any field).
		 */
		FieldElement() throw() : parent_field(NULL) {}
	/****************** Properties ******************/
		/* The field this element belongs to */
		const Field<T>& parent() const
		throw(UndefinedFieldException) {
			if (!parent_field)
				throw UndefinedFieldException();
			return *parent_field;
		}

	/****************** Copy ******************/
		FieldElement(const FieldElement<T>& e) throw();
		FieldElement<T>& operator=(const FieldElement<T>& e) throw();
		FieldElement<T>& operator=(const BigInt& i)
		throw(UndefinedFieldException);

	/****************** Arithmetics ******************/
		/* Binary operations */
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

		/* Self-incrementing binary operations. */
		void operator+=(const FieldElement<T>&)
			throw(NotInSameFieldException);
		void operator-=(const FieldElement<T>&)
			throw(NotInSameFieldException);
		void operator*=(const FieldElement<T>&)
			throw(NotInSameFieldException);
		void operator/=(const FieldElement<T>&)
			throw(NotInSameFieldException, DivisionByZeroException);

		/* Memory-efficient (NTL-like) binary operations.
		 * Aplly on the arguments and store in this.
		 */
		void sum(const FieldElement<T>& a, const FieldElement<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator+=(b); }
		void difference(const FieldElement<T>& a, const FieldElement<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator-=(b); }
		void product(const FieldElement<T>& a, const FieldElement<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator*=(b); }
		void division(const FieldElement<T>& a, const FieldElement<T>& b)
			throw(NotInSameFieldException, DivisionByZeroException)
		{ operator=(a); operator/=(b); }

		/* Unary operations */
		FieldElement<T> operator-() const throw() {
			FieldElement<T> tmp = *this;
			tmp.negate();
			return tmp;
		}
		FieldElement<T> inv() const throw(DivisionByZeroException) {
			FieldElement<T> tmp = *this;
			tmp.self_inv();
			return tmp;
		}
		FieldElement<T> operator^(const ZZ& i) const throw() {
			FieldElement<T> tmp = *this;
			tmp ^= i;
			return tmp;
		}
		FieldElement<T> operator^(const long i) const throw() {
			FieldElement<T> tmp = *this;
			tmp ^= i;
			return tmp;
		}
		/* Frobenius and iterated frobenius */
		FieldElement<T> frobenius() const throw() {
			FieldElement<T> tmp = *this;
			tmp.self_frobenius();
			return tmp;
		}
		FieldElement<T> frobenius(const long n) const throw() {
			FieldElement<T> tmp = *this;
			tmp.self_frobenius(n);
			return tmp;
		}
		/* Trace over the field F.
		 *
		 * throws : NotASubFieldException if this does not belong to
		 *          an overfield of F.
		 */
		FieldElement<T> trace(const Field<T>& F)
		const throw(NotASubFieldException) {
			FieldElement<T> tmp = *this;
			tmp.self_trace(F);
			return tmp;
		}
		/* Absolute trace over GF(p) */
		FieldElement<T> trace() const throw();
		/* n-th pseudotrace */
		FieldElement<T> pseudotrace(unsigned long n) const throw() {
			FieldElement<T> tmp = *this;
			tmp.self_pseudotrace(n);
			return tmp;
		}

		/* Self-incrementing Unary operations */
	 	void negate() throw();
		void self_inv() throw(DivisionByZeroException);
		void operator^=(const ZZ&) throw();
		void operator^=(const long) throw();
		void self_frobenius() throw();
		void self_frobenius(long) throw();
		void self_trace(const Field<T>& F) throw(NotASubFieldException);
		void self_trace() throw();
		void self_pseudotrace(unsigned long) throw();

	/****************** Minimal polynomials and Evaluation******************/
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

#include "../src/FieldElement.c++"
#include "../src/FE-Liftup-Pushdown.c++"
#include "../src/FE-Trace-Frob.c++"
#include "../src/Minpols.c++"

#endif /*FIELDELEMENT_H_*/
