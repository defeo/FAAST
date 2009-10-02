/*
	This file is part of the FAAST library.

	Copyright (c) 2009 Luca De Feo and Éric Schost.

	The most recent version of FAAST is available at http://www.lix.polytechnique.fr/~defeo/FAAST

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; see file COPYING. If not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
#ifndef FIELDPOLYNOMIAL_H_
#define FIELDPOLYNOMIAL_H_

#include "Exceptions.hpp"
#include <string>
#include <vector>

namespace FAAST {
	template <class T> class Field;

/****************** GCD ******************/
/* Find docs for these functions in the friends section of FieldElement */
	template <class T> FieldPolynomial<T>
	GCD(const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q)
	throw(NotInSameFieldException);

	template <class T> void
	XGCD(const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q,
			FieldPolynomial<T>& U, FieldPolynomial<T>& V, FieldPolynomial<T>& G)
	throw(NotInSameFieldException);

	template <class T> void
	HalfGCD(FieldPolynomial<T>& U0, FieldPolynomial<T>& V0,
			FieldPolynomial<T>& U1, FieldPolynomial<T>& V1,
			const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q,
			const long d)
	throw(NotInSameFieldException);
/****************** Class FieldPolynomial ******************/
	/**
	 * \ingroup Fields
	 * \brief An polynomial with coefficients over a finite field.
	 *
	 * Objects of this class represent polynomials over a finite field, as represented by the class Field.
	 * With the exception of the zero polynomial created by the \link FieldPolynomial() default constructor\endlink,
	 * any element has an unique \e parent \e field and binary operations can combine two elements only in one of the
	 * following two cases:
	 *  - the two polynomials have the same parent field,
	 *  - the parent element of one polynomial is the \link Field::primeField() prime field\endlink of the other's.
	 *
	 * Polynomials created through the \link FieldPolynomial() default constructor\endlink, as for example
	 * \code
	 * FieldPolynomial<T> P;
	 * \endcode
	 * have a special status as they don't belong to any specific field: their default value is 0 and they
	 * can be combined with any other polynomial. The result of a binary operation involving such a special
	 * polynomial is one of the following:
	 *  - A DivisionByZeroException if the operation is division and the divisor is the special 0 polynomial.
	 *  - A polynomial over the field \b K, if \b K is the parent field of the other polynomial.
	 *  - The special 0 polynomial if the other polynomial is the special 0 polynomial.
	 *  - An UndefinedFieldException if the other polynomial is not the special 0 polynomial and yet does
	 *    not belong to any field, as in this example
	 *    \code
	 *    FieldPolynomial<T> P;
	 *    return P + 1;
	 *    \endcode
	 *    See UndefinedFieldException for more details.
	 *
	 * The way the arithmetics of the field are actually implemented is
	 * given by the template parameter \a T that must be one of the \ref Infrastructures.
	 * Note that changing the Infrastructure may sensibly change the speed of your code.
	 *
	 * \tparam T An \ref Infrastructures "Infrastructure". It specfies which \NTL types will carry out
	 * the arithmetic operations.
	 *
	 * \see Field, FieldElement, UndefinedFieldException
	 */
	template <class T> class FieldPolynomial {

	friend class Field<T>;
	friend class FieldElement<T>;

	/****************** GCD ******************/
	/**
	 * \brief GCD of \a P and \a Q.
	 *
	 * \throw NotInSameFieldException If \a P and \a Q do not have the same \parent.
	 * \relates FieldPolynomial
	 */
	friend FieldPolynomial<T> GCD<T>(const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q)
	throw(NotInSameFieldException);
	/**
	 * \brief XGCD of \a P and \a Q.
	 *
	 * \param [in] P A polynomial.
	 * \param [in] Q A polynomial having the same \parent as P.
	 * \param [out] U A polynomial to hold the result.
	 * \param [out] V A polynomial to hold the result.
	 * \param [out] G The GCD of this polynomial and \a Q.
	 * \throw NotInSameFieldException If \a P and \a Q do not have the same \parent.
	 * \invariant at the end of the method the following relation holds:
	 * \code
	 * U*P + V*Q == G;
	 * \endcode
	 * \relates FieldPolynomial
	 */
	friend void XGCD<T>(const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q,
	FieldPolynomial<T>& U, FieldPolynomial<T>& V, FieldPolynomial<T>& G)
	throw(NotInSameFieldException);

	/**
	 * \brief Half GCD of \a P and \a Q.
	 *
	 * \param [out] U0 A polynomial to hold the result.
	 * \param [out] V0 A polynomial to hold the result.
	 * \param [out] U1 A polynomial to hold the result.
	 * \param [out] V1 A polynomial to hold the result.
	 * \param [in] P A polynomial.
	 * \param [in] Q A polynomial having the same \parent as P.
	 * \param [in] d A bound on the degree of the result.
	 * \throw NotInSameFieldException If \a P and \a Q do not have the same \parent.
	 * \throw BadParametersException If \f$d>\max(\deg P,\deg Q)\f$ or \f$d<0\f$.
	 * \invariant at the end of the method the following relation holds:
	 * \f[\left(\begin{array}{cc}U_0&V_0\\U_1&V_1\end{array}\right)
	 * \left(\begin{array}{c}P\\Q\end{array}\right) =
	 * \left(\begin{array}{c}R_j\\R_{j+1}\end{array}\right)\f]
	 * where \f$R_j\f$ and \f$R_{j+1}\f$ are the reminders in the XGCD computation of
	 * \a P and \a Q such that \f$\deg R_{j+1}\le\max(\deg P,\deg Q)-d<\deg\R_j\f$.
	 * \relates FieldPolynomial
	 */
	friend void
	HalfGCD<T>(FieldPolynomial<T>& U0, FieldPolynomial<T>& V0,
				FieldPolynomial<T>& U1, FieldPolynomial<T>& V1,
				const FieldPolynomial<T>& P, const FieldPolynomial<T>& Q,
				const long d)
	throw(NotInSameFieldException, BadParametersException);

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
		/** \brief The \NTL representation of this polynomial if the \parent is a prime field. */
		GFpX repBase;
		/** \brief The \NTL representation of this polynomial if the \parent is an extension field. */
		GFpEX repExt;
		/** \brief Whether the \parent of this polynomial is a prime or an extension field. */
		bool base;
		/** \brief The \parent. NULL if no \parent. */
		const Field<T>* parent_field;
	/** \endcond */


	public:
	/****************//** \name Constructors ******************/
	/** @{ */
		/**
		 * \brief Construct the special 0 polynomial.
		 *
		 * The special 0 polynomial has no coefficient field, yet it can be added, multiplied, etc.
		 * to any other FieldPolynomial. See the \link FieldPolynomial introduction \endlink for more details.
		 * If you want to construct the 0 polynomial of a specific field, use Field::zero() in
		 * conjunction with FieldPolynomial(const FieldElement<T>&) instead.
		 *
		 * \see UndefinedFieldException, Field::zero(), FieldPolynomial(const FieldElement<T>&).
		 */
		FieldPolynomial() throw() : parent_field(NULL) {}
	/** @} */
	/****************//** \name Properties ******************/
	/** @{ */
		/**
		 * \brief The \e parent \e field.
		 *
		 * With the exception of the zero polynomial created by the \link FieldPolynomial() default constructor\endlink,
		 * any polynomial has an unique \e parent \e field and binary operations can combine two polynomials only in one of the
		 * following two cases:
		 *  - the two polynomials have the same parent field,
		 *  - the parent element of one polynomial is the \link Field::primeField() prime field\endlink of the other's.
		 *
		 * \throw UndefinedFieldException If this polynomial has no coefficient field.
		 * \see FieldElement(), UndefinedFieldException.
		 */
		const Field<T>& parent() const
		throw(UndefinedFieldException) {
			if (!parent_field)
				throw UndefinedFieldException();
			return *parent_field;
		}
		/** \brief Degree of the polynomial.
		 *
		 * \return The degree of the polynomial or -1 if the polynomial is 0.
		 */
		long degree() const throw();
	/** @} */

	/****************//** \name Copy
	 * You can copy polynomials using assignment.
	 * @{
	 */
		FieldPolynomial(const FieldPolynomial<T>&) throw();
		FieldPolynomial<T>& operator=(const FieldPolynomial<T>&) throw();
		/* \brief Degree 0 polynomial with constant coefficient \a e. */
		FieldPolynomial(const FieldElement<T>& e) throw();
		/* \brief Degree 0 polynomial with constant coefficient \a e. */
		FieldPolynomial<T>& operator=(const FieldElement<T>& e) throw();
		/**
		 * \brief The degree 0 polynomial with constant coefficient \a i.
		 *
		 * The \parent of the polynomial does not change through the assignment.
		 *
		 * \return A reference to the result.
		 * \throw UndefinedFieldException If this is the \link FieldPolynomial() special 0 element \endlink
		 * and \a i is different from 0.
		 */
		FieldPolynomial<T>& operator=(const BigInt& i)
		throw(UndefinedFieldException);
	/** @} */

	/****************//** \name Access to coefficients ******************/
	/** @{ */
		/** \brief Store in \a e the <i>i</i>-th coefficient.
		 *
		 * \param [in] i A positive integer.
		 * \param [out] e A FieldElement to hold the result.
		 * \throw BadParametersException If \a i is negative.
		 */
		void getCoeff(const long i, FieldElement<T>& e)
		const throw(BadParametersException);
		/**
		 * \brief Set the <i>i</i>-th coefficient to \a e.
		 *
		 * If this polynomial is the special 0 polynomial, its \parent becomes
		 * the \parent of \a e after this call. If moreover \a e is the
		 * \link FieldElement::FieldElement() special 0 element\endlink, this method does nothing.
		 *
		 * \param [in] i A positive integer.
		 * \param [out] e A FieldElement having the same \parent as this polynomial or the
		 * \link FieldElement::FieldElement() special 0 element\endlink.
		 * \throw BadParametersException If \a i is negative.
		 * \throw NotInSameFieldException If the \link parent() parent fields\endlink of \a e and this
		 * polynomial are both defined and differ.
		 */
		void setCoeff(const long i, const FieldElement<T>& e)
		throw(NotInSameFieldException, BadParametersException);
		/**
		 * \brief Set the <i>i</i>-th coefficient to \a c.
		 *
		 * \param [in] i A positive integer.
		 * \param [out] c An integer.
		 * \throw BadParametersException If \a i is negative.
		 * \throw UndefinedFieldException If this is the special 0 polynomial and \a c is different
		 * from 0.
		 */
		void setCoeff(const long i, const BigInt& c)
		throw(UndefinedFieldException, BadParametersException);
		/**
		 * \brief Set the <i>i</i>-th coefficient to 1.
		 *
		 * \param [in] i A positive integer.
		 * \throw BadParametersException If \a i is negative.
		 * \throw UndefinedFieldException If this is the special 0 polynomial.
		 */
		void setCoeff(const long i)
		throw(UndefinedFieldException, BadParametersException);
	/** @} */

	/****************** Arithmetics ******************/
		/** \name Binary operators
		 * All binary operators throw a NotInSameFieldException if neither of this two conditions
		 * is satisfied:
		 *  - the two operands have the same \parent,
		 *  - the \parent of one operand is the prime field of the other's.
		 * @{
		 */
		FieldPolynomial<T> operator+(const FieldPolynomial<T>& e)
		const throw(NotInSameFieldException) {
			FieldPolynomial<T> tmp = *this;
			tmp += e;
			return tmp;
		}
		FieldPolynomial<T> operator-(const FieldPolynomial<T>& e)
		const throw(NotInSameFieldException) {
			FieldPolynomial<T> tmp = *this;
			tmp -= e;
			return tmp;
		}
		FieldPolynomial<T> operator*(const FieldPolynomial<T>& e)
			const throw(NotInSameFieldException) {
			FieldPolynomial<T> tmp = *this;
			tmp *= e;
			return tmp;
		}
		FieldPolynomial<T> operator/(const FieldPolynomial<T>& e)
			const throw(NotInSameFieldException, DivisionByZeroException) {
			FieldPolynomial<T> tmp = *this;
			tmp /= e;
			return tmp;
		}
		FieldPolynomial<T> operator%(const FieldPolynomial<T>& e)
			const throw(NotInSameFieldException, DivisionByZeroException) {
			FieldPolynomial<T> tmp = *this;
			tmp %= e;
			return tmp;
		}

		/* Self-incrementing binary operations. */
		void operator+=(const FieldPolynomial<T>&)
			throw(NotInSameFieldException);
		void operator-=(const FieldPolynomial<T>&)
			throw(NotInSameFieldException);
		void operator*=(const FieldPolynomial<T>&)
			throw(NotInSameFieldException);
		void operator/=(const FieldPolynomial<T>&)
			throw(NotInSameFieldException, DivisionByZeroException);
		void operator%=(const FieldPolynomial<T>&)
			throw(NotInSameFieldException, DivisionByZeroException);

		/** \brief Stores \a a + \a b in this polynomial. */
		void sum(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator+=(b); }
		/** \brief Stores \a a - \a b in this polynomial. */
		void difference(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator-=(b); }
		/** \brief Stores \a a * \a b in this polynomial. */
		void product(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator*=(b); }
		/** \brief Stores \a a / \a b in this polynomial. */
		void division(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException, DivisionByZeroException)
		{ operator=(a); operator/=(b); }
		/** \brief Stores \a a % \a b in this polynomial. */
		void mod(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException, DivisionByZeroException)
		{ operator=(a); operator%=(b); }

		/** @} */


		/** \name Unary operators */
		/** @{ */
		/** \brief Test divisibility of this polynomial by \a P. */
		bool divides(const FieldPolynomial<T>& P) const throw();
		FieldPolynomial<T> operator-() const throw() {
			FieldPolynomial<T> tmp = *this;
			tmp.negate();
			return tmp;
		}
		/** \brief Power. */
		FieldPolynomial<T> operator^(const long i) const throw() {
			FieldPolynomial<T> tmp = *this;
			tmp ^= i;
			return tmp;
		}
		FieldPolynomial<T> derivative() const throw() {
			FieldPolynomial<T> tmp = *this;
			tmp.self_derivative();
			return tmp;
		}
		/** \brief The monic polynomial obtained by dividing this polynomial by
		 * its leading coeffcient. */
		FieldPolynomial<T> monic() const throw() {
			FieldPolynomial<T> tmp = *this;
			tmp.normalize();
			return tmp;
		}
		/** \brief <i>p</i>-th power (frobenius morphism) of the coefficients.
		 *
		 * Apply the frobenius morphism to the coefficients of this polynomial.
		 * \see FieldElement::frobenius().
		 */
		FieldPolynomial<T> frobenius() const throw() {
			FieldPolynomial<T> tmp = *this;
			tmp.self_frobenius();
			return tmp;
		}
		/** \brief <i>p<sup>n</sup></i>-th power (iterated frobenius morphism) of the coefficients.
		 *
		 * Apply the iterated frobenius morphism to the coefficients of this polynomial.
		 * \see FieldElement::frobenius(const long) const.
		 */
		FieldPolynomial<T> frobenius(const long n) const throw() {
			FieldPolynomial<T> tmp = *this;
			tmp.self_frobenius(n);
			return tmp;
		}
		/** @} */

		/** \name Self-incrementing Unary operations
		 * These operators store the result of the operation into the polynomial itself.
		 * @{
		 */
		/** \brief Flip the sign of this polynomial. */
		void negate() throw();
		/** \brief \copybrief operator^() */
		void operator^=(const long) throw();
		/** \brief Derivative. */
		void self_derivative() throw();
		/** \brief Make this polynomial monic.
		 * \see monic().
		 */
		void normalize() throw();
		/** \brief \copybrief frobenius() const */
		void self_frobenius() throw();
		/** \brief \copybrief frobenius(const long) const */
		void self_frobenius(long) throw();
		/** @} */

	/****************//** \name Coercion of polynomials
	 * Every polynomial has an unique \parent, but its coefficients may belong to some subfield of its \parent
	 * or you may want to change its \parent to an overfield of its actual one. Coercion does
	 * exactly this. When the coericion is impossible because either there is no known embedding
	 * between the two fields or because the coefficients do not belong to the new field, an
	 * IllegalCoercionException is thrown.
	 * @{
	 */
	 	/**
	 	 * \brief Coerce to a scalar polynomial.
	 	 *
	 	 * Coerce to a polynomial with coefficients in F<sub>p</sub>.
	 	 *
	 	 * \return The newly created polynomial.
	 	 * \throw IllegalCoercionException If the polynomial is not a scalar polynomial.
	 	 * \invariant This is the same as doing
	 	 * \code
	 	 * *this >> parent().primeField();
	 	 * \endcode
	 	 */
		FieldPolynomial<T> toScalarPolynomial() const throw(IllegalCoercionException);
		/**
		 * \brief Coerce to the field \a F.
		 *
		 * \param [in] F A finite subfield or overfield of the \parent, containing the coefficients of
		 * the polynomial.
		 * \return The newly created polynomial.
		 * \throw IllegalCoercionException If no embedding is known between the \parent and \a F or
		 * if the coefficients do not belong to \a F.
		 */
		FieldPolynomial<T> operator>>(const Field<T>&) const throw(IllegalCoercionException);
		/**
		 * \brief Coerce to the field \a F and store the result in this polynomial.
		 *
		 * \param [in] F A finite subfield or overfield of the \parent, containing the coefficients of
		 * the polynomial.
		 * \throw IllegalCoercionException If no embedding is known between the \parent and \a F or
		 * if the coefficients do not belong to \a F.
		 */
		void operator>>=(const Field<T>& F) throw(IllegalCoercionException) {
			FieldPolynomial<T> tmp = *this >> F;
			*this = tmp;
		}
		/**
		 * \brief Test if this polynomial is coercible to \a F.
		 */
		bool isCoercible(const Field<T>&) const throw();
	/** @} */

	/****************//** \name Evaluation ******************/
	/** @{ */
		/**
		 * \brief The evaulation of this polynomial at \a e.
		 * \see FieldElement::evaluate(const FieldPolynomial<T>&, vector<FieldPolynomial<T> >&) const.
		 */
		FieldElement<T> evaluate(const FieldElement<T>& e,
		vector<FieldPolynomial<T> >& minpols)
		const throw(IllegalCoercionException, BadParametersException)
		{ return e.evaluate(*this, minpols); }

		/**
		 * \brief \copybrief evaluate(constFieldElement<T>&, vector<FieldPolynomial<T> >&) const
		 * \see FieldElement::evaluate(const FieldPolynomial<T>&) const.
		 */
		FieldElement<T> evaluate(const FieldElement<T>& e)
		const throw(IllegalCoercionException, BadParametersException)
		{ return e.evaluate(*this); }
	/** @} */
	/****************//** \name Predicates ******************/
	/** @{ */
		/**
		 * \brief Equality.
		 *
		 * This method does not try to coerce the polynomials to the same field to test equality.
		 * \throw NotInSameFieldException If the two polynomials do not have de same \parent.
		 */
		bool operator==(const FieldPolynomial<T>&) const throw(NotInSameFieldException);
		/** \brief Equality. */
		bool operator==(const FieldElement<T>&) const throw(NotInSameFieldException);
		/** \brief Equality. */
		bool operator==(const BigInt&) const throw();

		/**
		 * \brief Inequality.
		 *
		 * This method does not try to coerce the polynomials to the same field to test equality.
		 * \throw NotInSameFieldException If the two polynomials do not have de same \parent.
		 */
		bool operator!=(const FieldPolynomial<T>& e)
		const throw(NotInSameFieldException)
		{ return !(*this==e); }
		/** \brief Inequality. */
		bool operator!=(const FieldElement<T>& e)
		const throw(NotInSameFieldException)
		{ return !(*this==e); }
		/** \brief Inequality. */
		bool operator!=(const BigInt& i) const throw()
		{ return !(*this==i); }
		/** \brief Test to zero. */
		bool isZero() const throw() {
			return !parent_field ||
				(base ? IsZero(repBase) : IsZero(repExt));
		}
		/** \brief Test to one. */
		bool isOne() const throw() {
			return parent_field &&
				(base ? IsOne(repBase) : IsOne(repExt));
		}
		/** \brief Test if this polynomial has coefficients in F<sub>p</sub>.
		 * \see toScalarPolynomial().
		 */
		bool isScalarPolynomial() const throw();
	/** @} */

	/****************//** \name Access to the Infrastructure
	 * These methods let you access the internal \NTL representation of polynomials.
	 *
	 * \warning These methods are for advanced use only. Use them
	 * if you want to use an algorithm by you or in NTL that is not available for
	 * FieldPolynomial.
	 * @{
	 */
		/**
		 * \brief Get the representation of polynomials whose \parent is F<sub>p</sub>.
		 *
		 * \param [out] P An \NTL scalar polynomial to hold the result.
		 * \throw IllegalCoercionException If the \parent is not a prime field
		 * \note This method automatically switches the context to the \parent context.
		 * See Field::switchContext().
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink,
		 * Field::fromInfrastructure(), Field::switchContext().
		 */
		void toInfrastructure(GFpX& P) const throw(IllegalCoercionException);
		/**
		 * \brief Get the representation of polynomials whose \parent is an extension field.
		 *
		 * \param [out] P An \NTL polynomial to hold the result.
		 * \throw IllegalCoercionException If the \parent is a prime field
		 * \note This method automatically switches the context to the \parent context.
		 * See Field::switchContext().
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink,
		 * Field::fromInfrastructure(), Field::switchContext().
		 */
		void toInfrastructure(GFpEX& P) const throw(IllegalCoercionException);
	/** @} */

	/****************//** \name Printing ******************/
	/** @{ */
		/** \brief Print this polynomial to \a o */
		ostream& print(ostream&) const;
		/**
		 * \brief Print this element to \a o as a polynomial over its \parent.
		 *
		 * Print as a polynomial in the variable \a varPoly
		 * and print elements of the \parent as polynomials over F<sub>p</sub>
		 * in the variable \a varField. The coefficients of the polynomial are
		 * printed as through \link FieldElement::print(ostream&, const string&) const \c FieldElement::print(o,varField) \endlink.
		 *
		 * \param [in,out] o An output stream.
		 * \param [in] varPoly A variable name.
		 * \param [in] varField A variable name.
		 * \return A pointer to the modified stream.
		 * \see FieldElement::print(ostream&, const string&) const.
		 */
		ostream& print(ostream&, const string& varPoly, const string& varField) const;
		/**
		 * \brief Print this element to \a o as a polynomial over its \parent.
		 *
		 * Print as a polynomial in the variable \a varPoly
		 * and print elements of the \parent as multivariate polynomials over F<sub>p</sub>
		 * in the variables \a varsField. The coefficients of the polynomial are
		 * printed as through
		 * \link FieldElement::print(ostream&, const vector<string>&) const \c FieldElement::print(o,varsField) \endlink.
		 *
		 * \param [in,out] o An output stream.
		 * \param [in] varPoly A variable name.
		 * \param [in] varsField A vector of variable names.
		 * \return A pointer to the modified stream.
		 * \see FieldElement::print(ostream&, const vector<string>&) const.
		 */
		ostream& print(ostream&, const string& varPoly, const vector<string>& varsField) const;
	/** @} */
	/****************** Destructor ******************/
		~FieldPolynomial() throw() {}


	/*****************************************************/
	/****************** Private section ******************/
	/*****************************************************/

	private:
	/** \cond DEV */
	/****************//** \name Internal Constructors
	 * Construct a polynomial with given representation and \parent.
	 * Reserved for used by Field.
	 * @{
	 */
		FieldPolynomial(const Field<T>* p, const GFpX& PBase, const GFpEX& PExt, const bool b) throw() :
			repBase(PBase), repExt(PExt), base(b), parent_field(p) {}
		FieldPolynomial(const Field<T>* p, const GFpEX& P) throw() :
			repExt(P), base(false), parent_field(p) {}
		FieldPolynomial(const Field<T>* p, const GFpX& P) throw() :
			repBase(P), base(true), parent_field(p) {}
	/** @} */
	/****************** Utility Routines ******************/
		/**
		 * \brief Check if \a e has the same \parent as this polynomial.
		 *
		 * If the two have the same \parent, this method does nothing,
		 * otherwise it throws a NotInSameFieldException.
		 */
		void sameLevel(const FieldElement<T>& e) const
		throw(NotInSameFieldException) {
			if (parent_field != e.parent_field)
			throw NotInSameFieldException();
		}
		/**
		 * \brief Check if \a e has the same \parent as this polynomial.
		 *
		 * If the two have the same \parent, this method does nothing,
		 * otherwise it throws a NotInSameFieldException.
		 */
		void sameLevel(const FieldPolynomial<T>& e) const
		throw(NotInSameFieldException) {
			if (parent_field != e.parent_field)
			throw NotInSameFieldException();
		}
	/** \endcond */
	};
/****************** Printing ******************/
	/** \brief Print \a P to \a o.
	 * \relates FieldPolynomial
	 */
	template <class T> ostream&
	operator<<(ostream& o, const FieldPolynomial<T>& P) {
		return P.print(o);
	}
}

#endif /*FIELDPOLYNOMIAL_H_*/
