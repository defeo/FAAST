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
#ifndef FIELDELEMENT_H_
#define FIELDELEMENT_H_

#include "Exceptions.hpp"
#include <string>
#include <vector>

namespace FAAST {
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
	 * any element has an unique \e parent \e field and binary operations can combine two elements only in one of the
	 * following two cases:
	 *  - the two elements have the same parent field,
	 *  - the parent element of one element is the \link Field::primeField() prime field\endlink of the other's.
	 *
	 * Elements created through the \link FieldElement() default constructor\endlink, as for example
	 * \code
	 * FieldElement<T> elt;
	 * \endcode
	 * have a special status as they don't belong to any specific field: their default value is 0 and they
	 * can be combined with any other element. The result of a binary operation involving such a special
	 * element is one of the following:
	 *  - A DivisionByZeroException if the operation is division and the divisor is the special 0 element.
	 *  - An element of field \b K, if \b K is the parent field of the other element.
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
	 * If the \parent of \a e is the base field of an Artin-Schreier tower, then \a v is filled with the
	 * coefficients in F<sub>p</sub> of its univariate representation. Otherwise let
	 * \code
	 * x = e.parent().primitiveElement();
	 * \endcode
	 * This method fills the vector \a v with \ref Field::p "p" elements of
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
	 * \throw NoSubFieldException If F<sub>p</sub> is the \parent of \a e.
	 * \note The result always lies in the primitive tower (the stem), even if \a e does not.
	 * \invariant If the \parent of \a e is in the primitive tower (the stem), this is equivalent to
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
	 * If the \parent of the elements of \a v is the prime field, then \a e is the element whose univariate
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
	 * the \parent of the elements of \a v, this corresponds to convert \a v from the
	 * multivariate
	 * representation as an element of \b K[\a x] to the internal (univariate) representation.
	 * This routine implements the algorithm \c LiftUp of [\ref ISSAC "DFS '09", Section 4.4].
	 *
	 * \param [in] v A vector of elements all in the same field.
	 * \param [out] e An element satisfying condition (2).
	 * \throw NotInSameFieldException If the elements of \a v do not all
	 * have the same \parent.
	 * \throw NoOverFieldException If the \parent of the elements of \a v has no
	 * \link Field::overField() overfield\endlink.
	 * \note The result always lies in the primitive tower (the stem), even if the elements of \a v do not.
	 * \invariant If the \parent of the elements of \a v is in the primitive tower (the stem),
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
		/** \brief The \NTL representation of this element if the \parent is a prime field. */
		GFp repBase;
		/** \brief The \NTL representation of this element if the \parent is an extension field. */
		GFpE repExt;
		/** \brief Whether this element belongs to a prime or an extension field. */
		bool base;
		/** \brief The \parent. NULL if no \parent. */
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
		 * \brief The \e parent \e field.
		 *
		 * With the exception of the zero element created by the \link FieldElement() default constructor\endlink,
		 * any element has an unique \e parent \e field and binary operations can combine two elements only in one of the
		 * following two cases:
		 *  - the two elements have the same parent field,
		 *  - the parent element of one element is the \link Field::primeField() prime field\endlink of the other's.
		 *
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
		 * The \parent of the element does not change through the assignment.
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
		 * All binary operators throw a NotInSameFieldException if neither of this two conditions
		 * is satisfied:
		 *  - the two operands have the same \parent,
		 *  - the \parent of one operand is the prime field of the other's.
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
		 * This method is a generalization of the algorithm \c IterFrobenius of [\ref ISSAC "DFS '09", Section 5].
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
		 * This method is a generalization of the algorithm \c Pseudotrace of [\ref ISSAC "DFS '09", Section 5].
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
		/** \brief The minimal polynomial over F<sub>p</sub>. */
		FieldPolynomial<T> minimalPolynomial() const throw() {
			parent_field->switchContext();
			GFpX minpol;
			MinPolyMod(minpol, rep(repExt), GFpE::modulus());
			return parent_field->primeField().fromInfrastructure(minpol);
		}
		/**
		 * \brief The minimal polynomial over the field \a F.
		 *
		 * This method implements a yet unpublished algorithm to compute minimal polynomials
		 * in Artin-Schreier towers. It only works when \a F is a field of an Artin-Schreier
		 * tower as constructed in [\ref ISSAC "DFS '09", Section 3].
		 *
		 * \param [in] F A subfield of the \parent of this element.
		 * \return A polynomial over \a F being the minimal polynomial of this element.
		 * \throws NotASubFieldException If the \parent is not an extension field of \a F.
		 * \throws NotSupportedException If \a F is a prime field not being the base field of an Artin-Schreier
		 * tower. Use minimalPolynomial() instead.
		 */
		FieldPolynomial<T> minimalPolynomial(const Field<T>& F)
		const throw(NotASubFieldException, NotSupportedException) {
			vector<FieldPolynomial<T> > minpols;
			minimalPolynomials(F, minpols);
			return minpols[0];
		}
		/**
		 * \brief All the minimal polynomials up to the field \a F.
		 *
		 * The result is the same as doing
		 * \code
		 * res[0] = minimalPolynomial(F);
		 * res[1] = minimalPolynomial(F.overField());
		 * res[2] = minimalPolynomial(F.overField().overField());
		 * ...
		 * \endcode
		 * up to \c minimalPolynomial(parent()). The computation is more efficient, though.
		 *
		 * \param [in] F A subfield of the \parent field.
		 * \param [out] res The vector is filled with the minimal polynomials of this element over
		 * the intermediate extension fields of F. All previous data are discarded.
		 * \throws NotASubFieldException If the \parent is not an extension field of \a F.
		 * \throws NotSupportedException If \a F is a prime field not being the base field of an Artin-Schreier
		 * tower. Use minimalPolynomial() instead.
		 * \see minimalPolynomial(const Field<T>&) const.
		 */
		void minimalPolynomials(const Field<T>& F, vector<FieldPolynomial<T> >& res)
		const throw(NotASubFieldException, NotSupportedException);

		/**
		 * \brief \copybrief affineMinimalPolynomial(const Field<T>&, const FieldElement<T>&) const
		 *
		 * This is an optimized version of affineMinimalPolynomial(const Field<T>&, const FieldElement<T>&) const
		 * that lets you pass an additional parameter \a minpols containing some precomputed quantities.
		 *
		 * The optional parameter \a minpols must contain either the result of
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink,
		 * or must be an empty vector, in which
		 * case it is filled with the result of
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink.
		 * In any other case either a
		 * BadParametersException is thrown or the behaviour is undefined.
		 *
		 * This method should be preferred when you care about efficiency and you do many calls to
		 * affineMinimalPolynomial() and evaluate().
		 *
		 * \param [in] F A subfield of the \parent.
		 * \param [in] a An element contained in the field F[\a x], where \a x is this element.
		 * \param [in,out] minpols If this vector is empty, then
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink
		 * is called, otherwise it is assumed to contain the data computed by the former method call.
		 * \return The <i>a</i>-affine minimal polynomial over the field \a F.
		 * \throws NotASubFieldException If the \parent is not an extension field of \a F.
		 * \throws NotSupportedException If \a F is a prime field not being the base field of an Artin-Schreier
		 * tower.
		 * \throws NoSuchPolynomialException If no such polynomial exists. Equivalently, if
		 * \a a is not an element of F[\a x] where \a x is this element.
		 * \throws BadParametersException If \a minpols is not empty or does not contain the result of
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink.
		 * \see affineMinimalPolynomial(const Field<T>&, const FieldElement<T>&) const.
		 */
		FieldPolynomial<T> affineMinimalPolynomial(
		const Field<T>& F, const FieldElement<T>& a,
		vector<FieldPolynomial<T> >& minpols)
		const throw(NotASubFieldException, NoSuchPolynomialException,
		NotSupportedException, BadParametersException);

		/**
		 * \brief The <i>a</i>-affine minimal polynomial over the field \a F.
		 *
		 * That is the minimum degree polynomial \a P of \a F[X] such that
		 * \f[ P(x) = \mathtt{a}\mathrm{,} \f] where \a x is this
		 * element. Observe that this is the same as computing the interpolation polynomial such that
		 * \f[ P\bigl(\phi_F^n(x)\bigr) = \phi_F^n(\mathtt{a}) \quad \forall n\mathrm{,} \f]
		 * where \f$ \phi_F \f$ is the frobenius morphism fixing \a F.
		 * This implements a yet unpublished algorithm, it only works when \a F is a field of an Artin-Schreier
		 * tower as constructed in [\ref ISSAC "DFS '09", Section 3].
		 *
		 * \param [in] F A subfield of the \parent.
		 * \param [in] a An element contained in the field F[\a x], where \a x is this element.
		 * \return The <i>a</i>-affine minimal polynomial over the field \a F.
		 * \throws NotASubFieldException If the \parent is not an extension field of \a F.
		 * \throws NotSupportedException If \a F is a prime field not being the base field of an Artin-Schreier
		 * tower.
		 * \throws NoSuchPolynomialException If no such polynomial exists. Equivalently, if \a a is not
		 * an element of the field F[\a x], where \a x is this element.
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink.
		 */
		FieldPolynomial<T> affineMinimalPolynomial(
		const Field<T>& F, const FieldElement<T>& a)
		const throw(NotASubFieldException, NoSuchPolynomialException,
		NotSupportedException) {
			vector<FieldPolynomial<T> > minpols;
			return affineMinimalPolynomial(F,a,minpols);
		}

		/**
		 * \brief \copybrief evaluate(const FieldPolynomial<T>&) const
		 *
		 * This is an optimized version of evaluate(const FieldPolynomial<T>&) const
		 * that lets you pass an additional parameter \a minpols containing some precomputed quantities.
		 *
		 * The optional parameter \a minpols must contain either the result of
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink,
		 * or must be an empty vector, in which
		 * case it is filled with the result of
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink.
		 * In any other case either a
		 * BadParametersException is thrown or the behaviour is undefined.
		 *
		 * This method should be preferred when you care about efficiency and you do many calls to
		 * affineMinimalPolynomial() and evaluate().
		 *
		 * \param [in] P A polynomial with coefficients in a subfield or overfield of the \parent.
		 * \param [in,out] minpols If this vector is empty, then
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink
		 * is called, otherwise it is assumed to contain the data computed by the former method call.
		 * \return The evaluation of \a P at this element.
		 * \throws IllegalCoercionException If this element cannot be coerced to the coefficient field
		 * of \a P nor can \a P be coerced to the \parent.
		 * \throws BadParametersException If \a minpols is not empty or does not contain the result of
		 * \link minimalPolynomials(const Field<T>&, vector<FieldPolynomial<T> >&) const \c minimalPolynomials(F,minpols) \endlink.
		 * \invariant This is the same as \link FieldPolynomial::evaluate() \c P.evaluate(*this, minpols) \endlink.
		 * \see evaluate(const FieldPolynomial<T>&) const, FieldPolynomial::evaluate().
		 */
		FieldElement<T> evaluate(const FieldPolynomial<T>& P,
		vector<FieldPolynomial<T> >& minpols)
		const throw(IllegalCoercionException, BadParametersException);

		/**
		 * \brief The evaulation of \a P at this element.
		 *
		 * This implements a yet unpublished algorithm for Artin-Schreier
		 * towers constructed as in [\ref ISSAC "DFS '09", Section 3]. It uses Horner
		 * evaluation scheme for all the other cases.
		 *
		 * The optional parameter minpols must contain either the result of
		 * this.minimalPolynomials(P.parent,v), or must be an empty vector,
		 * in which case it is filled with the result of
		 * minimalPolynomials(P.parent,v).
		 *
		 * \param [in] P A polynomial with coefficients in a subfield or overfield of the \parent.
		 * \return The evaluation of \a P at this element.
		 * \throws IllegalCoercionException If this element cannot be coerced to the coefficient field
		 * of \a P nor can \a P be coerced to the \parent.
		 * \invariant This is the same as \link FieldPolynomial::evaluate() \c P.evaluate(*this) \endlink.
		 * \see FieldPolynomial::evaluate().
		 */
		FieldElement<T> evaluate(const FieldPolynomial<T>& P)
		const throw(IllegalCoercionException) {
			vector<FieldPolynomial<T> > minpols;
			return evaluate(P, minpols);
		}
	/** @} */

	/****************//** \name Coercion of elements
	 * Every element has an unique \parent, but it may belong to some subfield of its \parent
	 * or you may want to change its \parent to an overfield of its actual one. Coercion does
	 * exactly this. When the coericion is impossible because either there is no known embedding
	 * between the two fields or because the element does not belong to the new field, an
	 * IllegalCoercionException is thrown.
	 * @{
	 */
	 	/**
	 	 * \brief Coerce to a scalar.
	 	 *
	 	 * Coerce to an element of F<sub>p</sub>.
	 	 *
	 	 * \return The newly created element.
	 	 * \throw IllegalCoercionException If the element is not a scalar.
	 	 * \invariant This is the same as doing
	 	 * \code
	 	 * *this >> parent().primeField();
	 	 * \endcode
	 	 */
		FieldElement<T> toScalar() const throw(IllegalCoercionException);
		/**
		 * \brief Coerce to the field \a F.
		 *
		 * \param [in] F A finite subfield or overfield of the \parent, containing this element.
		 * \return The newly created element.
		 * \throw IllegalCoercionException If no embedding is known between the \parent and \a F or
		 * if this element does not belong to \a F.
		 */
		FieldElement<T> operator>>(const Field<T>& F) const
		throw(IllegalCoercionException) {
			FieldElement<T> tmp = *this;
			tmp >>= F;
			return tmp;
		}
		/**
		 * \brief Coerce to the field \a F and store the result in this element.
		 *
		 * \param [in] F A finite subfield or overfield of the \parent, containing this element.
		 * \throw IllegalCoercionException If no embedding is known between the \parent and \a F or
		 * if this element does not belong to \a F.
		 */
		void operator>>=(const Field<T>&) throw(IllegalCoercionException);
		/**
		 * \brief Test if this element is coercible to \a F.
		 */
		bool isCoercible(const Field<T>&) const throw();
	/** @} */

	/****************//** \name Predicates ******************/
	/** @{ */
		/**
		 * \brief Equality.
		 *
		 * This method does not try to coerce the elements to the same field to test equality.
		 * \throw NotInSameFieldException If the two elements do not have de same \parent.
		 */
		bool operator==(const FieldElement<T>&) const throw(NotInSameFieldException);
		/** \brief Equality. */
		bool operator==(const BigInt&) const throw();
		/**
		 * \brief Inequality.
		 *
		 * This method does not try to coerce the elements to the same field to test equality.
		 * \throw NotInSameFieldException If the two elements do not have de same \parent.
		 */
		bool operator!=(const FieldElement<T>& e) const throw(NotInSameFieldException)
		{ return !(*this==e); }
		/** \brief Inequality. */
		bool operator!=(const BigInt& i) const throw() { return !(*this==i); }
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
		/** \brief Test if this element belongs to F<sub>p</sub>.
		 * \see toScalar().
		 */
		bool isScalar() const throw();
	/** @} */
	/****************//** \name Access to the Infrastructure
	 * These methods let you access the internal \NTL representation of elements.
	 *
	 * \warning These methods are for advanced use only. Use them
	 * if you want to use an algorithm by you or in NTL that is not available for
	 * FieldElement.
	 * @{
	 */
		/**
		 * \brief Get the representation of elements whose \parent is F<sub>p</sub>.
		 *
		 * \param [out] e An \NTL scalar element to hold the result.
		 * \throw IllegalCoercionException If the \parent is not a prime field
		 * \note This method automatically switches the context to the \parent context.
		 * See Field::switchContext().
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink,
		 * Field::fromInfrastructure(), Field::switchContext().
		 */
		void toInfrastructure(GFp& e) const throw(IllegalCoercionException);
		/**
		 * \brief Get the representation of elements whose \parent is an extension field.
		 *
		 * \param [out] e An \NTL element to hold the result.
		 * \throw IllegalCoercionException If the \parent is a prime field
		 * \note This method automatically switches the context to the \parent context.
		 * See Field::switchContext().
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink,
		 * Field::fromInfrastructure(), Field::switchContext().
		 */
		void toInfrastructure(GFpE& e) const throw(IllegalCoercionException);
	/** @} */

	/****************//** \name Printing ******************/
	/** @{ */
		/** \brief Print this element to \a o */
		ostream& print(ostream& o) const;
		/**
		 * \brief Print this element to \a o as a polynomial over F<sub>p</sub> in the
		 * variable \a var.
		 */
		ostream& print(ostream& o, const string& var) const;
		/**
		 * \brief Print this element to \a o as a multivariate polynomial over
		 * F<sub>p</sub>.
		 *
		 * The number of variables in \a vars must be at least
		 * one plus the \link Field::ArtinSchreierHeight Artin-Schreier height\endlink
		 * of the \parent. The recursive descent along the tower is done via
		 * Field::toBivariate(), the list of the fields involved in the descent
		 * is internally stored and reproduces backwards the list of calls to
		 * Field::ArtinSchreierExtension that have created the \parent.
		 *
		 * \param [in,out] o An output stream.
		 * \param [in] vars A vector of variable names.
		 * \return A pointer to the modified stream.
		 * \throws BadParametersException If there is not enough variables
		 *          in \a vars.
		 * \todo This method sucks!
		 */
		ostream& print(ostream& o, const vector<string>& vars) const;
	/** @} */
	/****************** Destructor ******************/
		~FieldElement() throw() {}


	/*****************************************************/
	/****************** Private section ******************/
	/*****************************************************/

	private:
	/** \cond DEV */
	/****************//** \name Helpers for frobenius and trace ******************/
	/** @{ */
		/** \brief <i>p<sup>j</sup>d</i>-th iterated frobenius.
		 *
		 * This is the algorithm \c IterFrobenius of [\ref ISSAC "DFS '09"].
		 */
		void BigFrob(const long j);
		/** \brief <i>n</i>-th iterated frobenius, for \a n < \a d.
		 *
		 * This implements a naive algorithm for iterated frobenius.
		 * \todo Implement a faster algorithm using modular composition.
		 */
		void SmallFrob(const long n);
		/** \brief <i>p<sup>j</sup>d</i>-th pseudotrace.
		 *
		 * This is the algorithm \c Pseudotrace of [\ref ISSAC "DFS '09"].
		 */
		void BigPTrace(const long j);
		/**
		 * \brief Put in \a v all the <i>p<sup>i</sup>d</i> pseudotraces for 0 <= i <= j
		 *
		 * Useful for precomputing pseudotrace as in [\ref ISSAC "DFS '09", Section 5].
		 */
		void BigPTraceVector(vector<FieldElement<T> >& v, const long j) const;
		/** \brief <i>n</i>-th pseudotrace, for \a n < \a d.
		 *
		 * This implements a naive algorithm for pseudotrace.
		 * \todo Implement a faster algorithm using modular composition.
		 */
		void SmallPTrace(const long n);
	/** @} */

	/****************//** \name Internal Constructors
	 * Construct an element with given representation and \parent.
	 * Reserved for used by Field.
	 * @{
	 */
		FieldElement(const Field<T>* p, const GFp& PBase, const GFpE& PExt, const bool b) throw() :
			repBase(PBase), repExt(PExt), base(b), parent_field(p) {}
		FieldElement(const Field<T>* p, const GFpE& P) throw() :
			repExt(P), base(false), parent_field(p) {}
		FieldElement(const Field<T>* p, const GFp& P) throw() :
			repBase(P), base(true), parent_field(p) {}
	/** @} */
	/****************** Utility Routines ******************/
		/**
		 * \brief Check if \a e has the same \parent as this element.
		 *
		 * If the two elements have the same \parent, this method does nothing,
		 * otherwise it throws a NotInSameFieldException.
		 */
		void sameLevel(const FieldElement<T>& e) const
		throw(NotInSameFieldException) {
			if (parent_field != e.parent_field)
			throw NotInSameFieldException();
		}
	/** \endcond */
	};

/****************** Printing ******************/
	/** \brief Print \a e to \a o.
	 * \relates FieldElement
	 */
	template <class T> ostream&
	operator<<(ostream& o, const FieldElement<T>& e) {
		return e.print(o);
	}
}

#endif /*FIELDELEMENT_H_*/
