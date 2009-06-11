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
#ifndef FIELD_H_
#define FIELD_H_

#include "Exceptions.hpp"
#include "FieldElement.hpp"
#include "FieldPolynomial.hpp"
#include <memory>

namespace FAAST {

#ifdef FAAST_TIMINGS
	/**
	 * \brief This struct stores various timings related to precomputations.
	 *
	 * \see [\ref ISSAC "DFS '09"] for details on the computations.
	 *
	 * \note The preprocessor flag \c FAAST_TIMINGS has to be passed to the compiler
	 * in order to activate this feature.
	 */
	typedef struct TIMINGS {
		/** \brief The time spent precomputing the (2p-1)-th cyclotomic polynomial */
		double CYCLOTOMIC;
		/** \brief The time spent precomputing the pseudotraces. See [\ref ISSAC "DFS '09", Section 5]. */
		double PSEUDOTRACES;
		/** \brief The time spent precomputing the inverse of the derivative of the defining polynomial.
		 * See [\ref ISSAC "DFS '09", Section 4]. */
		double LIFTUP;
		/** \brief The time spent precomputing the coefficients of the trace form over the bivariate
		 * basis. See [\ref ISSAC "DFS '09", Section 4]. */
		double TRACEVEC;
		/** \brief The time spent precomputing irreducible polynomials for base fields. */
		double BUILDIRRED;
		/** \brief The time spent testing the irreducibility of polynomials defining base fields. */
		double IRREDTEST;
		/** \brief The time spent testing primality of cardinalities */
		double PRIMETEST;
		/** \brief The time spent precomputing and inverting the matrix of the application
		 * X<sup>p</sup> - X in base fields. See [\ref ISSAC "DFS '09", Section 6]. */
		double ARTINMATRIX;
		/** \brief The time spent computing the primitive tower. See [\ref ISSAC "DFS '09", Section 3]. */
		double BUILDSTEM;

		TIMINGS() : CYCLOTOMIC(0),
		PSEUDOTRACES(0),
		LIFTUP(0),
		TRACEVEC(0),
		BUILDIRRED(0),
		IRREDTEST(0),
		PRIMETEST(0),
		ARTINMATRIX(0),
		BUILDSTEM(0)
		{}
	} TIMINGS;
#endif


	/**
	 * \defgroup Fields Finite Field Arithmetics
	 * This module contains three types, namely FAAST::Field, FAAST::FieldElement and
	 * FAAST::FieldPolynomial, which represent, respectively, finite fields, elements
	 * of finite fields and polynomials with coefficient over finite fields.
	 *
	 * With some notable exception (see
	 * \link FAAST::FieldElement::FieldElement() FieldElement()\endlink and
	 * \link FAAST::FieldPolynomial::FieldPolynomial() FieldPolynomial()\endlink),
	 * each FAAST::FieldElement and FAAST::FieldPolynomial
	 * belongs to a FAAST::Field object called its \e parent \e field.
	 * Elements (and polynomials) having the same parent field can be
	 * combined (added, multiplied, etc.) freely, while limitations
	 * apply for combining elements belonging to two different fields; see
	 * FAAST::FieldElement and FAAST::FieldPolynomial.
	 *
	 * Elements
	 * (and polynomials) can be moved from one field to another when a morphism
	 * is known; see \ref Field_lattices.
	 *
	 */

	/**
	 * \ingroup Fields
	 * \brief A finite field.
	 *
	 * Objects of this class can only be built through the static instantiators createField() and
	 * can never be destroyed.
	 *
	 * The way the arithmetics of the field are actually implemented is
	 * given by the template parameter \a T that must be one of the \ref Infrastructures.
	 * Note that changing the Infrastructure may sensibly change the speed of your code.
	 *
	 * \tparam T An \ref Infrastructures "Infrastructure". It specfies which \NTL types will carry out
	 * the arithmetic operations.
	 *
	 * \todo Field objects are immortal. A better memory management involving garbage collection may be
	 * implemented one day.
	 */
	template <class T> class Field {

	friend class FieldElement<T>;
	friend void pushDown<T>(const FieldElement<T>& e, vector<FieldElement<T> >& v) throw(NoSubFieldException);
	friend void liftUp<T>(const vector<FieldElement<T> >& v, FieldElement<T>& e) throw(NotInSameFieldException, NoOverFieldException);

#ifdef FAAST_TIMINGS
	public:
		/**
		 * \brief \copybrief TIMINGS See TIMINGS.
		 */
		static TIMINGS TIME;
#endif

	/** \name Local types
	 * Local types defined in this class. They are aliases to simplify the access
	 * to the \ref Infrastructures "Infrastructure" \a T and its subtypes.
	 *
	 * \see \ref Infrastructures
	 * @{
	 */
	public:
		/** \brief A link to the \ref Infrastructures Infrastructure. */
		typedef T Infrastructure;

	private:
		typedef typename T::GFp         GFp;
		typedef typename T::MatGFp      MatGFp;
		typedef typename T::VecGFp      VecGFp;
		typedef typename T::GFpX        GFpX;
		typedef typename T::GFpE        GFpE;
		typedef typename T::GFpEX       GFpEX;
		typedef typename T::BigInt      BigInt;
		typedef typename T::Context     Context;
		typedef typename T::GFpXModulus GFpXModulus;
	/** @} */

	/** \cond DEV */
	/****************//** \name Data members for stem fields ******************/
	/** @{ */
		/** \brief The immediate subfield, if it is defined */
		mutable const Field<T>* subfield;
		/** \brief The immediate overfield, if it is defined */
		mutable const Field<T>* overfield;
		/** \brief  Infrastructure-dependent data to perform computations in
		 * the field. (e.g. GF2EContext in NTL)
		 */
		Context context;
		/** \brief  The generator over GF(p) */
		const auto_ptr<const FieldElement<T> > primitive;
		/** \brief  Precomputed pseudotraces */
		mutable vector<FieldElement<T> > pseudotraces;
		/** \brief  Lift-up precomputation */
		mutable auto_ptr<const FieldElement<T> > liftuphelper;
		/** \brief  The inverse matrix of the d-1 minor of the
		 * linear application X<sup>p</sup>-X
		 */
		mutable MatGFp artin;
		/** \brief  The line we took away from X<sup>p</sup>-X to make
		 * it invertible
		 */
		mutable long artinLine;
		/** \brief  Flag related to the construction of the extension */
		const bool plusone;
		/** \brief \copybrief plusone */
		const bool twopminusone;
		/** \brief  The (2 \a p - 1)th cyclotomic polynomial */
		mutable auto_ptr<const Context> Phi;
	/** @} */

	/****************//** \name Data members for non-stem fields ******************/
	/** @{ */
		/** \brief  The stem-field to which this one is isomorphic */
		const Field<T>* stem;
	/** @} */

	/****************//** \name Common data members ******************/
	/** @{ */
		/** \brief  The virtual subfield, if it is different from the real.
		 * It is given mainly for printing purposes
		 */
		const Field<T>* vsubfield;
		/** \brief  The generator over the subfield */
		const auto_ptr<const FieldElement<T> > gen;
		/** \brief  The element of the subfield such that this field is
		 * the Artin-Schreier extension defined by the polynomial
		 * \f$ X^p - X - \mathtt{alpha} \f$
		 */
		const auto_ptr<const FieldElement<T> > alpha;
	/** @} */
	/** \endcond */

	public:
	/****************//** \name Public data members ******************/
	/** @{ */
		/** \brief The characteristic of the field. */
		const BigInt p;
		/** \brief The degree over the prime field F<sub>p</sub>. */
		const long d;
		/** \brief The constructed Artin-Schreier height.
		 *
		 * This is the number of intermediate Artin-Schreier extensions over
		 * baseField() that have been constructed using the the techniques of
		 * [\ref ISSAC "DFS '09", Sections 3 and 6].
		 *
		 * \invariant The following formula is always true:
		 * \f$ [\mathtt{K} : \mathtt{K.baseField()}] = p^\mathtt{height} \f$
		 */
		const long height;
	/** @} */


	public:
	/****************//** \name Instantiators
	 * All instantiators are static. There's no constructor for
	 * Field objects.
	 * @{
	 */

		/**
		 * \brief Default instantiator, builds a field from \NTL's context.
		 *
		 * \param [in] test If false, do not perform primality and irreducibility
		 * tests.
		 * \return A reference to the newly created Field object.
		 *
		 * \throws NotPrimeException If the current \NTL modulus, as obtained by
		 * \code
		 * 	T::GFp::modulus();
		 * \endcode
		 * is not a prime (unless \a T is GF2_Algebra).
		 * \throws NotIrreducibleException If the current \NTL modulus, as obtained by
		 * \code
		 * 	T::GFpE::modulus();
		 * \endcode
		 * is not an irreducible polynomial.
		 *
		 * \warning If the current modulus has degree \a d divisible by the characteristic and if its
		 * (\a d - 1)-th coefficient is 0, then you won't be able to construct Artin-Schreier
		 * extensions through a call to ArtinSchreierExtension(). In fact in this case the generator
		 * of the field has trace 0 and the costruction of [\ref ISSAC "DFS '09"] doesn't work.
		 */
		static const Field<T>& createField(const bool test = true)
		throw (NotPrimeException, NotIrreducibleException);
		/**
		 * \brief Build a field from an irreducible polynomial \a P.
		 *
		 * Build the field F[X]/\a P(X). \verbatim T::GFp::modulus \endverbatim must
		 * be set accordingly, unless \a T is GF2_Algebra.
		 *
		 * \param [in] P An irreducible polynomial.
		 * \param [in] test If false, do not perform primality and irreducibility
		 * tests.
		 * \return A reference to the newly created Field object.
		 *
		 * \throws NotPrimeException If the current \NTL modulus, as obtained by
		 * \code
		 * 	T::GFp::modulus();
		 * \endcode
		 * is not a prime (unless \a T is GF2_Algebra).
		 * \throws NotIrreducibleException If \a P is not an irreducible polynomial.
		 *
		 * \warning If \a P has degree \a d divisible by the characteristic and if its
		 * (\a d - 1)-th coefficient is 0, then you won't be able to construct Artin-Schreier
		 * extensions through a call to ArtinSchreierExtension(). In fact in this case the generator
		 * of the field has trace 0 and the costruction of [\ref ISSAC "DFS '09"] doesn't work.
		 */
		static const Field<T>& createField(const GFpX& P, const bool test = true)
		throw (NotPrimeException, NotIrreducibleException);
		/**
		 * \brief Build the field F<sub>p<sup>d</sup></sub> using a default
		 * polynomial.
		 *
		 * Notice that this operation implicitely creates
		 * the field F<sub>p</sub> too.
		 *
		 * \param [in] p A prime number.
		 * \param [in] d A positive integer.
		 * \param [in] test If false, do not perform a primality test on \a p.
		 * \return A reference to the newly created Field object.
		 *
		 * \throws NotPrimeException If \a p is not prime
		 * \throws BadParametersException If \a d is less than one.
		 * \throws BadParametersException If \a T is GF2_Algebra and \a p is different from 2.
		 *
		 * \todo For the moment the default polynomial is generated randomly. When \a p divides \a d ,
		 * this prevents from
		 * building successive Artin-Schreier extensions if the (\a d - 1)-th
		 * coeffcient of the randomly generated polynomial is 0.  In fact in this case the generator
		 * of the field has trace 0 and the costruction of [\ref ISSAC "DFS '09"] doesn't work.
		 * It would be nice to have a better generation routine that discards such polynomials.
		 */
		static const Field<T>& createField
		(const BigInt& p, const long d = 1, const bool test = true)
		throw (NotPrimeException, BadParametersException);
	/** @} */

	/****************//** \name Artin-Schreier Extensions
	 * Build Artin-Schreier extensions as described in [\ref ISSAC "DFS '09"].
	 * @{
	 */
		/**
		 * \brief Build a primitive extension of degree \a p
		 * as in [\ref ISSAC "DFS '09", Section 3].
		 *
		 * \return A reference to the newly created Field object.
		 *
		 * \throws CharacteristicTooLargeException If \ref p is a multiprecision
		 * integer larger than the largest single precision integer.
		 * \throws NotSupportedException If \ref d is divisible by \ref p
		 * and the generator of this field has trace 0. See creteField() and [\ref ISSAC "DFS '09"]
		 * for necessary and sufficient conditions for this not to hold.
		 */
		const Field<T>& ArtinSchreierExtension() const
		throw (CharacteristicTooLargeException, NotSupportedException);
		/**
		 * \brief Build the splitting field of the polynomial \f$ X^p - X - \mathtt{alpha} \f$
		 * as in [\ref ISSAC "DFS '09", Section 6].
		 *
		 * The polynomial must be reducible (equivalently \a alpha must have trace 0).
		 * A primitive Artin-Schreier extension
		 * is constructed using ArtinSchreierExtension(), then Couveignes2000() is used in such
		 * an extension to find one
		 * of the roots and the isomorphic field is constructed. The computed root is used
		 * in pushDown() and liftUp() to navigate the tower as in
		 * [\ref ISSAC "DFS '09", Section 6.2].
		 *
		 * \return A reference to the newly created Field object.
		 *
		 * \see ArtinSchreierExtension(), Couveignes2000(), pushDown(), liftUp().
		 *
		 * \throws CharacteristicTooLargeException If \ref p is a multiprecision
		 * integer larger than the largest single precision integer.
		 * \throws NotIrreducibleException If \a alpha has trace 0.
		 * \throws NotSupportedException If the primitive field cannot be created. See
		 * ArtinSchreierExtension() for details
		 * \throws IllegalCoercionException If \a alpha cannot be coerced to an element of this field.
		 */
		const Field<T>& ArtinSchreierExtension(const FieldElement<T>& alpha)
		const throw (CharacteristicTooLargeException, NotIrreducibleException,
				NotSupportedException, IllegalCoercionException);
		/**
		 * \brief Finds a root of the polynomial \f$ X^p - X - \mathtt{alpha} \f$.
		 *
		 * It uses the algorithm described in [\ref Cou00 "Couveignes '00"] and
		 * [\ref ISSAC "DFS '09", Section 6.1].
		 *
		 * \throws IllegalCoercionException if \a alpha cannot be
		 * 		coerced to this field.
		 * \throws IsIrreducibleException If the polynomial is irreducible (equivalently, if
		 * 		\a alpha has trace 0).
		 */
		FieldElement<T> Couveignes2000(const FieldElement<T>& alpha)
		const throw(IllegalCoercionException, IsIrreducibleException);
	/** @} */

	/****************//** \name Properties ******************/
	/** @{ */
		/**
		 * \brief \copybrief p
		 *
		 * \return The same value as \ref p.
		 */
		BigInt characteristic() const throw () { return p; }
		/**
		 * \brief \copybrief d
		 *
		 * \return The same value as \ref d.
		 */
		long degree() const throw () { return d; }
		/** \brief The cardinality of the field. */
		ZZ cardinality() const throw ();
		/**
		 * \brief \copybrief height
		 *
		 * \return The same value as \ref height.
		 */
		long ArtinSchreierHeight() const throw () { return height; }
		/**
		 * \brief The polynomial with coefficients in subField() that has been
		 * used to generate this extension.
		 *
		 * Or X - 1 if this is a prime field.
		 *
		 * \see [\ref ISSAC "DFS '09"].
		 *
		 * \invariant This is the minimal polynomial of generator().
		 */
		FieldPolynomial<T> generatingPolynomial() const throw();
		/**
		 * \brief The polynomial with coefficients in F<sub>p</sub> used to
		 * represent elements of this field.
		 *
		 * Or X - 1 if this is a prime field.
		 *
		 * \see [\ref ISSAC "DFS '09"].
		 *
		 * \invariant This is the minimal polynomial of primitiveElement().
		 */
		FieldPolynomial<T> primitivePolynomial() const throw();
	/** @} */

	/****************//** \name Field Elements
	 * Routines to create elements of the field.
	 * @{
	 */
		/**
		 * \brief The element \a i mod \a p
		 *
		 * \param [in] i An integer.
		 * \return An element of this field.
		 */
		FieldElement<T> scalar(const BigInt& i) const throw ();
		/** \brief The zero element of this field */
		FieldElement<T> zero() const throw () { return scalar(0); }
		/** \brief The identity element of this field */
		FieldElement<T> one() const throw () { return scalar(1); }
		/**
		 * \brief The generator over subField().
		 *
		 * \return A root of generatingPolynomial().
		 */
		FieldElement<T> generator() const throw () { return *gen; }
		/**
		 * \brief The generator over F<sub>p</sub>
		 *
		 * \return A root of primitivePolynomial().
		 *
		 * \invariant This is the same as generator() for primitive fields built using
		 * ArtinSchreierExtension().
		 */
		FieldElement<T> primitiveElement() const throw ()
		{ return *stem->primitive >> *this; }
		/** \brief A random element of the field. */
		FieldElement<T> random() const throw ();
	/** @} */

	/****************//** \name Access to the Infrastructure
	 * These methods let you use the internal \NTL representation of elements
	 * to build a FieldElement or a FieldPolynomial.
	 *
	 * \warning These methods are for advanced use only. Use them
	 * if you want to use an algorithm by you or in NTL that is not available for
	 * FieldElement or FieldPolynomial.
	 * @{
	 */
		/**
		 * \brief Set the current context to this field's context
		 *
		 * \NTL's context holds information about the modulus and the characteristic
		 * of a finite field. By calling this method you set the current \NTL context
		 * to this field's context, so that any subsequent operation on \NTL types
		 * such as T::GFpX or T::GFpE will use that context.
		 *
		 * You usually shouldn't be concerned about this method as the library takes
		 * care of switching the context for you when needed. The only time when you
		 * have to explicitly call it is when you want to use NTL types outside
		 * of the library and then transform the result to a FieldElement or
		 * FieldPolynomial through a call to a fromInfrastructure() method.
		 *
		 * \warning Be aware that the current context is undefined after any
		 * call to a function of this library.
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink
		 */
		void switchContext() const throw();
		/**
		 * \brief Build an element of this field from an \NTL type.
		 *
		 * Returns a new element of this field having \a e as representation.
		 *
		 * \param [in] e An \NTL scalar element.
		 * \return The newly created element.
		 * \warning \a e must have been created in the context of this field. See switchContext().
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink
		 */
		FieldElement<T> fromInfrastructure(const GFp& e) const throw();
		/**
		 * \brief Build an element of this field from an \NTL type.
		 *
		 * Returns a new element of this field having \a e as representation.
		 *
		 * \param [in] e an \NTL element.
		 * \return The newly created element.
		 * \throws IllegalCoercionException If this field is a prime field only scalar elements
		 * can belong to it. Use fromInfrastructure(const GFp&) const instead.
		 * \warning \a e must have been created in the context of this field. See switchContext().
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink
		 */
		FieldElement<T> fromInfrastructure(const GFpE& e) const throw(IllegalCoercionException);
		/**
		 * \brief Build an polynomial with coefficients in this field from an \NTL type.
		 *
		 * Returns a new polynomial over this field having \a P as representation.
		 *
		 * \param [in] P an \NTL scalar polynomial.
		 * \return The newly created element.
		 * \warning \a P must have been created in the context of this field. See switchContext().
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink
		 */
		FieldPolynomial<T> fromInfrastructure(const GFpX& P) const throw();
		/**
		 * \brief Build an polynomial with coefficients in this field from an \NTL type.
		 *
		 * Returns a new polynomial over this field having \a P as representation.
		 *
		 * \param [in] P an \NTL polynomial.
		 * \return The newly created element.
		 * \throws IllegalCoercionException If this field is a prime field only scalar polynomials
		 * can belong to it. Use fromInfrastructure(const GFpX&) const instead.
		 * \warning \a P must have been created in the context of this field. See switchContext().
		 * \see \link using_infrastructure.c++ using_infrastructure.c++ \endlink
		 */
		FieldPolynomial<T> fromInfrastructure(const GFpEX& P) const throw(IllegalCoercionException);
	/** @} */

	/****************//** \name Field lattice navigation
	 * These routines permit to move around in the lattice of fields
	 * created by calls to ArtinSchreierExtension() as described in
	 * section \ref Field_lattices.
	 * @{
	 */
		/** \brief The prime field F<sub>p</sub> of this field. */
		const Field<T>& primeField() const throw();
		/**
		 * \brief The base field of the Artin-Schreier tower.
		 *
		 * The Artin-Schreier height 0 subfield of this field.
		 * \return A reference to the base field.
		 * \see ArtinSchreierExtension(), height, ArtinSchreierHeight(), \ref Field_lattices.
		 */
		const Field<T>& baseField() const throw();
		/**
		 * \brief The immediate subfield in the primitive tower.
		 *
		 * If this field belongs to the primitive tower (the stem),
		 * then its immediate subfield is returned. Otherwise
		 * stemField().subField() is returned.
		 *
		 * \return A reference to the subfield.
		 * \throws NoSubFieldException If this is a prime field.
		 * \see ArtinSchreierExtension(), \ref Field_lattices.
		 */
		const Field<T>& subField() const throw(NoSubFieldException) {
#ifdef FAAST_DEBUG
			if (!stem) throw FAASTException("No stem in stemField().");
#endif
			if (!stem->subfield) throw NoSubFieldException();
			return *(stem->subfield);
		}
		/**
		 * \brief The immediate overfield in the primitive tower.
		 *
		 * If this field belongs to the primitive tower (the stem),
		 * then its immediate overfield is returned. Otherwise
		 * stemField().overField() is returned.
		 *
		 * \return A reference to the overfield.
		 * \throws NoOverFieldException If no overfield as been constructed
		 * through a call to ArtinSchreierExtension()
		 * \see ArtinSchreierExtension(), \ref Field_lattices.
		 */
		const Field<T>& overField() const throw(NoOverFieldException) {
#ifdef FAAST_DEBUG
			if (!stem) throw FAASTException("No stem in stemField().");
#endif
			if (!stem->overfield) throw NoOverFieldException();
			return *(stem->overfield);
		}
		/**
		 * \brief The field on the primitive tower isomorphic to this
		 * field.
		 *
		 * If this field belongs to the primitive tower, then it returns itself.
		 * Otherwise it returns the field in the primitive tower (the stem) that
		 * is isomorphic to this field. The isomorphism has been computed through
		 * Couveignes2000(const FieldElement<T>&) const as described in
		 * ArtinSchreierExtension(const FieldElement<T>&) const.
		 *
		 * \return A reference to the stem field.
		 * \see ArtinSchreierExtension(), Couveignes2000, \ref Field_lattices
		 */
		const Field<T>& stemField() const throw() {
#ifdef FAAST_DEBUG
			if (!stem) throw FAASTException("No stem in stemField().");
#endif
			return *stem;
		}
	/** @} */


	/****************//** \name Applying the isomorphism
	 * These routines implement the algorithms of [\ref ISSAC "DFS '09", Section 6.2]
	 * that convert elements written on the univariate basis of the primitive Artin-Schreier
	 * tower to and from the multivariate basis of any isomorphic tower.
	 * @{
	 */
		/**
		 * \brief Convert the univariate representation of \a e to the
		 * multivariate representation over this field.
		 *
		 * If this is a prime field, then \a v is filled with the
		 * coefficients in F<sub>p</sub> of its univariate representation. Otherwise let
		 * \code
		 * x = e.parent().generator();
		 * \endcode
		 * This method fills the vector \a v with \ref p elements of this field such that
		 * \f{equation}{
		 * \mathtt{e} = \mathtt{v[0]} + \mathtt{v[1]}*\mathtt{x} + ...
		 * 		+ \mathtt{v[p-1]}*\mathtt{x}^{p-1}
		 * 		\mathrm{.}
		 * \f}
		 *
		 * Let \b K be this field, this corresponds to convert \a e from its internal (univariate)
		 * representation to the bivariate representation as an element of \b K[\a x].
		 * A repeated application of this method implements \c ApplyInverse of
		 * [\ref ISSAC "DFS '09", Section 6.2].
		 *
		 * \param [in] e An element of any field isomorphic to overField().
		 * \param [out] v A vector of elements of this field that satisfies condition (1).
		 * \throw IllegalCoercionException If the field \a e belongs to
		 *         is not isomorphic to overField().
		 * \invariant When \a e belongs to a field in the primitive tower (the stem),
		 * this is equivalent to
		 * \link FieldElement::pushDown() \c pushDown(e, v) \endlink and then coerce all the contents of \a v to this field.
		 * \see pushDown().
		 */
		void toBivariate(const FieldElement<T>& e, vector<FieldElement<T> >& v) const
			 throw(IllegalCoercionException);

		/**
		 * \brief Convert the multivariate representation of \a v to the
		 * univariate representation of this field.
		 *
		 * If the elements of \a v belong to the prime field, then \a e is the element whose univariate
		 * representation has \a v as cofficients. Otherwise let
		 * \code
		 * x = generator();
		 * \endcode
		 * This method stores in \a e an element of this field such that
		 * \f{equation}{
		 * \mathtt{e} = \mathtt{v[0]} + \mathtt{v[1]}*\mathtt{x} + ...
		 * 		+ \mathtt{v[p-1]}*\mathtt{x}^{p-1}
		 * 		\mathrm{,}
		 * \f}
		 * If \a v is too short, it is filled with zeros. If \a v is too
		 * long, the unnecessary elements are ignored.
		 *
		 * Let \b K be the field containing the elements of \a v, this corresponds to convert
		 * \a v from the multivariate
		 * representation as an element of \b K[\a x] to the internal (univariate) representation
		 * of this field.
		 * A repeated application of this method implements
		 * \c ApplyIsomorphism of [\ref ISSAC "DFS '09", Section 6.2].
		 *
		 * \param [in] v A vector of elements all belonging to a field isomorphic to subField().
		 * \param [out] e An element of this field satisfying condition (2).
		 * \throw NotInSameFieldException If the elements of \a v do not
		 *         belong all to the same field.
		 * \throw IllegalCoercionException If the field the elements of \a v belong to
		 *         is not isomorphic to subField().
		 * \invariant When this field is in the primitive tower (the stem), this is equivalent to
		 * coerce all the contents of \a v to the stem and then
		 * \link FieldElement::liftUp() \c liftUp(v, e) \endlink.
		 * \see liftUp().
		 */
		void toUnivariate(const vector<FieldElement<T> >& v, FieldElement<T>& e) const
			throw(NotInSameFieldException, IllegalCoercionException);
	/** @} */

	/****************//** \name Predicates ******************/
	/** @{ */
		/** \brief Equality.
		 *
		 * \note Comparison on the address of the object.
		 */
		bool operator==(const Field<T>& F) const throw () { return this==&F; }
		/** \brief Inequality.
		 *
		 * \note Comparison on the address of the object.
		 */
		bool operator!=(const Field<T>& F) const throw () { return !(*this==F); }
		/**
		 * \brief This field is isomorphic to F.
		 * \return \c true only if the isomorphism has been computed.
		 * \see stemField().
		 */
		bool isIsomorphic(const Field<T>& F) const throw () { return stem==F.stem; }
		/**
		 * \brief This field is contained in \a F.
		 * \return \c true only if the inclusion has been computed.
		 * \see subField().
		 */
		bool isSubFieldOf(const Field<T>& F) const throw ();
		/**
		 * \brief This field contains \a F
		 * \return \c true only if the inclusion has been computed.
		 * \see overField()
		 */
		bool isOverFieldOf(const Field<T>& F) const throw ()
		{ return isIsomorphic(F) || F.isSubFieldOf(*this); }
		/**
		 * \brief This is a prime field.
		 */
		bool isPrimeField() const throw ()
		{ return stem->subfield == NULL; }
		/** This is the base field of an ArtinSchreier tower */
		bool isBaseField() const throw ()
		{ return !stem->overfield || stem->overfield->height == 1; }
	/** @} */
	/****************//** \name Printing ******************/
	/** @{ */
		/** \brief Print details about the field to \a o */
		ostream& print(ostream& o) const;
	/** @} */
	/****************** Destructor ******************/
		~Field() throw (FAASTException)
		{ throw FAASTException("Destroying fields is no good."); }


	/*****************************************************/
	/****************** Private section ******************/
	/*****************************************************/
	/** \cond DEV */
	private:
	/****************//** \name Access to precomputed values
	 * These methods grant acces to precomputed values. They compute the values
	 * on-demand when called the first time.
	 * @{ */
		const FieldElement<T>& getPseudotrace(const long i) const;
		const FieldElement<T>& getLiftup() const;
		const MatGFp& getArtinMatrix() const;
		const Context& getCyclotomic() const;
	/** @} */

	/****************** Copy prohibited ******************/
		void operator=(const Field<T>&);
		Field(const Field<T>&);

	/****************** Couveignes 2000 subroutines ******************/
		/**
		 * \brief Couveignes' algorithm as in [\ref ISSAC "DFS '09", Section 6.1].
		 *
		 * This routine is internally called by Couveignes2000().
		 *
		 * \pre Assumes Tr(\a alpha) = 0.
	 	 */
		void couveignes00(FieldElement<T>& res, const FieldElement<T>& alpha) const;

	/****************//** \name Constructors ******************/
	/** @{ */
		/** \brief Construct a field with specified parameters */
		Field<T> (
			const Field<T>* sub,
			const Field<T>* over,
			const Context& ctxt,
			const FieldElement<T>* pri,
			const vector<FieldElement<T> >& pseudo,
			const FieldElement<T>* liftup,
			const MatGFp& mat,
			const long line,
			const bool pluso,
			const bool twopminuso,
			const Context* Ph,
			const Field<T>* st,
			const Field<T>* vsub,
			const FieldElement<T>* g,
			const FieldElement<T>* a,
			const BigInt& cha,
			const long deg,
			const long h
		) throw() :
		subfield(sub), overfield(over),
		context(ctxt),
		primitive(pri),
		pseudotraces(pseudo),
		liftuphelper(liftup),
		artin(mat), artinLine(line),
		plusone(pluso), twopminusone(twopminuso),
		Phi(Ph),
		stem(st), vsubfield(vsub),
		gen(g), alpha(a),
		p(cha), d(deg), height(h)
		{}
		/** \brief Construct a base fields */
		Field<T> (
			const Field<T>* sub,
			const Context& ctxt,
			const GFpE& pri,
			const BigInt& cha,
			const long deg,
			const GFpE& g
		) throw() :
		subfield(sub), overfield(NULL),
		context(ctxt),
		primitive(new FieldElement<T>(this, pri)),
		pseudotraces(),
		liftuphelper(),
		artin(), artinLine(-1),
		plusone(false), twopminusone(false),
		Phi(),
		stem(this), vsubfield(NULL),
		gen(new FieldElement<T>(this, g)),
		alpha(),
		p(cha), d(deg), height(0)
		{}
		/** \brief Construct F<sub>p</sub> */
		Field<T> (
			const Context& ctxt,
			const GFp& pri,
			const BigInt& cha
		) throw() :
		subfield(NULL), overfield(NULL),
		context(ctxt),
		primitive(new FieldElement<T>(this, pri)),
		pseudotraces(),
		liftuphelper(),
		artin(), artinLine(-1),
		plusone(false), twopminusone(false),
		Phi(),
		stem(this), vsubfield(NULL),
		gen(new FieldElement<T>(this, pri)),
		alpha(),
		p(cha), d(1), height(0)
		{}
		/** \brief Construct a field of the primitive tower (the stem). */
		Field<T> (
			const Field<T>* sub,
			const Context& ctxt,
			const GFpE& pri,
			const bool po,
			const bool tpmo,
			const BigInt& cha,
			const long deg,
			const long h,
			const FieldElement<T>* aleph,
			const Field<T>* vsub = NULL
		) throw() :
		subfield(sub), overfield(NULL),
		context(ctxt),
		primitive(new FieldElement<T>(this, pri)),
		pseudotraces(),
		liftuphelper(),
		artin(), artinLine(-1),
		plusone(po), twopminusone(tpmo),
		Phi(),
		stem(this), vsubfield(vsub),
		gen(new FieldElement<T>(this, pri)),
		alpha(aleph),
		p(cha), d(deg), height(h)
		{}
		/** \brief Construct a generic field (not on the stem) */
		Field<T> (
			const Field<T>* st,
			const FieldElement<T>& gen,
			const FieldElement<T>* aleph,
			const Field<T>* vsub
		) throw() :
		subfield(NULL), overfield(NULL),
		context(),
		primitive(NULL),
		pseudotraces(),
		liftuphelper(),
		artin(), artinLine(-1),
		plusone(), twopminusone(),
		Phi(),
		stem(st), vsubfield(vsub),
		gen(new FieldElement<T>(this, gen.repBase, gen.repExt, gen.base)),
		alpha(aleph),
		p(st->p), d(st->d), height(st->height)
		{}
	/** @} */
	/** \endcond */
	};

	/****************** Printing ******************/
	/** \brief Print details about \a F to \a o
	 * \relates Field
	 */
	template <class T> ostream& operator<<(ostream& o, const Field<T>& F) {
		return F.print(o);
	}
}

#endif /*FIELD_H_*/
