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
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/GF2XFactoring.h>
#include <string>
#include <sstream>

namespace FAAST {
/****************** Constructors ******************/
	/* All constructors are static. There's no way to directly
	 * create a Field object. Field objects are permanent and
	 * they cannot be deleted after creation. They live in their
	 * own lattice structure.
	 */

	/* Default constructor, builds a field from some default value
	 * (i.e. the context for NTL).
	 *
	 * If test is false, do not perform primality and irreducibility
	 * tests
	 *
	 * throws : NotPrimeException, NotIrreducibleException if the
	 *          default makes no sense
	 */
	template <class T> const Field<T>& Field<T>::createField(const bool test)
	throw (NotPrimeException, NotIrreducibleException) {
		// retrieve all the informations on the modulus
		BigInt p = GFp::modulus();
		const GFpXModulus& P = GFpE::modulus();
		long d = deg(P);
		Context context; context.p.save(); context.P.save();
		// test primality
#ifdef FAAST_TIMINGS
		TIME.PRIMETEST = -GetTime();
#endif
		if ( p <= long(1) || (test && !ProbPrime(p)) ) {
			throw NotPrimeException();
		}
#ifdef FAAST_TIMINGS
		TIME.PRIMETEST += GetTime();
#endif
		// test irreducibility
#ifdef FAAST_TIMINGS
		TIME.IRREDTEST = -GetTime();
#endif
		if (d > 1 && test && !DetIrredTest(P)) {
			throw NotIrreducibleException();
		}
#ifdef FAAST_TIMINGS
		TIME.IRREDTEST += GetTime();
#endif
		// build GF(p^d)
		if (d >= 2) {
			// build GF(p)
			GFp one; one = 1;
			Field<T>* baseField = new Field<T>(context, one, p);
			// compute the generator of this field
			GFpE primitive; conv(primitive, GFpX(1,1));
			// build the field
			Field<T>* K = new Field<T>(baseField, context, primitive,
			                           p, d, primitive);
			// connect the base field
			baseField->overfield = K;
			return *K;
		}
		// build GF(p)
		else {
			GFp primitive; primitive = 1;
			return *(new Field<T>(context, primitive, p));
		}
	}

	template<> const Field<GF2_Algebra>& Field<GF2_Algebra>::createField
	(const bool test) throw (NotPrimeException, NotIrreducibleException) {
		// retrieve all the informations on the modulus
		long p = 2;
		const GFpXModulus& P = GFpE::modulus();
		long d = deg(P);
		Context context; context.P.save();
		// test irreducibility
#ifdef FAAST_TIMINGS
		TIME.IRREDTEST = -GetTime();
#endif
		if (d > 1 && test && !IterIrredTest(P)) {
			throw NotIrreducibleException();
		}
#ifdef FAAST_TIMINGS
		TIME.IRREDTEST += GetTime();
#endif

		// build GF(p^d)
		if (d >= 2) {
			// build GF(p)
			GFp one; one = 1;
			Field<GF2_Algebra>* baseField = new Field<GF2_Algebra>(context, one, p);
			// compute the generator of this field
			GFpE primitive; conv(primitive, GFpX(1,1));
			// build the field
			Field<GF2_Algebra>* K =
				new Field<GF2_Algebra>(baseField, context, primitive,
				                        p, d, primitive);
			// connect the base field
			baseField->overfield = K;
			return *K;
		}
		// build GF(p)
		else {
			GFp primitive; primitive = 1;
			return *(new Field<GF2_Algebra>(context, primitive, p));
		}
	}


	/* Build a field from an irreducible polynomial P.
	 *
	 * If test is false, do not perform primality and irreducibility
	 * tests
	 *
	 * throws : NotPrimeExeption if the polynomial isn't defined
	 *          over a prime field
	 * throws : NotIrreducibleException if P is not irreducible
	 */
	template<class T> const Field<T>& Field<T>::createField
	(const GFpX& P, const bool test)
	throw (NotPrimeException, NotIrreducibleException) {
		GFpE::init(P);
		return createField(test);
	}

	/* Build the field GF(p^d) using some default
	 * polynomial.
	 * Notice that this operation implicitely creates
	 * the field GF(p) too.
	 *
	 * If test is false, do not perform primality and irreducibility
	 * tests
	 *
	 * throws : NotPrimeException if p is not prime
	 * throws : ASException if d less than one
	 */
	template <class T> const Field<T>& Field<T>::createField
	(const BigInt& p, const long d, const bool test)
	throw (NotPrimeException, BadParametersException) {
		if (d < 1) {
			throw BadParametersException("Cannot create an extension field of negative degree.");
		}
		if (p <= long(1)) {
			throw NotPrimeException();
		}
#ifdef FAAST_TIMINGS
		TIME.PRIMETEST = -GetTime();
#endif
		if (test && !ProbPrime(p)) {
			throw NotPrimeException();
		}
#ifdef FAAST_TIMINGS
		TIME.PRIMETEST += GetTime();
#endif

		GFp::init(p);
		GFpX P;
#ifdef FAAST_TIMINGS
		TIME.BUILDIRRED = -GetTime();
#endif
		if (d >= 2) BuildIrred(P, d);
		else SetX(P);
#ifdef FAAST_TIMINGS
		TIME.BUILDIRRED += GetTime();
#endif
		return createField(P, false);
	}

	template<> const Field<GF2_Algebra>& Field<GF2_Algebra>::createField
	(const BigInt& p, const long d, const bool test)
	throw (NotPrimeException, BadParametersException) {
		if (d < 1) {
			throw BadParametersException("Cannot create an extension field of negative degree.");
		}
		if (p != 2) {
			stringstream s;
			s << "GF2 does not support characteristic " << p << ".";
			throw BadParametersException(s.str().c_str());
		}

		GFpX P;
#ifdef FAAST_TIMINGS
		TIME.BUILDIRRED = -GetTime();
#endif
		if (d >= 2) BuildIrred(P, d);
		else SetX(P);
#ifdef FAAST_TIMINGS
		TIME.BUILDIRRED += GetTime();
#endif
		return createField(P, false);
	}

/****************** Properties ******************/
	template <class T> ZZ Field<T>::cardinality() const throw () {
		switchContext();
		if (d == 1) return to_ZZ(p);
		else return GFpE::cardinality();
	}

	/* The minimal polynomial over the immediate subfield of the
	 * element returned by generator.
	 */
	template <class T> FieldPolynomial<T> Field<T>::generatingPolynomial()
	const throw() {
		// if this is GF(p) or a base field, this is the same as
		// the primitive polynomial
		if (height == 0) return primitivePolynomial();
		// else, it is X^p - X - alpha
		else {
			switchContext();
			FieldPolynomial<T> res(-(*alpha));
			res.setCoeff(1, -1);
			res.setCoeff(to_long(p));
			return res;
		}
	}

	/* The minimal polynomial over GF(p) of the element returned
	 * by primitiveElement().
	 */
	template <class T> FieldPolynomial<T> Field<T>::primitivePolynomial()
	const throw() {
		switchContext();
		// if this is GF(p), return X-1
		if (d == 1) {
			FieldPolynomial<T> res = -one();
			res.setCoeff(1);
			return res;
		}
		// else, return the modulus
		else return FieldPolynomial<T>(this, GFpE::modulus().val());
	}

/****************** Field Elements ******************/
	/* Constructs the element i times 1 */
	template <class T> FieldElement<T> Field<T>::scalar(const BigInt& i)
	const throw () {
		switchContext();
		if (d == 1) {
			GFp x; x = i;
			return FieldElement<T>(this, x);
		} else {
			GFpE x; x = i;
			return FieldElement<T>(this, x);
		}
	}
	/* A random element of the field */
	template <class T> FieldElement<T> Field<T>::random() const throw () {
		switchContext();
		if (d == 1) {
			GFp x; NTL::random(x);
			return FieldElement<T>(this, x);
		} else {
			GFpE x; NTL::random(x);
			return FieldElement<T>(this, x);
		}
	}

/****************** Infrastructure ******************/
	/* Use these methods only if you are sure of what you do ! */

	/* Build elements from infrastracture */
	template <class T> FieldElement<T> Field<T>::fromInfrastructure(const GFp& x)
	const throw() {
		switchContext();
		if (d == 1) return FieldElement<T>(this, x);
		else {
			GFpE X; X = x;
			return FieldElement<T>(this, X);
		}
	}

	template <class T> FieldElement<T> Field<T>::fromInfrastructure(const GFpE& x)
	const throw(IllegalCoercionException) {
		if (d == 1) throw IllegalCoercionException();
		return FieldElement<T>(this, x);
	}

	template <class T> FieldPolynomial<T> Field<T>::fromInfrastructure(const GFpX& x)
	const throw() {
		switchContext();
		if (d == 1) return FieldPolynomial<T>(this, x);
		else {
			GFpEX X; conv(X, x);
			return FieldPolynomial<T>(this, X);
		}
	}

	template <class T> FieldPolynomial<T> Field<T>::fromInfrastructure(const GFpEX& x)
	const throw(IllegalCoercionException) {
		if (d == 1) throw IllegalCoercionException();
		return FieldPolynomial<T>(this, x);
	}

	/* Set the context to work in this field */
	template<class T> void Field<T>::switchContext() const throw() {
		stem->context.p.restore();
		stem->context.P.restore();
	}
	template<> void Field<GF2_Algebra>::switchContext() const throw() {
		stem->context.P.restore();
	}

/****************** Field lattice navigation ******************/
	/* The field GF(p) */
	template <class T> const Field<T>& Field<T>::primeField() const throw() {
		if (!stem->subfield) return *stem;
		const Field<T>* result = stem->subfield;
		while (result->subfield) result = result->subfield;
		return *result;
	}
	/* The base field of the tower */
	template <class T> const Field<T>& Field<T>::baseField() const throw() {
		const Field<T>* result = stem;
		while (result->height != 0) result = result->subfield;
		return *result;
	}

/****************** Comparison ******************/
	/* There is inclusion between two fields only if the inclusion
	 * has actually been computed.
	 */
	template <class T> bool Field<T>::isSubFieldOf(const Field<T>& f)
	const throw () {
		const Field<T>* stemf = f.stem;
#ifdef FAAST_DEBUG
		if (!stemf || !stem) throw ASException("Error : Stem is NULL.");
#endif
		while (stemf != stem && (stemf = stemf->subfield));
		return stemf == stem;
	}

/****************** Printing ******************/
	template <class T> ostream& Field<T>::print(ostream& o) const {
		o << "Finite field GF(" << p;
		if (d > 1) o << "^" << d;
		return o << ")";
	}
}
