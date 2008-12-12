#include "Types.hpp"

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2XFactoring.h>
#include <string>
#include <sstream>

namespace AS {
	
	template <class T> typename T::MatGFp artinMatrix(const long p, const long d,
	const typename T::GFpXModulus& P) {
		typedef typename T::MatGFp      MatGFp;
		typedef typename T::GFpX        GFpX;

		// we just want the d-1 minor on the bottom right
		MatGFp appl; appl.SetDims(d-1, d-1);
		GFpX Xp = PowerXMod(p, P);
		GFpX column(0,1);
		for (long i = 1 ; i < d ; i++) {
			MulMod(column, column, Xp, P);
			for (int j = 1 ; j < d ; j++) {
				appl[j-1][i-1] = coeff(column, j);
			}
			appl[i-1][i-1]--;
		}
		//invert
		return inv(appl);		
	}

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
	template<> Field<ZZ_p_Algebra>& Field<ZZ_p_Algebra>::createField
	(bool test)	throw (NotPrimeException, NotIrreducibleException) {
		// retrieve all the informations on the modulus
		long p = to_long(GFp::modulus());
		GFpXModulus P = GFpE::modulus();
		long d = deg(P);
		Context context; context.p.save(); context.P.save();
		// test primality
		if ( p <= 1 || (test && !ProbPrime(p)) ) {
			throw NotPrimeException();
		}
		// test irreducibility
		if (d > 1 && test && !DetIrredTest(P)) {
			throw NotIrreducibleException();
		}
		
		// build GF(p^d)
		if (d >= 2) {
			// build GF(p)
			GFp one; one = 1;
			Field<ZZ_p_Algebra>* baseField = new Field<ZZ_p_Algebra>(context, one, p);
			// compute the generator of this field
			GFpE primitive; conv(primitive, GFpX(1,1));
			// compute the inverse matrix of X^p-X
			GFpXModulus Pmod; build(Pmod, P);
			MatGFp artin = artinMatrix<ZZ_p_Algebra>(p,d,Pmod);
			// build the field
			Field<ZZ_p_Algebra>* K =
				new Field<ZZ_p_Algebra>(baseField, context, primitive,
				                        artin, p, d, primitive);
			// connect the base field
			baseField->overfield = K;
			return *K;
		}
		// build GF(p)
		else {
			GFp primitive; primitive = 1;
			return *(new Field<ZZ_p_Algebra>(context, primitive, p));
		}
	}

	template<> Field<GF2_Algebra>& Field<GF2_Algebra>::createField
	(bool test)	throw (NotPrimeException, NotIrreducibleException) {
		// retrieve all the informations on the modulus
		long p = 2;
		GFpXModulus P = GFpE::modulus();
		long d = deg(P);
		Context context; context.P.save();
		// test irreducibility
		if (d > 1 && test && !IterIrredTest(P)) {
			throw NotIrreducibleException();
		}
		
		// build GF(p^d)
		if (d >= 2) {
			// build GF(p)
			GFp one; one = 1;
			Field<GF2_Algebra>* baseField = new Field<GF2_Algebra>(context, one, p);
			// compute the generator of this field
			GFpE primitive; conv(primitive, GFpX(1,1));
			// compute the inverse matrix of X^p-X
			GFpXModulus Pmod; build(Pmod, P);
			MatGFp artin = artinMatrix<GF2_Algebra>(p,d,Pmod);
			// build the field
			Field<GF2_Algebra>* K =
				new Field<GF2_Algebra>(baseField, context, primitive,
				                        artin, p, d, primitive);
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
	template<class T> Field<T>& Field<T>::createField
	(const GFpX& P, bool test)
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
	template<> Field<ZZ_p_Algebra>& Field<ZZ_p_Algebra>::createField
	(long p, long d, bool test)
	throw (NotPrimeException, BadParametersException) {
		if (d < 1) {
			throw BadParametersException("Cannot create an extension field of negative degree.");
		}
		if (p <= 1) {
			throw NotPrimeException();
		}
		if (test && !ProbPrime(p)) {
			throw NotPrimeException();
		}
		
		GFp::init(to_ZZ(p));
		GFpX P;
		if (d >= 2) BuildIrred(P, d);
		else SetX(P);
		return createField(P, false);
	}
	
	template<> Field<GF2_Algebra>& Field<GF2_Algebra>::createField
	(long p, long d, bool test)
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
		if (d >= 2) BuildIrred(P, d);
		else SetX(P);
		return createField(P, false);
	}

/****************** Field Extensions ******************/
	/* Build a default extension of degree p over this field. 
	 */
//	Field<T> ArtinSchreierExtension() const throw ();
	/* Build the splitting field of the polynomial
	 * 			X^p - X - alpha
	 * over this field. This may or may not be an extension of
	 * degree p depending if the polynomial is irreducible.
	 */
//	Field<T> ArtinSchreierExtension(const FieldElement<T>& alpha) const throw ();
	
/****************** Properties ******************/
	template <class T> typename T::BigInt Field<T>::cardinality() const throw () {
		BigInt c;
		power(c, p, d);
		return c;
	}
	
/****************** Field Elements ******************/
	/* Constructs the element i times 1 */
//	FieldElement<T> scalar(const long i) const throw ();
	/* A random element of the field */
//	FieldElement<T> random() const throw ();
	/* Interface with infrastructure. Use this only if you are
	 * sure of what you do !
	 */
//	FieldElement<T> fromInfrastructure(const GFp&) const throw();
//	FieldElement<T> fromInfrastructure(const GFpE&) const throw(IllegalCoercionException);
//	FieldPolynomial<T> fromInfrastructure(const GFpX&) const throw();
//	FieldPolynomial<T> fromInfrastructure(const GFpEX&) const throw(IllegalCoercionException);

/****************** Field lattice navigation ******************/
	/* The field GF(p) */
//	Field<T> baseField() const throw();

/****************** Level embedding ******************/
	/* Push the element e down to this field and store
	 * the result in v.
	 * 
	 * throw : IllegalCoercionException if the field e belongs to
	 *         is not the immediate overfield of this.
	 */
//	void pushDown(FieldElement<T>& e, vector<FieldElement<T> >& v) const
//		 throw(IllegalCoercionException);

	/* Lift the elements in v up to this field and store the result in e.
	 * If v is too short, it is filled with zeros. If v is too
	 * long, the unnecessary elements are ignored.
	 * 
	 * throw : NotInSameFieldException if the elements of v do not
	 *         belong all to the same field.
	 * throw : IllegalCoercionException if the field e belongs to
	 *         is not the immediate subfield of this.
	 */
//	void liftUp(vector<FieldElement<T> >& v, FieldElement<T>& e) const
//		throw(NotInSameFieldException, IllegalCoercionException);

/****************** Comparison ******************/
	/* There is inclusion between two fields only if the inclusion
	 * has actually been computed.
	 */
//	bool isSubFieldOf(const Field<T>& f) const throw ();
//	bool isOverFieldOf(const Field<T>& f) const throw ();
/****************** Printing ******************/
//	ostream& print(ostream&) const;
}
