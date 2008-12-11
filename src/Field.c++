#include "Types.hpp"

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2XFactoring.h>

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
	/* Build the field GF(p^d) using some default
	 * polynomial.
	 * Notice that this operation implicitely creates
	 * the field GF(p) too.
	 * 
	 * throws : NotPrimeException if p is not prime
	 * throws : ASException if d less than one
	 */
	template<> Field<ZZ_p_Algebra>::Field(long cha, long deg) throw (ASException) :
	subfield(NULL), overfield(NULL), vsubfield(NULL), stem(NULL), 
	plusone(false), twopplusone(false), 
	p(cha), d(deg), height(0) {
		if (d < 1) {
			throw ASException("Cannot create an extension field of negative degree.");
		} else if (d >= 2) {
			// create GF(p)
			subfield = new Field<ZZ_p_Algebra>(p);
			subfield->overfield = this;
			// generate an irreducible polynomial and save the context
			context.p.save();
			GFpX P; BuildIrred(P, d);
			GFpE::init(P);
			context.P.save();
			// compute the inverse matrix of X^p-X
			GFpXModulus Pmod = GFpE::modulus();
			artin = artinMatrix<ZZ_p_Algebra>(p,d,Pmod);
		} else {
			// if modulus isn't a prime, stop
			if (p <= 1 || !ProbPrime(p)) {
				throw NotPrimeException();
			}
			GFp::init(to_ZZ(p));
			context.p.save();
		}
		// compute the primitive element
		primitive = FieldElement<ZZ_p_Algebra>(this, GFp(to_ZZ_p(1)), to_ZZ_pE(GFpX(1,1)), d==1);
		gen = primitive;
	}

	template<> Field<GF2_Algebra>::Field(long cha, long deg) throw (ASException) :
	subfield(NULL), overfield(NULL), vsubfield(NULL), stem(NULL), 
	plusone(false), twopplusone(false), 
	p(2), d(deg), height(0) {
		if (d < 1) {
			throw ASException("Cannot create an extension field of negative degree.");
		} else if (d >= 2) {
			// create GF(p)
			subfield = new Field<GF2_Algebra>(p);
			subfield->overfield = this;
			// generate an irreducible polynomial and save the context
			GFpX P; BuildIrred(P, d);
			GFpE::init(P);
			context.P.save();
			// compute the inverse matrix of X^p-X
			GFpXModulus Pmod = GFpE::modulus();
			artin = artinMatrix<GF2_Algebra>(p,d,Pmod);
		}
		// else create GF(p)  (nothing to do)
		// compute the primitive element
		primitive = FieldElement<GF2_Algebra>(this, GFp(to_GF2(1)), to_GF2E(GFpX(1,1)), d==1); 
		gen = primitive;
	}

	/* Default constructor, builds a field from some default value
	 * (i.e. the context for NTL).
	 * 
	 * throws : NotPrimeException, NotIrreducibleException if the
	 *          default makes no sense
	 */
	template<> Field<ZZ_p_Algebra>::Field() throw (NotPrimeException, NotIrreducibleException) :
	subfield(NULL), overfield(NULL), vsubfield(NULL), stem(NULL), 
	plusone(false), twopplusone(false), 
	p(to_long(GFp::modulus())), d(max(1,GFpE::degree())), height(0) {
		// if there's a modulus that makes sense
		if (d >= 2) {
			// create GF(p)
			subfield = new Field<ZZ_p_Algebra>(p);
			subfield->overfield = this;
			// if the modulus isn't irreducible
			if (!DetIrredTest(GFpE::modulus())) {
				throw NotIrreducibleException();
			}
			// store moduli
			context.p.save();
			context.P.save();
			// compute the inverse matrix of X^p-X
			GFpXModulus P = GFpE::modulus();
			artin = artinMatrix<ZZ_p_Algebra>(p,d,P);
		}
		// else create GF(p)
		else {
			// if modulus isn't a prime, stop
			if (p <= 1 || !ProbPrime(p)) {
				throw NotPrimeException();
			}
			context.p.save();
		}
		// compute the primitive element
		primitive = FieldElement<ZZ_p_Algebra>(this, GFp(to_ZZ_p(1)), to_ZZ_pE(GFpX(1,1)), d==1); 
		gen = primitive;
	}

	template<> Field<GF2_Algebra>::Field() throw (NotPrimeException, NotIrreducibleException) :
	subfield(NULL), overfield(NULL), vsubfield(NULL), stem(NULL), 
	plusone(false), twopplusone(false), 
	p(2), d(max(1,GFpE::degree())), height(0) {
		// if there's a modulus that makes sense
		if (d >= 2) {
			// create GF(p)
			subfield = new Field<GF2_Algebra>(p);
			subfield->overfield = this;
			// if the modulus isn't irreducible
			if (!IterIrredTest(GFpE::modulus())) {
				throw NotIrreducibleException();
			}
			// store moduli
			context.P.save();
			// compute the inverse matrix of X^p-X
			GFpXModulus P = GFpE::modulus();
			artin = artinMatrix<GF2_Algebra>(p,d,P);
		}
		// else create GF(p)  (nothing to do)
		// compute the primitive element
		primitive = FieldElement<GF2_Algebra>(this, GFp(to_GF2(1)), to_GF2E(GFpX(1,1)), d==1); 
		gen = primitive;
	}

	/* Build a field from an irreducible polynomial P.
	 * 
	 * throws : NotIrreducibleException if P is not irreducible
	 */
//	Field(const Poly& P) throw (NotIrreducibleException);

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
//	long characteristic() const throw ();
//	long degree() const throw ();
//	BigInt cardinality() const throw ();
	/* The number of Artin-Schreier extensions from the base
	 * field.
	 */
//	long ArtinSchreierHeight() const throw ();
	
/****************** Field Elements ******************/
	/* Constructs the element i times 1 */
//	FieldElement<T> scalar(const long i) const throw ();
//	FieldElement<T> zero() const throw () { return scalar(0); }
//	FieldElement<T> one() const throw () { return scalar(1); }
	/* The generator over the immediately preceding subfield */
//	FieldElement<T> generator() const throw ();
	/* The generator over GF(p) */
//	FieldElement<T> primitiveElement() const throw ();
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
	/* Two fields are the same only if they are the same
	 * object.
	 */
//	bool operator==(const Field<T>& f) const throw ();
//	bool operator!=(const Field<T>& f) const throw () { return !this==f; }
	/* Two fields are isomorphic if the isomorphism between them
	 * has actually been computed
	 */
//	bool isIsomporphic(const Field<T>& f) const throw ();
	/* There is inclusion between two fields only if the inclusion
	 * has actually been computed.
	 */
//	bool isSubFieldOf(const Field<T>& f) const throw ();
//	bool isOverFieldOf(const Field<T>& f) const throw ();
/****************** Printing ******************/
//	ostream& print(ostream&) const;
/****************** Destructor ******************/
//	~Field() throw ();


/*****************************************************/
/****************** Private section ******************/
/*****************************************************/

/****************** Copy prohibited ******************/
//	void operator=(const Field<T>&);
//	Field(const Field<T>&);
	
/****************** Internal Constructors ******************/
	/* For constructing field extensions */
	//Field<T> () throw();
}
