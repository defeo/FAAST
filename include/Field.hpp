#ifndef FIELD_H_
#define FIELD_H_

#include "Exceptions.hpp"
#include "FieldElement.hpp"
#include "FieldPolynomial.hpp"
#include <memory>

namespace AS {
	template <class T> class Field {
	private:
		typedef typename T::GFp         GFp;
		typedef typename T::MatGFp      MatGFp;
		typedef typename T::GFpX        GFpX;
		typedef typename T::GFpE        GFpE;
		typedef typename T::GFpEX       GFpEX;
		typedef typename T::BigInt      BigInt;
		typedef typename T::Context     Context;
		typedef typename T::GFpXModulus GFpXModulus;
		
	/****************** Members for the stem ******************/
		/* The immediate subfield and overfield, if they're defined */
		const Field<T>* subfield;
		const Field<T>* overfield;
		/* Infrastructure-dependent data to perform computations in
		 * the field. (e.g. GF2EContext in NTL)
		 */
		Context context;
		/* The generator over GF(p) */
		const auto_ptr<const FieldElement<T> > primitive;
		/* Precomputed pseudotraces */
		const auto_ptr<const vector<FieldElement<T> > > pseudotraces;
		/* Lift-up precomputation */
		const auto_ptr<const FieldElement<T> > liftuphelper;
		/* The inverse matrix of the d-1 minor of the
		 * linear application X^p-X
		 */
		const MatGFp artin;
		/* The line we took away from X^p-X to make
		 * it invertible
		 */
		const long artinLine;
		/* Flags related to the construction of the extension */
		const bool plusone;
		const bool twopminusone;
		/* The (2p-1)th cyclotomic polynomial */
		mutable auto_ptr<const Context> Phi;

	/****************** Members for non-stem fields ******************/
		/* The stem-field to which this one is isomorphic */
		const Field<T>* stem;
		
	/****************** Common members ******************/
		/* The virtual subfield, if it is different from the real.
		 * It is given mainly for printing purposes
		 */
		const Field<T>* vsubfield;
		/* the characteristic */
		const BigInt p;
		/* the degree */
		const long d;
		/* The Artin-Schreier height */
		const long height;
		/* The generator over the subfield */
		const auto_ptr<const FieldElement<T> > gen;
		/* The element of the subfield such that this field is
		 * the Artin-Schreier extension defined by the polynomial
		 * 			X^p - X - alpha
		 */
		const auto_ptr<const FieldElement<T> > alpha;
		

	public:
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
		static const Field<T>& createField(const bool test = true)
		throw (NotPrimeException, NotIrreducibleException);
		/* Build a field from an irreducible polynomial P.
		 * 
		 * If test is false, do not perform primality and irreducibility
		 * tests
		 * 
		 * throws : NotPrimeExeption if the polynomial isn't defined
		 *          over a prime field
		 * throws : NotIrreducibleException if P is not irreducible
		 */
		static const Field<T>& createField(const GFpX& P, const bool test = true)
		throw (NotPrimeException, NotIrreducibleException);
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
		static const Field<T>& createField
		(const BigInt& p, const long d = 1, const bool test = true)
		throw (NotPrimeException, BadParametersException);
	
	/****************** Field Extensions ******************/
		/* Build a default extension of degree p over this field. 
		 */
		const Field<T>& ArtinSchreierExtension() const
		throw (CharacteristicTooLargeException, NotSupportedException);
		/* Build the splitting field of the polynomial
		 * 			X^p - X - alpha
		 * over this field. This may or may not be an extension of
		 * degree p depending if the polynomial is irreducible.
		 */
		const Field<T>& ArtinSchreierExtension(const FieldElement<T>& alpha)
		const throw (CharacteristicTooLargeException, NotSupportedException);
		
	/****************** Properties ******************/
		BigInt characteristic() const throw () { return p; }
		long degree() const throw () { return d; }
		ZZ cardinality() const throw ();
		/* The number of Artin-Schreier extensions from the base
		 * field.
		 */
		long ArtinSchreierHeight() const throw () { return height; }
		/* The minimal polynomial over the immediate subfield of the
		 * element returned by generator.
		 */
		FieldPolynomial<T> generatingPolynomial() const throw();
		/* The minimal polynomial over GF(p) of the element returned
		 * by primitiveElement().
		 */
		FieldPolynomial<T> primitivePolynomial() const throw();
		
	/****************** Field Elements ******************/
		/* Constructs the element i times 1 */
		FieldElement<T> scalar(const BigInt& i) const throw ();
		FieldElement<T> zero() const throw () { return scalar(0); }
		FieldElement<T> one() const throw () { return scalar(1); }
		/* The generator over the immediately preceding subfield */
		FieldElement<T> generator() const throw () { return *gen; }
		/* The generator over GF(p) */
		FieldElement<T> primitiveElement() const throw () { return *primitive; } 
		/* A random element of the field */
		FieldElement<T> random() const throw ();
		
	/****************** Infrastructure ******************/
		/* Use these methods only if you are sure of what you do ! */

		/* Build elements from infrastracture */
		FieldElement<T> fromInfrastructure(const GFp&) const throw();
		FieldElement<T> fromInfrastructure(const GFpE&) const throw(IllegalCoercionException);
		FieldPolynomial<T> fromInfrastructure(const GFpX&) const throw();
		FieldPolynomial<T> fromInfrastructure(const GFpEX&) const throw(IllegalCoercionException);
		/* Set the context to work in this field */
		void switchContext() const throw();
	
	/****************** Field lattice navigation ******************/
		/* The field GF(p) */
		const Field<T>& baseField() const throw();

	/****************** Level embedding ******************/
		/* Push the element e down to this field and store
		 * the result in v.
		 * 
		 * throw : IllegalCoercionException if the field e belongs to
		 *         is not the immediate overfield of this.
		 */
		void pushDown(FieldElement<T>& e, vector<FieldElement<T> >& v) const
			 throw(IllegalCoercionException);

		/* Lift the elements in v up to this field and store the result in e.
		 * If v is too short, it is filled with zeros. If v is too
		 * long, the unnecessary elements are ignored.
		 * 
		 * throw : NotInSameFieldException if the elements of v do not
		 *         belong all to the same field.
		 * throw : IllegalCoercionException if the field e belongs to
		 *         is not the immediate subfield of this.
		 */
		void liftUp(vector<FieldElement<T> >& v, FieldElement<T>& e) const
			throw(NotInSameFieldException, IllegalCoercionException);

	/****************** Comparison ******************/
		/* Two fields are the same only if they are the same
		 * object.
		 */
		bool operator==(const Field<T>& f) const throw () { return this==&f; }
		bool operator!=(const Field<T>& f) const throw () { return !(*this)==f; }
		/* Two fields are isomorphic if the isomorphism between them
		 * has actually been computed
		 */
		bool isIsomorphic(const Field<T>& f) const throw () { return stem==f.stem; }
		/* There is inclusion between two fields only if the inclusion
		 * has actually been computed.
		 */
		bool isSubFieldOf(const Field<T>& f) const throw ();
		bool isOverFieldOf(const Field<T>& f) const throw ()
		{ return f.isSubFieldOf(*this); }
	/****************** Printing ******************/
		ostream& print(ostream&) const;
	/****************** Destructor ******************/
		~Field() throw (ASException)
		{ throw ASException("Destroying fields is no good."); }


	/*****************************************************/
	/****************** Private section ******************/
	/*****************************************************/
	
	private:
	/****************** Copy prohibited ******************/
		void operator=(const Field<T>&);
		Field(const Field<T>&);
		
	/****************** Internal Constructors ******************/
		/* Construct a field with specified parameters */
		Field<T> (
			const Field<T>* sub,
			const Field<T>* over,
			const Context& ctxt,
			const FieldElement<T>* pri,
			const vector<FieldElement<T> >* pseudo,
			const FieldElement<T>* liftup,
			const MatGFp& mat,
			const long line,
			const bool pluso,
			const bool twopminuso,
			const Context& Ph,
			const Field<T>* st,
			const Field<T>* vsub,
			const BigInt& cha,
			const long deg,
			const long h,
			const FieldElement<T>* g,
			const FieldElement<T>* a
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
		p(cha), d(deg), height(h),
		gen(g), alpha(a)
		{}
		/* Private constructor for base fields */
		Field<T> (
			const Field<T>* sub,
			const Context& ctxt,
			const GFpE& pri,
			const MatGFp& mat,
			const long line,
			const BigInt& cha,
			const long deg,
			const GFpE& g
		) throw() :
		subfield(sub), overfield(),
		context(ctxt),
		primitive(new FieldElement<T>(this, pri)),
		pseudotraces(),
		liftuphelper(),
		artin(mat), artinLine(line),
		plusone(false), twopminusone(false),
		Phi(),
		stem(this), vsubfield(),
		p(cha), d(deg), height(0),
		gen(new FieldElement<T>(this, g)),
		alpha()
		{}
		/* Private constructor for GF(p) */
		Field<T> (
			const Context& ctxt,
			const GFp& pri,
			const BigInt& cha
		) throw() :
		subfield(), overfield(),
		context(ctxt),
		primitive(new FieldElement<T>(this, pri)),
		pseudotraces(),
		liftuphelper(),
		artin(), artinLine(),
		plusone(false), twopminusone(false),
		Phi(),
		stem(this), vsubfield(),
		p(cha), d(1), height(0),
		gen(new FieldElement<T>(this, pri)),
		alpha()
		{}

	};

	/****************** Printing ******************/
	template <class T> ostream& operator<<(ostream& o, const Field<T>& f) {
		return f.print(o);
	}
}



#include "../src/Field.c++"
#include "../src/FieldAlgorithms.c++"

#endif /*FIELD_H_*/
