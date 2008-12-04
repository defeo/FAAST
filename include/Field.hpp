#ifndef FIELD_H_
#define FIELD_H_

#include "Exceptions.hpp"
#include "FieldElement.hpp"

namespace AS {
	template <class T> class Field {
	private:
		typedef typename T::GFp         GFp;
		typedef typename T::MatGFp      MatGFp;
		typedef typename T::Poly        Poly;
		typedef typename T::ModPoly     ModPoly;
		typedef typename T::PolyModPoly PolyModPoly;
		typedef typename T::BigInt      BigInt;
		typedef typename T::Context     Context;
		
	/****************** Members for the stem ******************/
		/* The immediate subfield and overfield, if they're defined */
		Field<T>* subfield;
		Field<T>* overfield;
		/* Infrastructure-dependent data to perform computations in
		 * the field. (e.g. GF2EContext in NTL)
		 */
		Context context;
		/* The generator over GF(p) */
		const FieldElement<T> primitive;
		/* Precomputed pseudotraces */
		vector<const FieldElement<T> > pseudotraces;
		/* Lift-up precomputation */
		const FieldElement<T> liftuphelper;
		/* The inverse matrix of the linear application X^p-X */
		const MatGFp artin;
		/* Flags related to the construction of the extension */
		const bool plusone;
		const bool cube;

	/****************** Members for non-stem fields ******************/
		/* The stem-field to which this one is isomorphic */
		Field<T>* stem;
		
	/****************** Common members ******************/
		/* The virtual subfield, if it is different from the real.
		 * It is given mainly for printing purposes
		 */
		Field<T>* vsubfield;
		/* the characteristic */
		const long p;
		/* the degree */
		const long d;
		/* The Artin-Schreier height */
		const long height;
		/* The generator over the subfield */
		const FieldElement<T> gen;
		/* The element of the subfield such that this field is
		 * the Artin-Schreier extension defined by the polynomial
		 * 			X^p - X - alpha
		 */
		const FieldElement<T> alpha;
		

	public:
	/****************** Constructors ******************/
		/* Default constructor, builds a field from some default value
		 * (i.e. the context for NTL).
		 * 
		 * throws : NotIrreducibleException if the default makes no sense
		 */
		Field() throw (NotIrreducibleException);
		/* Build a field from an irreducible polynomial P.
		 * 
		 * throws : NotIrreducibleException if P is not irreducible
		 */
		Field(const Poly& P) throw (NotIrreducibleException);
		/* Build the field GF(p^d) using some default
		 * polynomial.
		 * Notice that this operation implicitely creates
		 * the field GF(p) too.
		 * 
		 * throws : NotPrimeException if p is not prime
		 */
		Field(long p, long d = 1) throw (NotPrimeException);
	
	/****************** Field Extensions ******************/
		/* Build a default extension of degree p over this field. 
		 */
		Field<T> ArtinSchreierExtension() const throw ();
		/* Build the splitting field of the polynomial
		 * 			X^p - X - alpha
		 * over this field. This may or may not be an extension of
		 * degree p depending if the polynomial is irreducible.
		 */
		Field<T> ArtinSchreierExtension(const FieldElement<T>& alpha) const throw ();
		
	/****************** Properties ******************/
		long characteristic() const throw ();
		long degree() const throw ();
		BigInt cardinality() const throw ();
		/* The number of Artin-Schreier extensions from the base
		 * field.
		 */
		long ArtinSchreierHeight() const throw ();
		
	/****************** Field Elements ******************/
		/* Constructs the element i times 1 */
		FieldElement<T> scalar(const long i) const throw ();
		FieldElement<T> zero() const throw () { return scalar(0); }
		FieldElement<T> one() const throw () { return scalar(1); }
		/* The generator over the immediately preceding subfield */
		FieldElement<T> generator() const throw ();
		/* The generator over GF(p) */
		FieldElement<T> primitiveElement() const throw ();
		/* A random element of the field */
		FieldElement<T> random() const throw ();
	
	/****************** Field lattice navigation ******************/
		/* The field GF(p) */
		Field<T> baseField() const throw();

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
		bool operator==(const Field<T>& f) const throw ();
		bool operator!=(const Field<T>& f) const throw () { return !this==f; }
		/* Two fields are isomorphic if the isomorphism between them
		 * has actually been computed
		 */
		bool isIsomporphic(const Field<T>& f) const throw ();
		/* There is inclusion between two fields only if the inclusion
		 * has actually been computed.
		 */
		bool isSubFieldOf(const Field<T>& f) const throw ();
		bool isOverFieldOf(const Field<T>& f) const throw ();
	/****************** Printing ******************/
		ostream& print(ostream&) const;
	/****************** Destructor ******************/
		~Field() throw ();


	/*****************************************************/
	/****************** Private section ******************/
	/*****************************************************/
	
	private:
	/****************** Copy prohibited ******************/
		void operator=(const Field<T>&);
		Field(const Field<T>&);
		
	/****************** Internal Constructors ******************/
		/* For constructing field extensions */
		//Field<T> () throw();
	};
}

#endif /*FIELD_H_*/
