#ifndef FIELD_H_
#define FIELD_H_

template <class T> class Field {
public:
	typedef typename T::GFp         GFp;
	typedef typename T::Poly        Poly;
	typedef typename T::ModPoly     ModPoly;
	typedef typename T::PolyModPoly PolyModPoly;
	typedef typename T::BigInt      BigInt;
	

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
	 * 
	 * throws : NotPrimeException if p is not prime
	 */
	Field(long p, long d = 1) throw (NotPrimeException);

/****************** Field Extensions ******************/
	/* Build a default extension of degree p over this field. 
	 */
	Field<T> ArtinSchreierExtension() const;
	/* Build the splitting field of the polynomial
	 * 			X^p - X - alpha
	 * over this field. This may or may not be an extension of
	 * degree p depending if the polynomial is irreducible.
	 */
	Field<T> ArtinSchreierExtension(const FieldElement<T>& alpha) const;
	
/****************** Properties ******************/
	long characteristic() const;
	long degree() const;
	BigInt cardinality() const;
	/* The number of Artin-Schreier extensions from the base
	 * field.
	 */
	long ArtinSchreierHeight() const;
	
/****************** Field Elements ******************/
	/* Constructs the element i times 1 */
	FieldElement<T> scalar(const long i) const;
	FieldElement<T> zero() const { return scalar(0); }
	FieldElement<T> one() const { return scalar(1); }
	/* The generator over the immediately preceding subfield */
	FieldElement<T> generator() const;
	/* The generator over GF(p) */
	FieldElement<T> primitiveElement() const;
	/* A random element of the field */
	FieldElement<T> random() const;

/****************** Comparison ******************/
	/* Two fields are the same if they are defined by the same
	 * polynomial over GF(p).
	 */
	bool operator==(const Field<T> f) const;
	bool operator!=(const Field<T> f) const { return !this==f; }
	/* Two fields are isomorphic if the isomorphism between them
	 * has actually been computed
	 */
	bool isIsomporphic(const Field<T> f) const;
	/* There is inclusion between two fields only if the inclusion
	 * has actually been computed.
	 */
	bool isSubFieldOf(const Field<T> f) const;
	bool isOverFieldOf(const Field<T> f) const;
/****************** Printing ******************/
	ostream& print(ostream&) const;
/****************** Destructor ******************/
	~Field();
/****************** Copy prohibited ******************/
private:
	void operator=(const Field<T>&);
};


#endif /*FIELD_H_*/
