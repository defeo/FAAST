#ifndef FIELDELEMENT_H_
#define FIELDELEMENT_H_

#include "Exceptions.hpp"
#include <string>
#include <vector>

namespace AS {
	template <class T> class Field;
	
	template <class T> class FieldElement {
	
	friend class Field<T>;
	
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
		Field<T> parent() const throw(UndefinedFieldException)
		{ return parent_field; }
		
	/****************** Copy ******************/
//		FieldElement(const FieldElement<T>&) throw();
//		FieldElement<T>& operator=(const FieldElement<T> &) throw();
		
	/****************** Arithmetics ******************/
		/* Binary operations */
		FieldElement<T> operator+(const FieldElement<T>&)
			const throw(NotInSameFieldException);
		FieldElement<T> operator-(const FieldElement<T>&)
			const throw(NotInSameFieldException);
		FieldElement<T> operator*(const FieldElement<T>&)
			const throw(NotInSameFieldException);
		FieldElement<T> operator/(const FieldElement<T>&)
			const throw(NotInSameFieldException, DivisionByZeroException);
			
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
		FieldElement<T> operator-() const throw();
		FieldElement<T> inv() const throw(DivisionByZeroException);
		FieldElement<T> operator^(const ZZ&) const throw();
		FieldElement<T> operator^(const long) const throw();
		/* Frobenius and iterated frobenius */
		FieldElement<T> frobenius() const throw();
		FieldElement<T> frobenius(const long) const throw();
		/* Trace over the field F.
		 * 
		 * throws : NotASubFieldException if this does not belong to
		 *          an overfield of F.
		 */
		FieldElement<T> trace(const Field<T> F) const throw(NotASubFieldException);
		/* Absolute trace over GF(p) */
		FieldElement<T> trace() const throw();
		/* n-th pseudotrace */
		FieldElement<T> pseudotrace(const long n) const throw();
		
		/* Self-incrementing Unary operations */
	 	void negate() throw();
		void self_inv() throw(DivisionByZeroException);
		void operator^=(const ZZ&) throw();
		void operator^=(const long) throw();
		void self_frobenius() throw();
		void self_frobenius(const long) throw();
		void self_trace(const Field<T> F) throw(NotASubFieldException);
		void self_trace() throw();
		void self_pseudotrace(const long n) throw();
		
	/****************** Coercion of elements ******************/
		FieldElement<T> operator>>(const Field<T>&) const 
			throw(IllegalCoercionException);
		void operator>>=(const Field<T>&) throw(IllegalCoercionException);
		bool isCoercible(const Field<T>&) const throw();
		
	/****************** Comparison ******************/
		bool operator==(const FieldElement<T>&);
		bool operator==(const BigInt&);
		bool operator!=(const FieldElement<T>& e) { return !this==e; }
		bool operator!=(const BigInt& i) { return !this==i; }
		bool IsZero() const throw();
		bool IsOne() const throw();

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
		 * GF(p). The number of variables in vars must match one
		 * plus the Artin-Schreier height of the field the element
		 * belongs to.
		 * 
		 * throws : ASException if there's not enough variables
		 *          in var
		 */
		ostream& print(ostream&, vector<const string>& vars) const;
	/****************** Destructor ******************/
		~FieldElement() throw() {}


	/*****************************************************/
	/****************** Private section ******************/
	/*****************************************************/
	
	private:
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

	};

/****************** Level embedding ******************/
	/* Push the element e down along the stem and store
	 * the result in v.
	 * 
	 * throw : NoSubFieldException if e belongs to GF(p)
	 */
	template <class T>
	void pushDown(FieldElement<T>& e, vector<FieldElement<T> >& v)
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
	void liftUp(vector<FieldElement<T> >& v, FieldElement<T>& e)
	throw(NotInSameFieldException, NoOverFieldException);
}

#endif /*FIELDELEMENT_H_*/
