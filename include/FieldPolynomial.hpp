#ifndef FIELDPOLYNOMIAL_H_
#define FIELDPOLYNOMIAL_H_

#include "Exceptions.hpp"
#include <string>
#include <vector>

namespace AS {
	template <class T> class Field;
	
	template <class T> class FieldPolynomial {
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
		
		
	/****************** Members ******************/
		/* The representation of this element */
		GFpX repBase;
		GFpEX repExt;
		bool base;
		/* The field this element belongs to */
		const Field<T>* parent_field;


	public:
	/****************** Constructors ******************/
		/* Constructor by default.
		 * Constructs the 0 polynomial (over any field).
		 */
		FieldPolynomial() throw() : parent_field(NULL) {}
	/****************** Properties ******************/
		/* The field this polynomial is defined over */
		Field<T> parent() const throw(UndefinedFieldException)
		{ return parent_field; }
		/* Degree. Returns -1 if the polynomial is 0. */
		long degree() const throw();
		
	/****************** Copy ******************/
		FieldPolynomial(const FieldPolynomial<T>&) throw();
		FieldPolynomial<T>& operator=(const FieldPolynomial<T>&) throw();
		/* Starting only from the constant coefficient */
		FieldPolynomial(const FieldElement<T>& e) throw();
		FieldPolynomial<T>& operator=(const FieldElement<T>& e) throw();
		
	/****************** Coefficients ******************/
		void getCoeff(const long, FieldElement<T>&) const throw();
		void setCoeff(const long, const FieldElement<T>&) throw();
		void setCoeff(const long i, const BigInt& c)
		throw(UndefinedFieldException)
		{ setCoeff(i, parent_field->scalar(c)); }
		void setCoeff(const long i)
		throw(UndefinedFieldException)
		{ setCoeff(i, parent_field->one()); }
	
	/****************** Arithmetics ******************/
		/* Binary operations */
		FieldPolynomial<T> operator+(const FieldPolynomial<T>&)
			const throw(NotInSameFieldException);
		FieldPolynomial<T> operator-(const FieldPolynomial<T>&)
			const throw(NotInSameFieldException);
		FieldPolynomial<T> operator*(const FieldPolynomial<T>&)
			const throw(NotInSameFieldException);
		FieldPolynomial<T> operator/(const FieldPolynomial<T>&)
			const throw(NotInSameFieldException, DivisionByZeroException);
		FieldPolynomial<T> operator%(const FieldPolynomial<T>&)
			const throw(NotInSameFieldException, DivisionByZeroException);
			
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

		/* Memory-efficient (NTL-like) binary operations.
		 * Aplly on the arguments and store in this.
		 */
		void sum(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException) 
		{ operator=(a); operator+=(b); }
		void difference(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator-=(b); }
		void product(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException)
		{ operator=(a); operator*=(b); }
		void division(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException, DivisionByZeroException)
		{ operator=(a); operator/=(b); }
		void mod(const FieldPolynomial<T>& a, const FieldPolynomial<T>& b)
			throw(NotInSameFieldException, DivisionByZeroException)
		{ operator=(a); operator%=(b); }
			
		/* Unary operations */
		bool divides(const FieldPolynomial<T>&) const throw();
		FieldPolynomial<T> operator-() const throw();
		FieldPolynomial<T> operator^(const long) const throw();
		FieldPolynomial<T> derivative() const throw();
		FieldPolynomial<T> monic() const throw();
		/* Frobenius and iterated frobenius over the coefficients */
		FieldPolynomial<T> frobenius() const throw();
		FieldPolynomial<T> frobenius(const long) const throw();
		
		/* Self-incrementing Unary operations */
		void negate() throw();
		void operator^=(const long) throw();
		void self_derivative() throw();
		void normalize() throw();
		void self_frobenius() throw();
		void self_frobenius(const long) throw();
		
	/****************** Coercion of elements ******************/
		FieldPolynomial<T> operator>>(const Field<T>&) const 
			throw(IllegalCoercionException);
		void operator>>=(const Field<T>&) throw(IllegalCoercionException);
		bool isCoercible(const Field<T>&) const throw();
		
	/****************** Comparison ******************/
		bool operator==(const FieldPolynomial<T>&) const throw(NotInSameFieldException);
		bool operator!=(const FieldPolynomial<T>& e)
		const throw(NotInSameFieldException)
		{ return !this==e; }
		bool isZero() const throw();
		bool isOne() const throw();

	/****************** Infrastructure ******************/
		/* Interface with infrastructure. Use this only if you are
		 * sure of what you do !
		 */
		void toInfrastructure(GFpX& ) const throw(IllegalCoercionException);
		void toInfrastructure(GFpEX& ) const throw(IllegalCoercionException);

	/****************** Printing ******************/
		ostream& print(ostream&) const;
		/* Print the element as a polynomial over the base field in the
		 * variable varPoly. varField is used to print the element of the field.
		 */
		ostream& print(ostream&, const string& varPoly, const string& varField) const;
		/* Print the element as a polynomial over the base field in the
		 * variable varPoly.
		 * varsField is used to print the element of the field on the
		 * multivariate basis. The number of variables in vars must match
		 * one plus the Artin-Schreier height of the base field.
		 * 
		 * throws : ASException if there's not enough variables
		 *          in var
		 */
		ostream& print(ostream&, const string& varPoly, vector<const string>& varsField) const;
	/****************** Destructor ******************/
		~FieldPolynomial() throw() {}


	/*****************************************************/
	/****************** Private section ******************/
	/*****************************************************/
	
	private:
	/****************** Internal Constructors ******************/
		/* Construct an element with given representation and parent.
		 * Reserved for use by Field<T>
		 */
		FieldPolynomial(const Field<T>* p, const GFpX& PBase, const GFpEX& PExt, const bool b) throw() :
			parent_field(p), repBase(PBase), repExt(PExt), base(b) {}
		FieldPolynomial(const Field<T>* p, const GFpEX& P) throw() :
			parent_field(p), repExt(P), base(false) {}
		FieldPolynomial(const Field<T>* p, const GFpX& P) throw() :
			parent_field(p), repBase(P), base(true) {}

	};
}

#endif /*FIELDPOLYNOMIAL_H_*/
