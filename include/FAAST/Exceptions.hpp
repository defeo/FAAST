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
#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <exception>

namespace FAAST {
	/**
	 * \defgroup Exceptions Exceptions
	 * These are the exceptions thrown by the functions of FAAST.
	 * @{
	 */

	/**
	 * \brief A generic exception.
	 */
	class FAASTException : public std::exception {
	private:
		const char* message;
	public:
		FAASTException() {}
		FAASTException(const char* m) : message(m) {}
		virtual const char* what() const throw() {
			return message;
		}
	};

	/**
	 * \brief The polynomial is not irreducible.
	 */
	class NotIrreducibleException : public FAASTException {
	public:
		NotIrreducibleException() : FAASTException("NotIrreducibleException") {}
	};
	/**
	 * \brief The polynomial is irreducible.
	 */
	class IsIrreducibleException : public FAASTException {	public:
		IsIrreducibleException() : FAASTException("IsIrreducibleException") {}
	};
	/**
	 * \brief The number is not a prime.
	 */
	class NotPrimeException : public FAASTException {	public:
		NotPrimeException() : FAASTException("NotPrimeException") {}
	};
	/**
	 * \brief No overfield is known for the given field.
	 */
	class NoOverFieldException : public FAASTException {	public:
		NoOverFieldException() : FAASTException("NoOverFieldException") {}
	};
	/**
	 * \brief No subfield is known for the given field.
	 */
	class NoSubFieldException : public FAASTException {	public:
		NoSubFieldException() : FAASTException("NoSubFieldException") {}
	};
	/**
	 * \brief The field is not a subfield of the given field.
	 */
	class NotASubFieldException : public FAASTException {	public:
		NotASubFieldException() : FAASTException("NotASubFieldException") {}
	};
	/**
	 * \brief The two elements do not belong to the same field.
	 */
	class NotInSameFieldException : public FAASTException {	public:
		NotInSameFieldException() : FAASTException("NotInSameFieldException") {}
	};
	/**
	 * \brief Division by zero.
	 */
	class DivisionByZeroException : public FAASTException {
	public:
		DivisionByZeroException() : FAASTException("DivisionByZeroException") {}
	};
	/**
	 * \brief The element cannot be coerced to the given field.
	 */
	class IllegalCoercionException : public FAASTException {	public:
		IllegalCoercionException() : FAASTException("IllegalCoercionException") {}
	};
	/**
	 * \brief There is no polynomial with the given property.
	 */
	class NoSuchPolynomialException : public FAASTException {	public:
		NoSuchPolynomialException() : FAASTException("NoSuchPolynomialException") {}
	};
	/**
	 * \brief The operation is not supported (yet?)
	 */
	class NotSupportedException : public FAASTException {
	public:
		NotSupportedException() : FAASTException("NotSupportedException") {}
		NotSupportedException(const char* m) : FAASTException(m) {}
	};
	/**
	 * \brief The element has no parent field specified.
	 *
	 * You are trying to make a field-dependent operation on an element
	 * for which no parent field has been specified yet.
	 *
	 * A typical mistake leading to this exception is trying to assign
	 * or sum a scalar to a FieldElement before having made the FieldElement
	 * member of any Field.
	 * \code
	 * FieldElement<T> x;
	 * x += 1;
	 * \endcode
	 * In this example, after the creation \c x is 0,
	 * but it does not belong
	 * to any field. You can add \c x to any element of
	 * any field, or even to
	 * another generic 0 element, but you cannot add \c x
	 * to a scalar since the
	 * library wouldn't know modulo which integer to reduce the scalar.
	 *
	 * A correct example, assuming \c K has been defined
	 * to be a field, would be
	 * \code
	 * FieldElement<T> x = K.zero();
	 * x += 1;
	 * \endcode
	 * or, equivalently,
	 * \code
	 * FieldElement<T> x;
	 * x += K.scalar(1);
	 * \endcode
	 */
	class UndefinedFieldException : public FAASTException {
	public:
		UndefinedFieldException() : FAASTException("This is the 0 element of any field.") {}
	};
	/**
	 * \brief The function does not accept such parameters.
	 */
	class BadParametersException : public FAASTException {
	public:
		BadParametersException() : FAASTException("BadParametersException") {}
		BadParametersException(const char* m) : FAASTException(m) {}
	};
	/**
	 * \brief The characteristic is larger than what AS can handle.
	 */
	class CharacteristicTooLargeException : public FAASTException {
	public:
		CharacteristicTooLargeException() : FAASTException("Do you really think it is wise to do computational\n Artin-Schreier theory in characteristic\n >= 2^(your machine word length) ?!") {}
	};

	/**
	 * @}
	 */
}


#endif /*EXCEPTIONS_H_*/
