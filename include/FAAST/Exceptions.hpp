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
	 * These are the exceptions thrown by the functions of AS.
	 * @{
	 */

	/**
	 * \brief A generic exception.
	 */
	class ASException : public std::exception {
	private:
		const char* message;
	public:
		ASException() {}
		ASException(const char* m) : message(m) {}
		virtual const char* what() const throw() {
			return message;
		}
	};

	/**
	 * \brief The polynomial is not irreducible.
	 */
	class NotIrreducibleException : public ASException {
	public:
		virtual const char* what() const throw() {
			return "NotIrreducibleException";
		}
	};
	/**
	 * \brief The polynomial is irreducible.
	 */
	class IsIrreducibleException : public ASException {	public:
		virtual const char* what() const throw() {
			return "IsIrreducibleException";
		}
	};
	/**
	 * \brief The number is not a prime.
	 */
	class NotPrimeException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NotPrimeException";
		}
	};
	/**
	 * \brief No overfield is known for the given field.
	 */
	class NoOverFieldException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NoOverFieldException";
		}
	};
	/**
	 * \brief No subfield is known for the given field.
	 */
	class NoSubFieldException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NoSubFieldException";
		}
	};
	/**
	 * \brief The field is not a subfield of the given field.
	 */
	class NotASubFieldException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NotASubFieldException";
		}
	};
	/**
	 * \brief The two elements do not belong to the same field.
	 */
	class NotInSameFieldException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NotInSameFieldException";
		}
	};
	/**
	 * \brief Division by zero.
	 */
	class DivisionByZeroException : public ASException {	public:
		virtual const char* what() const throw() {
			return "DivisionByZeroException";
		}
	};
	/**
	 * \brief The element cannot be coerced to the given field.
	 */
	class IllegalCoercionException : public ASException {	public:
		virtual const char* what() const throw() {
			return "IllegalCoercionException";
		}
	};
	/**
	 * \brief There is no polynomial with the given property.
	 */
	class NoSuchPolynomialException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NoSuchPolynomialException";
		}
	};
	/**
	 * \brief The operation is not supported (yet?)
	 */
	class NotSupportedException : public ASException {
	public:
		NotSupportedException() : ASException("NotSupportedException") {}
		NotSupportedException(const char* m) : ASException(m) {}
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
	class UndefinedFieldException : public ASException {
	public:
		virtual const char* what() const throw() {
			return "This is the 0 element of any field.";
		}
	};
	/**
	 * \brief The function does not accept such parameters.
	 */
	class BadParametersException : public ASException {
	public:
		BadParametersException() : ASException("BadParametersException") {}
		BadParametersException(const char* m) : ASException(m) {}
	};
	/**
	 * \brief The characteristic is larger than what AS can handle.
	 */
	class CharacteristicTooLargeException : public ASException {
	public:
		virtual const char* what() const throw() {
			return "Do you really think it is wise to do computational\n Artin-Schreier theory in characteristic\n >= 2^(your machine word length) ?!";
		}
	};

	/**
	 * @}
	 */
}


#endif /*EXCEPTIONS_H_*/
