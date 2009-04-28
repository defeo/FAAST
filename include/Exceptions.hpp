#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <exception>

namespace AS {
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

	class NotIrreducibleException : public ASException {
	public:
		virtual const char* what() const throw() {
			return "NotIrreducibleException";
		}
	};
	class IsIrreducibleException : public ASException {	public:
		virtual const char* what() const throw() {
			return "IsIrreducibleException";
		}
	};
	class NotPrimeException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NotPrimeException";
		}
	};
	class NoOverFieldException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NoOverFieldException";
		}
	};
	class NoSubFieldException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NoSubFieldException";
		}
	};
	class NotASubFieldException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NotASubFieldException";
		}
	};
	class NotInSameFieldException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NotInSameFieldException";
		}
	};
	class DivisionByZeroException : public ASException {	public:
		virtual const char* what() const throw() {
			return "DivisionByZeroException";
		}
	};
	class IllegalCoercionException : public ASException {	public:
		virtual const char* what() const throw() {
			return "IllegalCoercionException";
		}
	};
	class NoSuchPolynomialException : public ASException {	public:
		virtual const char* what() const throw() {
			return "NoSuchPolynomialException";
		}
	};
	class NotSupportedException : public ASException {
	public:
		NotSupportedException() : ASException("NotSupportedException") {}
		NotSupportedException(const char* m) : ASException(m) {}
	};
	class UndefinedFieldException : public ASException {
	public:
		virtual const char* what() const throw() {
			return "This is the 0 element of any field.";
		}
	};
	class BadParametersException : public ASException {
	public:
		BadParametersException() : ASException("BadParametersException") {}
		BadParametersException(const char* m) : ASException(m) {}
	};
	class CharacteristicTooLargeException : public ASException {
	public:
		virtual const char* what() const throw() {
			return "Do you really think it is wise to do computational\n Artin-Schreier theory in characteristic\n >= 2^(your machine word length) ?!";
		}
	};
}


#endif /*EXCEPTIONS_H_*/
