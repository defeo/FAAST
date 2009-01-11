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
	
	class NotIrreducibleException : public ASException {};
	class IsIrreducibleException : public ASException {};
	class NotPrimeException : public ASException {};
	class NoOverFieldException : public ASException {};
	class NoSubFieldException : public ASException {};
	class NotASubFieldException : public ASException {};
	class NotInSameFieldException : public ASException {};
	class DivisionByZeroException : public ASException {};
	class IllegalCoercionException : public ASException {};
	class NotSupportedException : public ASException {
	public:
		NotSupportedException() {}
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
		BadParametersException() {}
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
