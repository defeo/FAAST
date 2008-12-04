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
	class UndefinedFieldException : public ASException {
		virtual const char* what() const throw() {
			return "This is the 0 element of any field.";
		}
	};
}


#endif /*EXCEPTIONS_H_*/
