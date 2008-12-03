#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <exception>

namespace AS {
	class ASException : public std::exception {};
	class NotIrreducibleException : public ASException {};
	class IsIrreducibleException : public ASException {};
	class NotPrimeException : public ASException {};
	class NoOverFieldException : public ASException {};
	class NoSubFieldException : public ASException {};
	class NotInSameFieldException : public ASException {};
}


#endif /*EXCEPTIONS_H_*/
