AM_CPPFLAGS = -DFAAST_TIMINGS -I$(srcdir)/../include

lib_LTLIBRARIES = libfaast.la
libfaast_la_SOURCES = explicit_instantiation.c++ Types.hpp \
	Couveignes2000.hpp FE-Liftup-Pushdown.hpp FE-Trace-Frob.hpp \
	Field.hpp FieldAlgorithms.hpp FieldElement.hpp FieldPolynomial.hpp \
	FieldPrecomputations.hpp GCD.hpp Minpols.hpp utilities.hpp NTLhacks.hpp
libfaast_la_LDFLAGS = -versioninfo 1:0:0
