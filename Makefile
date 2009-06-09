##############################################################################
##                            Configuration section                         ##
## Edit this variables if you installed NTL, GMP or GF2X in some exotic     ##
## directory, if you compiled NTL as a static library or if you want to     ##
## install the library in a different location.                             ##
##############################################################################

# The paths were NTL, GMP and gf2x libraries are to be found.
NTLLIBPATH = # -L$HOME/lib

# The paths were NTL, GMP and gf2x headers are to be found.
NTLINCPATH = # -I$HOME/include  

# The libraries to link to. If you link to ntl as a static library, you must
# uncomment the end of the line. Leave this line as is if you link to ntl as
# a shared library.
NTLLIB = -lntl # -lgmp -lgf2x -lm

PREFIX = /usr
##############################################################################
##                          End of configuration section                    ##
##############################################################################


SHELL = /bin/bash
CC = g++
OPT = -g -Wall -DFAAST_TIMINGS

# FAAST files
ROOT = .
SRC = $(ROOT)/src
TESTDIR = $(ROOT)/test
INC = $(ROOT)/include
BIN = $(ROOT)/bin/$(shell hostname)-$(shell uname -m)-$(CC)$(shell $(CC) -dumpversion)

# Parameters for gcc
IPATH=-I$(INC) $(NTLINCPATH)
LIB=$(NTLLIBPATH) $(NTLLIB)

# Objects being prerequisites for every build
COMMONOBJS := faast.o
COMMONOBJS := $(COMMONOBJS:%=$(BIN)/%)

########## targets
# the default target
.PHONY: all
all: createbin library

# Script to create the bin directory
createbin:
	if [[ ! -d $(BIN) ]] ; then mkdir -p $(BIN); fi

########## binaries

TEST := test.o
TEST := $(TEST:%=$(BIN)/%)
$(BIN)/test: $(COMMONOBJS) $(TEST)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TEST) $(LIB) -o $(BIN)/test

TESTISO := testIso.o
TESTISO := $(TESTISO:%=$(BIN)/%)
$(BIN)/testIso: $(COMMONOBJS) $(TESTISO)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TESTISO) $(LIB) -o $(BIN)/testIso

TESTSTEM := testStem.o
TESTSTEM := $(TESTSTEM:%=$(BIN)/%)
$(BIN)/testStem: $(COMMONOBJS) $(TESTSTEM)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TESTSTEM) $(LIB) -o $(BIN)/testStem

TESTTF := testTraceFrob.o
TESTTF := $(TESTTF:%=$(BIN)/%)
$(BIN)/testTraceFrob: $(COMMONOBJS) $(TESTTF)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TESTTF) $(LIB) -o $(BIN)/testTraceFrob

TESTLE := testLE.o
TESTLE := $(TESTLE:%=$(BIN)/%)
$(BIN)/testLE: $(COMMONOBJS) $(TESTLE)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TESTLE) $(LIB) -o $(BIN)/testLE

TESTCYCLO := testCyclotomic.o
TESTCYCLO := $(TESTCYCLO:%=$(BIN)/%)
$(BIN)/testCyclotomic: $(COMMONOBJS) $(TESTCYCLO)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TESTCYCLO) $(LIB) -o $(BIN)/testCyclotomic

TESTTMUL := testTmul.o
TESTTMUL := $(TESTTMUL:%=$(BIN)/%)
$(BIN)/testTmul: $(COMMONOBJS) $(TESTTMUL)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TESTTMUL) $(LIB) -o $(BIN)/testTmul

########## object files

$(BIN)/faast.o: $(SRC)/explicit_instantiation.c++ $(INC)/faast.hpp \
	$(SRC)/Couveignes2000.c++ $(SRC)/FE-Liftup-Pushdown.c++ $(SRC)/FE-Trace-Frob.c++ \
	$(SRC)/Field.c++ $(SRC)/FieldAlgorithms.c++ $(SRC)/FieldElement.c++ $(SRC)/FieldPolynomial.c++ \
	$(SRC)/FieldPrecomputations.c++ $(SRC)/Minpols.c++ $(SRC)/utilities.c++ $(SRC)/NTLhacks.c++ \
	$(SRC)/Types.c++
	$(CC) -c $(OPT) $(IPATH) $(SRC)/explicit_instantiation.c++ -o $(BIN)/faast.o

########## object files for test routines

$(BIN)/test.o: $(TESTDIR)/test.c++ $(INC)/faast.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/test.c++ -o $(BIN)/test.o

$(BIN)/testIso.o: $(TESTDIR)/testIso.c++ $(INC)/faast.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testIso.c++ -o $(BIN)/testIso.o

$(BIN)/testStem.o: $(TESTDIR)/testStem.c++ $(INC)/faast.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testStem.c++ -o $(BIN)/testStem.o

$(BIN)/testTraceFrob.o: $(TESTDIR)/testTraceFrob.c++ $(INC)/faast.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testTraceFrob.c++ -o $(BIN)/testTraceFrob.o

$(BIN)/testLE.o: $(TESTDIR)/testLE.c++ $(INC)/faast.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testLE.c++ -o $(BIN)/testLE.o

$(BIN)/testCyclotomic.o: $(TESTDIR)/testCyclotomic.c++ $(INC)/faast.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testCyclotomic.c++ -o $(BIN)/testCyclotomic.o

$(BIN)/testTmul.o: $(TESTDIR)/testTmul.c++ $(INC)/FAAST/Tmul.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testTmul.c++ -o $(BIN)/testTmul.o

########## other files

$(INC)/faast.hpp: $(INC)/FAAST/Exceptions.hpp $(INC)/FAAST/Field.hpp $(INC)/FAAST/Types.hpp
	touch $(INC)/faast.hpp

$(INC)/FAAST/Types.hpp: $(INC)/FAAST/NTLhacks.hpp
	touch $(INC)/FAAST/Types.hpp

$(INC)/FAAST/Field.hpp: $(INC)/FAAST/Exceptions.hpp $(INC)/FAAST/FieldElement.hpp $(INC)/FAAST/FieldPolynomial.hpp
	touch $(INC)/FAAST/Field.hpp

$(INC)/FAAST/FieldElement.hpp: $(INC)/FAAST/Exceptions.hpp
	touch $(INC)/FAAST/FieldElement.hpp

$(INC)/FAAST/FieldPolynomial.hpp: $(INC)/FAAST/Exceptions.hpp
	touch $(INC)/FAAST/FieldPolynomial.hpp

$(INC)/FAAST/utilities.hpp: $(INC)/FAAST/Types.hpp $(INC)/FAAST/Exceptions.hpp
	touch $(INC)/FAAST/utilities.hpp

$(SRC)/FieldAlgorithms.c++: $(INC)/FAAST/utilities.hpp
	touch $(SRC)/FieldAlgorithms.c++

$(SRC)/FE-Liftup-Pushdown.c++: $(INC)/FAAST/utilities.hpp $(INC)/FAAST/Tmul.hpp
	touch $(SRC)/FE-Liftup-Pushdown.c++

$(SRC)/FE-Trace-Frob.c++: $(INC)/FAAST/utilities.hpp
	touch $(SRC)/FE-Trace-Frob.c++

$(SRC)/FieldPrecomputations.c++: $(INC)/FAAST/utilities.hpp
	touch $(SRC)/FieldPrecomputations.c++

######################################################################
########## Other targets

.PHONY: library
library: $(COMMONOBJS)
	ar rcs libfaast.a $(BIN)/faast.o

.PHONY: install
install: library
	cp libfaast.a $(PREFIX)/lib/
	cp -r include/* $(PREFIX)/include/

.PHONY: examples
examples: createbin $(BIN)/test $(BIN)/testIso $(BIN)/testStem \
	$(BIN)/testTraceFrob $(BIN)/testLE $(BIN)/testCyclotomic \
	$(BIN)/testTmul

.PHONY: doc
doc:
	doxygen doxy.conf

.PHONY: doc-dev
doc-dev:
	(cat doxy.conf; \
	echo "ENABLED_SECTIONS=DEV"; \
	echo "OUTPUT_DIRECTORY=doc-dev"; \
	echo "PREDEFINED = FAAST_TIMINGS"; \
	echo "SHOW_USED_FILES = YES"; \
	echo "SHOW_FILES = YES") | doxygen -

.PHONY: package
package:
	tar czf faast-$Name:  $.tgz .

.PHONY: clean
clean:
	rm -f $(BIN)/*.o
