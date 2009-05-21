SHELL = /bin/bash
CC = g++
OPT = -g -Wall -DAS_TIMINGS

# Linked libraries
# For cvs compatibilty reasons, library files must be in /usr/local/lib,
# /usr/lib, /lib or some other directory specified by the $LIBRARY_PATH
# environment variable
GMPLIB =
NTLLIB = -lntl 
MLIB =
CPROFLIB =

# Project's files
ROOT = .
SRC = $(ROOT)/src
TESTDIR = $(ROOT)/test
INC = $(ROOT)/include
BIN = $(ROOT)/bin/$(shell hostname)-$(shell uname -m)-$(CC)$(shell $(CC) -dumpversion)

# Parameters for gcc
# For cvs compatibility reasons, headers of linked libraries must be in
# /usr/include, /usr/local/include or some other directory specified by the
# $CPATH environment variable
IPATH=-I$(INC)
LIB=$(NTLLIB) $(GMPLIB) $(MLIB) $(ZLIB) $(CPROFLIB)

# Objects being prerequisites for every build
COMMONOBJS := Artin-Schreier.o
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

$(BIN)/Artin-Schreier.o: $(SRC)/explicit_instantiation.c++ $(INC)/Artin-Schreier.hpp \
	$(SRC)/Couveignes2000.c++ $(SRC)/FE-Liftup-Pushdown.c++ $(SRC)/FE-Trace-Frob.c++ \
	$(SRC)/Field.c++ $(SRC)/FieldAlgorithms.c++ $(SRC)/FieldElement.c++ $(SRC)/FieldPolynomial.c++ \
	$(SRC)/FieldPrecomputations.c++ $(SRC)/Minpols.c++ $(SRC)/utilities.c++ $(SRC)/NTLhacks.c++ \
	$(SRC)/Types.c++
	$(CC) -c $(OPT) $(IPATH) $(SRC)/explicit_instantiation.c++ -o $(BIN)/Artin-Schreier.o

########## object files for test routines

$(BIN)/test.o: $(TESTDIR)/test.c++ $(INC)/Artin-Schreier.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/test.c++ -o $(BIN)/test.o

$(BIN)/testIso.o: $(TESTDIR)/testIso.c++ $(INC)/Artin-Schreier.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testIso.c++ -o $(BIN)/testIso.o

$(BIN)/testStem.o: $(TESTDIR)/testStem.c++ $(INC)/Artin-Schreier.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testStem.c++ -o $(BIN)/testStem.o

$(BIN)/testTraceFrob.o: $(TESTDIR)/testTraceFrob.c++ $(INC)/Artin-Schreier.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testTraceFrob.c++ -o $(BIN)/testTraceFrob.o

$(BIN)/testLE.o: $(TESTDIR)/testLE.c++ $(INC)/Artin-Schreier.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testLE.c++ -o $(BIN)/testLE.o

$(BIN)/testCyclotomic.o: $(TESTDIR)/testCyclotomic.c++ $(INC)/Artin-Schreier.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testCyclotomic.c++ -o $(BIN)/testCyclotomic.o

$(BIN)/testTmul.o: $(TESTDIR)/testTmul.c++ $(INC)/AS/Tmul.hpp
	$(CC) -c $(OPT) $(IPATH) $(TESTDIR)/testTmul.c++ -o $(BIN)/testTmul.o

########## other files

$(INC)/Artin-Schreier.hpp: $(INC)/AS/Exceptions.hpp $(INC)/AS/Field.hpp $(INC)/AS/Types.hpp
	touch $(INC)/Artin-Schreier.hpp

$(INC)/AS/Types.hpp: $(INC)/AS/NTLhacks.hpp
	touch $(INC)/AS/Types.hpp

$(INC)/AS/Field.hpp: $(INC)/AS/Exceptions.hpp $(INC)/AS/FieldElement.hpp $(INC)/AS/FieldPolynomial.hpp
	touch $(INC)/AS/Field.hpp

$(INC)/AS/FieldElement.hpp: $(INC)/AS/Exceptions.hpp
	touch $(INC)/AS/FieldElement.hpp

$(INC)/AS/FieldPolynomial.hpp: $(INC)/AS/Exceptions.hpp
	touch $(INC)/AS/FieldPolynomial.hpp

$(INC)/AS/utilities.hpp: $(INC)/AS/Types.hpp $(INC)/AS/Exceptions.hpp
	touch $(INC)/AS/utilities.hpp

$(SRC)/FieldAlgorithms.c++: $(INC)/AS/utilities.hpp
	touch $(SRC)/FieldAlgorithms.c++

$(SRC)/FE-Liftup-Pushdown.c++: $(INC)/AS/utilities.hpp $(INC)/AS/Tmul.hpp
	touch $(SRC)/FE-Liftup-Pushdown.c++

$(SRC)/FE-Trace-Frob.c++: $(INC)/AS/utilities.hpp
	touch $(SRC)/FE-Trace-Frob.c++

$(SRC)/FieldPrecomputations.c++: $(INC)/AS/utilities.hpp
	touch $(SRC)/FieldPrecomputations.c++

######################################################################
########## Obsolete

# Artin-Schreier.hpp.gch: $(INC)/Artin-Schreier.hpp
#	$(CC) -x c++-header -c $(OPT) $(IPATH) $(INC)/Artin-Schreier.hpp -o Artin-Schreier.hpp.gch

######################################################################
########## Other targets

.PHONY: library
library: $(COMMONOBJS)
	ar rcs libArtin-Schreier.a $(BIN)/Artin-Schreier.o

.PHONY: test
test: createbin $(BIN)/test $(BIN)/testIso $(BIN)/testStem \
	$(BIN)/testTraceFrob $(BIN)/testLE $(BIN)/testCyclotomic \
	$(BIN)/testTmul

.PHONY: doc
doc:
	doxygen doxy.conf

.PHONY: doc-dev
doc-dev:
	(cat doxy.conf; echo "ENABLED_SECTIONS=DEV"; echo "OUTPUT_DIRECTORY=doc-dev";) | doxygen -

.PHONY: clean
clean:
	rm -f $(BIN)/*.o
