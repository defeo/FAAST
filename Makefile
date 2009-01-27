SHELL = /bin/bash
CC = g++
OPT = -g -pg -Wall -DAS_DEBUG=2 -DAS_TIMINGS

# Linked libraries
# For cvs compatibilty reasons, library files must be in /usr/local/lib,
# /usr/lib, /lib or some other directory specified by the $LIBRARY_PATH
# environment variable
GMPLIB = -lgmp
NTLLIB = -lntl_p 
MLIB = -lm
CPROFLIB =

# Project's files
ROOT = .
SRC = $(ROOT)/src
INC = $(ROOT)/include
BIN = $(ROOT)/bin/$(shell hostname)-$(shell uname -m)-$(CC)$(shell $(CC) -dumpversion)

# Parameters for gcc
# For cvs compatibility reasons, headers of linked libraries must be in
# /usr/include, /usr/local/include or some other directory specified by the
# $CPATH environment variable
IPATH=-I$(INC)
LIB=$(NTLLIB) $(GMPLIB) $(MLIB) $(ZLIB) $(CPROFLIB)

# Objects being prerequisites for every build
COMMONOBJS := 
COMMONOBJS := $(COMMONOBJS:%=$(BIN)/%)

# Script to create the bin directory
createbin:
	if [[ ! -d $(BIN) ]] ; then mkdir -p $(BIN); fi

########## binaries

TEST := test.o
TEST := $(TEST:%=$(BIN)/%)
$(BIN)/test: $(COMMONOBJS) $(TEST)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TEST) $(LIB) -o $(BIN)/test

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

########## .o files

$(BIN)/test.o: $(SRC)/test.c++ $(INC)/Types.hpp $(INC)/Field.hpp
	$(CC) -c $(OPT) $(IPATH) $(SRC)/test.c++ -o $(BIN)/test.o

$(BIN)/testStem.o: $(SRC)/testStem.c++ $(INC)/Types.hpp $(INC)/Field.hpp
	$(CC) -c $(OPT) $(IPATH) $(SRC)/testStem.c++ -o $(BIN)/testStem.o

$(BIN)/testTraceFrob.o: $(SRC)/testTraceFrob.c++ $(INC)/Types.hpp $(INC)/Field.hpp
	$(CC) -c $(OPT) $(IPATH) $(SRC)/testTraceFrob.c++ -o $(BIN)/testTraceFrob.o

$(BIN)/testLE.o: $(SRC)/testLE.c++ $(INC)/Types.hpp $(INC)/Field.hpp
	$(CC) -c $(OPT) $(IPATH) $(SRC)/testLE.c++ -o $(BIN)/testLE.o

$(BIN)/testCyclotomic.o: $(SRC)/testCyclotomic.c++ $(INC)/Types.hpp \
	$(INC)/Field.hpp $(INC)/utilities.hpp
	$(CC) -c $(OPT) $(IPATH) $(SRC)/testCyclotomic.c++ -o $(BIN)/testCyclotomic.o

$(BIN)/testTmul.o: $(SRC)/testTmul.c++ $(INC)/Tmul.hpp
	$(CC) -c $(OPT) $(IPATH) $(SRC)/testTmul.c++ -o $(BIN)/testTmul.o

########## other files

$(INC)/Field.hpp: $(INC)/Exceptions.hpp $(INC)/FieldElement.hpp \
	$(INC)/FieldPolynomial.hpp $(SRC)/Field.c++ $(SRC)/FieldAlgorithms.c++ \
	$(SRC)/FieldPrecomputations.c++ $(SRC)/Couveignes2000.c++
	touch $(INC)/Field.hpp

$(INC)/FieldElement.hpp: $(INC)/Exceptions.hpp \
	$(SRC)/FieldElement.c++ $(SRC)/FE-Liftup-Pushdown.c++ \
	$(SRC)/FE-Trace-Frob.c++
	touch $(INC)/FieldElement.hpp

$(INC)/FieldPolynomial.hpp: $(INC)/Exceptions.hpp $(SRC)/FieldPolynomial.c++
	touch $(INC)/FieldPolynomial.hpp

$(INC)/utilities.hpp: $(SRC)/utilities.c++
	touch $(INC)/utilities.hpp

$(SRC)/Field.c++: $(INC)/Types.hpp
	touch $(SRC)/Field.c++

$(SRC)/FieldElement.c++: $(INC)/Types.hpp
	touch $(SRC)/FieldElement.c++

$(SRC)/FieldPolynomial.c++: $(INC)/Types.hpp
	touch $(SRC)/FieldPolynomial.c++

$(SRC)/FieldAlgorithms.c++: $(INC)/utilities.hpp
	touch $(SRC)/FieldAlgorithms.c++

$(SRC)/FE-Liftup-Pushdown.c++: $(INC)/utilities.hpp $(INC)/Tmul.hpp
	touch $(SRC)/FE-Liftup-Pushdown.c++

$(SRC)/FE-Trace-Frob.c++: $(INC)/utilities.hpp
	touch $(SRC)/FE-Trace-Frob.c++

$(SRC)/FieldPrecomputations.c++: $(INC)/utilities.hpp
	touch $(SRC)/FieldPrecomputations.c++

######################################################################
.PHONY: all
all: createbin $(BIN)/test $(BIN)/testStem $(BIN)/testTraceFrob \
	$(BIN)/testLE $(BIN)/testCyclotomic $(BIN)/testTmul

.PHONY: now
now: createbin $(BIN)/test

.PHONY: clean
clean:
	rm -f $(BIN)/*.o
