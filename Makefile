SHELL = /bin/bash
CC = g++
OPT=-g -pg -Wall -DAS_DEBUG=3 -DAS_TIMINGS

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
	if [[ ! -d $(BIN) ]] ; then mkdir $(BIN); fi

########## binaries

TEST := test.o
TEST := $(TEST:%=$(BIN)/%)
$(BIN)/test: $(COMMONOBJS) $(TEST)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TEST) $(LIB) -o $(BIN)/test

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

$(BIN)/testCyclotomic.o: $(SRC)/testCyclotomic.c++ $(INC)/Types.hpp \
	$(INC)/Field.hpp $(INC)/utilities.hpp
	$(CC) -c $(OPT) $(IPATH) $(SRC)/testCyclotomic.c++ -o $(BIN)/testCyclotomic.o

$(BIN)/testTmul.o: $(SRC)/testTmul.c++ $(INC)/Tmul.hpp
	$(CC) -c $(OPT) $(IPATH) $(SRC)/testTmul.c++ -o $(BIN)/testTmul.o

########## other files

$(INC)/Field.hpp: $(INC)/Exceptions.hpp $(INC)/FieldElement.hpp \
	$(INC)/FieldPolynomial.hpp $(SRC)/Field.c++ $(SRC)/FieldAlgorithms.c++
	touch $(INC)/Field.hpp

$(INC)/FieldElement.hpp: $(INC)/Exceptions.hpp \
	$(SRC)/FieldElement.c++ $(SRC)/FieldElementAlgorithms.c++ 
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

$(SRC)/FieldElementAlgorithms.c++: $(INC)/utilities.hpp $(INC)/Tmul.hpp
	touch $(SRC)/FieldElementAlgorithms.c++

######################################################################
.PHONY: all
all: createbin $(BIN)/test $(BIN)/testCyclotomic $(BIN)/testTmul

.PHONY: clean
clean:
	rm -f $(BIN)/*.o
