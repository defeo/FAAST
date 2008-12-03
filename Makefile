SHELL = /bin/bash
CC = g++
OPT=-g -pg -Wall -DAS_DEBUG=1

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
test: $(COMMONOBJS) $(TEST)
	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TEST) $(LIB) -o $(BIN)/test

########## .o files

$(BIN)/test.o: $(SRC)/test.c++ $(INC)/Types.h
	$(CC) -c $(OPT) $(IPATH) $(SRC)/test.c++ -o $(BIN)/test.o

########## other files

#$(INC)/options.h: $(INC)/crvell/twistinvariant.h $(INC)/crvell.h $(INC)/crvell/classic.h \
#	$(INC)/artintower.h $(INC)/artintower/ntl-based/tower.h
	

######################################################################
.PHONY: all
all: createbin test

.PHONY: clean
clean:
	rm -f $(BIN)/*.o
