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
BIN = $(ROOT)/bin/$(shell uname -m)

# Parameters for gcc
# For cvs compatibility reasons, headers of linked libraries must be in
# /usr/include, /usr/local/include or some other directory specified by the
# $CPATH environment variable
IPATH=-I$(INC)
LIB=$(NTLLIB) $(GMPLIB) $(MLIB) $(ZLIB) $(CPROFLIB)

# Objects being prerequisites for every build
COMMONOBJS := 
COMMONOBJS := $(COMMONOBJS:%=$(BIN)/%)

########## binaries

#TESTCRV := testcrv.o montgomery.o io.o $(NTLFLDOBJS)
#TESTCRV := $(TESTCRV:%=$(BIN)/%)
#testcrv: $(COMMONOBJS) $(TESTCRV)
#	$(CC) $(OPT) $(IPATH) $(COMMONOBJS) $(TESTCRV) $(LIB) -o $(BIN)/testcrv

########## .o files

#$(BIN)/testcrv.o: $(SRC)/testcrv.c++ $(INC)/crvell.h $(INC)/options.h $(INC)/io.h
#	$(CC) -c $(OPT) $(IPATH) $(SRC)/testcrv.c++ -o $(BIN)/testcrv.o

########## other files

#$(INC)/options.h: $(INC)/crvell/twistinvariant.h $(INC)/crvell.h $(INC)/crvell/classic.h \
#	$(INC)/artintower.h $(INC)/artintower/ntl-based/tower.h


######################################################################
.PHONY: all
all: 

.PHONY: clean
clean:
	rm -f $(BIN)/*.o
