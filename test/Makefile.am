AM_CPPFLAGS = -DFAAST_TIMINGS -I$(srcdir)/../include

check_PROGRAMS = test testCyclotomic testGCD testIso testLE testStem testTraceFrob
test_SOURCES = test.c++
test_LDADD = ../src/libfaast.la
testCyclotomic_SOURCES = testCyclotomic.c++
testCyclotomic_LDADD = ../src/libfaast.la
testGCD_SOURCES = testGCD.c++
testGCD_LDADD = ../src/libfaast.la
testIso_SOURCES = testIso.c++
testIso_LDADD = ../src/libfaast.la
testLE_SOURCES = testLE.c++
testLE_LDADD = ../src/libfaast.la
testNTLmul_SOURCES = testNTLmul.c++
testNTLmul_LDADD = ../src/libfaast.la
testStem_SOURCES = testStem.c++
testStem_LDADD = ../src/libfaast.la
testTraceFrob_SOURCES = testTraceFrob.c++
testTraceFrob_LDADD = ../src/libfaast.la

TESTS_ENVIRONMENT = $(SHELL) < /dev/null
TESTS = $(check_PROGRAMS)
