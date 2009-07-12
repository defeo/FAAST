AC_INIT([FAAST], [0.2.0])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CANONICAL_TARGET
AC_LANG([C++])

AM_INIT_AUTOMAKE([-Wall -Werror foreign dist-bzip2])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_CHECK_PROG([doxy], [doxygen], [yes], [no])
if test x$doxy != xyes ; then
   AC_MSG_WARN([Cannot find doxygen. You won't be able to generate the documentation.])
fi
AM_CONDITIONAL([DOXY], [test x$doxy = xyes])

LT_INIT

# Checks for libraries.
AC_ARG_ENABLE([ntl-static], [AS_HELP_STRING([--enable-ntl-static], 
			[link -lgf2x -lgmp -lm (disabling this might speed up the compilation but may also break the linking process) @<:@default=yes@:>@])])
if test x$enable_ntl_static != xno ; then
   AC_CHECK_LIB([m], [cos])
   AC_CHECK_LIB([gf2x], [gf2x_mul])
   AC_CHECK_LIB([gmp], [printf])
fi
AC_CHECK_LIB([ntl], [printf])
if ! echo $LIBS | grep -q ntl ; then
   AC_MSG_ERROR([Cannot find NTL. FAAST requires NTL, install it first.])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([floor])

AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_FILES([doxy.conf Makefile include/Makefile src/Makefile
			   test/Makefile doxygen/Makefile examples/Makefile])
AC_OUTPUT
