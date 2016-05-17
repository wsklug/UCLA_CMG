#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([voom], [1.0], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE([foreign])
##AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_RANLIB


# Checks for libraries.
dnl AC_CHECK_LIB(efence, malloc)
VOOM_CHECK_BLITZ
VOOM_CHECK_TVMET

# Checks for header files.
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

AC_SUBST(blitz_includes)
AC_SUBST(blitz_libraries)
AC_SUBST(tvmet_includes)

# Checks for library functions.
AC_OUTPUT([Makefile])

