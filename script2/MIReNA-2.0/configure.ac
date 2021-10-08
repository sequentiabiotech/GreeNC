#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([MIReNA], [1.0] [anthony.mathelier@gmail.com])
#AC_CONFIG_SRCDIR([asml/asml.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_MAKE_SET

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/file.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([sqrt strchr strdup strrchr strstr])

AC_CONFIG_FILES([Makefile
                 folding/Makefile
                 asml/Makefile])

AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_MAKE_SET
AM_PROG_CC_C_O

AC_OUTPUT
