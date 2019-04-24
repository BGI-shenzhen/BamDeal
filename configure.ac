#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT(BamDeal, 0.22, hewm2008@qq.com)
AM_INIT_AUTOMAKE(BamDeal,0.22)
AC_CONFIG_SRCDIR([BamDeal.cpp])
AC_CONFIG_HEADER([config.h])
# Checks for programs
AC_PROG_AWK
AC_PROG_CXX([g++])
AC_PROG_CC([gcc])
AC_PROG_CPP
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main],,AC_MSG_ERROR([You need zlib >= 1.2.3 to build $PACKAGE]))

# FIXME: Replace `main' with a function in `-lhts':
AC_CHECK_LIB([hts], [main],,AC_MSG_ERROR([You need Pre-install Library libhts samtools/htslib >= 1.5 to build $PACKAGE]))



# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([zlib.h],[],[AC_MSG_ERROR([You need include htslib/zlib.h to build $PACKAGE])],[])
AC_CHECK_HEADERS([htslib/sam.h htslib/kseq.h zlib.h htslib/kstring.h stdlib.h string.h unistd.h])
AC_CHECK_HEADERS([htslib/kseq.h],[],[AC_MSG_ERROR([You need include htslib/kseq.h to build $PACKAGE])],[])
AC_CHECK_HEADERS([htslib/sam.h],[],[AC_MSG_ERROR([You need include htslib/sam.h to build $PACKAGE])],[])
AC_CHECK_HEADERS([htslib/kstring.h],[],[AC_MSG_ERROR([You need include htslib/kstring.h to build $PACKAGE])],[])


# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T



# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset regcomp sqrt strdup strncasecmp strtol])




AC_CONFIG_FILES([Makefile])
AC_OUTPUT
