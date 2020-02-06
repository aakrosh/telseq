#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([TelSeq], [0.0.3], [zhihao.ding@gmail.com])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([Telseq/telseq.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for header files.
AC_LANG([C++])

AC_ARG_WITH(seqlib,
        AS_HELP_STRING([--with-seqlib=DIR], [specify directory containing the seqlib library]),
        [ if test -d "$withval" ; then
            bamtools_ldflags="-L$withval/lib"
    		bamtools_include="-I$withval -I$withval/fermi-lite -I$withval/htslib -I$withval/bwa -I$withval/SeqLib"
		else
          	AC_MSG_ERROR([Cannot stat directory $withval])
      	fi ]
        )

# Set compiler flags.
#AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Werror -Wno-unknown-pragmas")
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Wno-unknown-pragmas -g -O0")
#AC_SUBST(CXXFLAGS, "-O3")
#AC_SUBST(CFLAGS, "-O3")

AC_SUBST(CXXFLAGS, "-std=c++11")
AC_SUBST(CFLAGS, "-std=c++11")

AC_SUBST(CPPFLAGS, "$CPPFLAGS $bamtools_include $gzstream_include")
AC_SUBST(LDFLAGS, "$bamtools_ldflags $gzstream_ldflags $LDFLAGS")

AC_MSG_NOTICE(["CPPFLAGS is $CPPFLAGS"])
AC_MSG_NOTICE(["LDFLAGS is $LDFLAGS"])

# We always need to specify to link in seqlib
AC_SUBST(LIBS, "$LIBS $withval/lib/libhts.a $withval/lib/libbwa.a $withval/lib/libfml.a $withval/lib/libseqlib.a -lz -lbz2 -llzma")

#AC_MSG_NOTICE(["LIBS is $LIBS"])

# Make sure the seqlib headers can be found
AC_CHECK_HEADER([SeqLib/BamReader.h],,[AC_MSG_ERROR(
[The SeqLib library must be installed. \
You can specify its path with the --with-seqlib=PATH option])],-)

AC_CONFIG_FILES([Makefile
		Util/Makefile
    Telseq/Makefile])

AC_OUTPUT
