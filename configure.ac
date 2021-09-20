AC_INIT([RASPA],[2.0])
AC_CONFIG_SRCDIR([src/main.c])
AM_INIT_AUTOMAKE

AC_SEARCH_LIBS(ddot_,blas)
AC_SEARCH_LIBS(dsyevd_,lapack,HAVE_LAPACK="yes")

if test "$HAVE_LAPACK" = yes; then
  AC_DEFINE(HAVE_LAPACK, 0, [Define if you have LAPACK library.])
fi

AC_SEARCH_LIBS(fftw_plan_dft_r2c_1d,fftw3,HAVE_FFTW3="yes")

if test "$HAVE_FFTW3" = yes; then
  AC_DEFINE(HAVE_FFTW3, 0, [Define if you have FFTW library.])
fi

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77

dnl Checks for header files.
m4_warn([obsolete],
[The preprocessor macro `STDC_HEADERS' is obsolete.
  Except in unusual embedded environments, you can safely include all
  ISO C90 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
AC_PROG_EGREP

AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([string.h])
AC_CHECK_HEADERS([fftw3.h])

LT_INIT
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_FILES([Makefile \
src/Makefile \
forcefield/Makefile \
forcefield/ExampleMOFsForceField/Makefile \
forcefield/ExampleZeolitesForceField/Makefile \
forcefield/ExampleMoleculeForceField/Makefile \
forcefield/Dubbeldam2007FlexibleIRMOF-1/Makefile \
forcefield/Dubbeldam2007FlexibleIRMOF-10/Makefile \
forcefield/Dubbeldam2007FlexibleIRMOF-16/Makefile \
forcefield/CoreShellSchroderSauer/Makefile \
forcefield/Nicholas/Makefile \
molecules/Makefile \
molecules/ExampleDefinitions/Makefile \
framework/Makefile \
framework/Dubbeldam2007FlexibleIRMOF-1/Makefile \
framework/Dubbeldam2007FlexibleIRMOF-10/Makefile \
framework/Dubbeldam2007FlexibleIRMOF-16/Makefile \
framework/CoreShellSchroderSauer/Makefile \
framework/Demontis/Makefile \
framework/DemontisModified/Makefile \
framework/Nicholas/Makefile \
framework/NicholasModified/Makefile \
structures/Makefile \
structures/zeolites/Makefile \
structures/mofs/Makefile \
structures/minerals/Makefile \
structures/ceramics/Makefile \
structures/zeolites/cif/Makefile \
structures/zeolites/block/Makefile \
structures/zeolites/ions/Makefile \
structures/mofs/cif/Makefile  \
structures/mofs/block/Makefile \
structures/mofs/ions/Makefile \
structures/minerals/cif/Makefile  \
structures/ceramics/cif/Makefile  \
])
AC_OUTPUT
