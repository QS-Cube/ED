#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([QS3], [1.0.0],[qsss.inquiry{at}gmail.com])
AC_CONFIG_SRCDIR([source/main.f90])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign])

# Checks for compilers
AC_PROG_CC

FCFLAGS_PRESERVE=$FCFLAGS
AC_LANG(Fortran)
AC_PROG_FC(ifort gfortran)
FCFLAGS=$FCFLAGS_PRESERVE
# Restore FCFLAGS, because AC_PROG_FC() sets FCFLAGS to '-g -O2', '-O2' or '-g' undesirably.

## FCFLAGS:
case "$FC" in
  gfortran*)
    test -n "$FCFLAGS" || FCFLAGS="-fopenmp -Ofast"
    ;;
  ifort)
   test -n "$FCFLAGS" || FCFLAGS="-O3 -ipo -no-prec-div -fp-model fast=2 -xHost -parallel -qopenmp"
    ;;
esac

## LAPACK 
AC_ARG_WITH(lapack,[AS_HELP_STRING([--with-lapack=ARG],[LAPACK library.
                                     ARG must be one of
                                     lapack (default),
                                     mkl (Intel MKL).])])
case "$with_lapack" in
  lapack | yes | "")
    AC_CHECK_LIB(blas,  dsymv,,[AC_CHECK_LIB(f77blas,  dsymv,,AC_MSG_ERROR([You need  BLAS  library. See config.log.]))])
    AC_CHECK_LIB(lapack,dsyevx,,AC_MSG_ERROR([You need LAPACK library. See config.log.]))
    ;;
  mkl)
    case "$FC" in
      gfortran*)
        case "$LIBS" in
     	  *mkl_gf_lp64*)
       	    AC_MSG_NOTICE([lapack... It seems that there is -lmkl_gf_lp64 already in LIBS.])
       	    ;;
      	  *)
     	    AC_MSG_NOTICE([lapack... Add -Wl,--start-group -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group -lgomp -lpthread -lm -ldl at the tail of LIBS.])
            LIBS="$LIBS -Wl,--start-group -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group -lgomp -lpthread -lm -ldl"
      	    ;;
        esac
        ;;
      ifort)
        case "$LIBS" in
          *mkl_intel_lp64*)
            AC_MSG_NOTICE([lapack... It seems that there is -lmkl_intel_lp64 already in LIBS.])
            ;;
          *)
            AC_MSG_NOTICE([lapack... Add -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -lgomp -lpthread -lm -ldl at the tail of LIBS.])
            LIBS="$LIBS -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -lgomp -lpthread -lm -ldl"
            ;;
        esac
        ;;
    esac
    ;;
  *)
    AC_MSG_ERROR([Unknown LAPACK library: $with_lapack])
    ;;
esac

AC_CONFIG_FILES([Makefile
          source/Makefile
  source_only_u1/Makefile
   source_DSFLan/Makefile])
AC_OUTPUT
