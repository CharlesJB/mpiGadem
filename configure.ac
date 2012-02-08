AC_PREREQ(2.62)
AC_INIT([rGADEM], [1.1.1], [http://www.rglab.org])

AC_CONFIG_SRCDIR([DESCRIPTION])

# determine RHOME
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

## pick up compiler as will be used by R CMD INSTALL/SHLIB
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`
CFLAGS=`"${RBIN}" CMD config CFLAGS`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP

# Pick up compiler flags to enable OpenMP
AC_OPENMP

dnl # Checks for libraries.
dnl # FIXME: Replace `main' with a function in `-lm':
dnl AC_CHECK_LIB([m], [pow])

# Checks for header files.
AC_SUBST(OPENMP_CFLAGS)
AC_ARG_ENABLE([parallel],
  AS_HELP_STRING([--enable-parallel@<:@=ARG@:>@],
    [Use GCD (ARG=gcd) or OpenMP (ARG=omp) to speed up calculation where possible (enabled by default). 
Not specifying ARG will choose the best method available on the platform. 
Use --disable-parallel (or ARG=no) to keep the code sequential even if GCD or OpenMP are present.]))

# AC_CHECK_PKG_OPENMP([true-action], [false-action])
# ---------------------
# Returns whether it is possible to use OpenMP from an R package
# As recommended by Simon Urbanek on R-Devel
AC_DEFUN([AC_CHECK_PKG_OPENMP], 
  [dnl
    AC_MSG_CHECKING([whether OpenMP will work in a package])
    ac_pkg_openmp=no
    AS_IF(
      [
        AS_IF([test -n "${OPENMP_CFLAGS}"],
          [
            AC_LANG_CONFTEST(
              [AC_LANG_PROGRAM(
                [[#include <omp.h>]],
                [[ return omp_get_num_threads (); ]])])
            PKG_CFLAGS="${OPENMP_CFLAGS}" PKG_LIBS="${OPENMP_CFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
          ])
        AC_MSG_RESULT([${ac_pkg_openmp}])
        test "x${ac_pkg_openmp}" = xyes], 
      [$1], 
      [$2])])


dnl
dnl If --enable_parallel=gcd: try GCD or fail
dnl
AS_IF([test "x${enable_parallel}" = xgcd],
  [AC_CHECK_HEADERS([dispatch/dispatch.h],
    :,
    [AC_MSG_FAILURE(
      [--enable_parallel=gcd was given, but the test for dispatch.h failed])
    ])
  ],
dnl
dnl If --enable_parallel=omp: try OpenMp or fail
dnl
  [AS_IF([test "x${enable_parallel}" = xomp],
    [AC_CHECK_PKG_OPENMP(
      :,
      [AC_MSG_FAILURE(
        [--enable_parallel=omp was given, but the test for OpenMP failed])
      ])
    ],
dnl
dnl If --enable_parallel=yes: first try GCD, then OMP, then fail
dnl
    [AS_IF([test "x${enable_parallel}" = xyes],
      [AC_CHECK_HEADERS([dispatch/dispatch.h],
        :,
        [AC_CHECK_PKG_OPENMP(
          :,
          [AC_MSG_FAILURE(
            [--enable_parallel was given, but the test for GCD and OpenMP failed])
          ])
        ])
      ],
dnl
dnl If --enable_parallel=no: use sequential code
dnl
      [AS_IF([test "x${enable_parallel}" = xno],
        [AC_MSG_NOTICE([no parallelism])],
dnl
dnl By default, first try GCD, then OpenMp, then use sequential code
dnl
        [AS_IF([test "x${enable_parallel}" = x],
          [AC_CHECK_HEADERS([dispatch/dispatch.h],
            :,
            [AC_CHECK_PKG_OPENMP])
          ],
          [AC_MSG_FAILURE(
            [Unrecognized argument to --enable_parallel: '${enable_parallel}'])
          ])
        ])
      ])
    ])
  ])

AS_IF([test "x${ac_pkg_openmp}" = xyes],
  [AC_DEFINE([HAVE_OPENMP])]
)
dnl AC_HEADER_STDC
dnl AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h])


dnl # Checks for typedefs, structures, and compiler characteristics.
dnl AC_C_CONST
dnl AC_TYPE_SIZE_T

dnl AC_CHECK_SIZEOF(long)

dnl AM_CONDITIONAL(LARGERINT, [test "$ac_cv_sizeof_long" = "8"])

dnl # Checks for library functions.
dnl AC_CHECK_FUNCS([pow fabs qsort floor rint sqrt])
dnl AC_FUNC_MALLOC
dnl AC_FUNC_REALLOC

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
