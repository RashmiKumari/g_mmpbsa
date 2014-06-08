
AC_INIT(g_mmpbsa, version-1.1)
AC_CONFIG_HEADER(config.h)
AC_MSG_NOTICE([g_mmpbsa: A tool for MMPBSA with GROMACS])
AM_INIT_AUTOMAKE([g_mmpbsa], [1.1])


AC_CANONICAL_HOST
dnl # LT_INIT([disable-shared])
AC_PROG_EGREP
AC_PROG_LN_S


cclist="gcc"
f77list="gfortran g95 g77"

AC_PROG_CC(${cclist})
AC_PROG_F77(${f77list})

dnl # -------------------------------------------------------------------------
dnl # Test to see if the compiler combination is valid for gcc/g77 or
dnl # gcc/gfortran/g95. Since code compiled with gcc 4.0 is only link
dnl # compatible with gfortran/g95. If they are using g77 then they need to
dnl # use gcc < 4.0 - D. Gohara 8/30/2006
dnl # -------------------------------------------------------------------------
if test ["$CC" = "gcc"]; then
	gcc_version=`gcc --version`
	gcc_complete_version=`echo ${gcc_version} | awk '{printf("%s",$3)}'`
	gcc_major_version=`echo ${gcc_complete_version} | awk '{printf("%i",$1)}'`

	if test ["$F77" = "g77"]; then
		if test ["${gcc_major_version}" -ge 4]; then
			AC_MSG_ERROR([Error, your GCC version is: ${gcc_complete_version} which is not compatible])
			AC_MSG_ERROR([with g77. You need to either use gcc < 4.0 as your C/CC/CXX compiler])
			AC_MSG_ERROR([or gfortran/g95 as your Fortran compiler ... exiting.])
		fi
	fi
fi

# ------------
# Math library
# ------------
AC_CHECK_LIB(m,sqrt,,AC_MSG_ERROR([No math library found]))



dnl #----------------------------------------------
dnl #To Check the Path of GROMACS Headers Directory
dnl #----------------------------------------------
AC_MSG_CHECKING([for --with-gmx-include option])
AC_ARG_WITH(
                [gmx-include],
                AS_HELP_STRING([--with-gmx-include=PATH], [full path to the GROMACS Include Directory]),
                [gmx_include_path=$withval]
        )
if test x"$gmx_include_path" != "xno" -a -n "$gmx_include_path"; then
	AC_MSG_RESULT([${gmx_include_path}])
else
	AC_MSG_RESULT([no])
	AC_MSG_NOTICE([WARNING: Path to GROMACS Inculde directory not provided])
	AC_MSG_NOTICE([checking for GROMACS Inlcude Directory in standard location])
	AC_CHECK_FILE(/usr/local/include/gromacs,[a="yes"],[a="no"])
	if test x"$a" = "xyes"; then
        	gmx_include_path="/usr/local/include"
	fi

	if test x"$a" = "xno" ; then
        	AC_CHECK_FILE(${HOME}/include/gromacs,[b="yes"],[b="no"])
	fi

	if test x"$b" = "xyes"; then
        	gmx_include_path="${HOME}/include"
	fi

	if test [ x"$a" = "xno" ] -a  [ x"$b" = "xno" ] ; then
		AC_MSG_ERROR([Please Specify the path to GROMACS Inlcude Directory
                e.g. --with-gmx-include=/path/to/include])
        fi
fi

gmx_include_path=${gmx_include_path}/gromacs
AC_CHECK_FILES([${gmx_include_path}/statutil.h 
		${gmx_include_path}/typedefs.h
		${gmx_include_path}/smalloc.h
		${gmx_include_path}/copyrite.h
		${gmx_include_path}/vec.h
		${gmx_include_path}/tpxio.h
		${gmx_include_path}/rmpbc.h
		${gmx_include_path}/xvgr.h
		${gmx_include_path}/readinp.h
		${gmx_include_path}/warninp.h
		${gmx_include_path}/string2.h ],[a="yes"],[a=""])

if test -z "$a"; then
	AC_MSG_ERROR([Could not find GROMACS Header files in ${gmx_include_path}])
fi

CPPFLAGS="$CPPFLAGS -I${gmx_include_path}"


dnl ------------------------------------
dnl To check the Path of GROMACS Library
dnl -----------------------------------
AC_MSG_CHECKING([for --with-gmx-lib option])
AC_ARG_WITH(
                [gmx-lib],
                AS_HELP_STRING([--with-gmx-lib=PATH], [full path to the GROMACS Library Directory]),
                [gmx_lib_path=$withval]
        )

if test x"$gmx_lib_path" != "xno" -a -n "$gmx_lib_path"; then
	AC_MSG_RESULT([${gmx_lib_path}])

dnl       libgmx.a
	AC_CHECK_FILE(${gmx_lib_path}/libgmx.a,[a="yes"],[a=""])
	if test -z "$a"; then
        	AC_MSG_ERROR([Could not find ${gromacs_lib_path}/libgmx.a library file])
	fi

dnl       libgmxana.a
	AC_CHECK_FILE(${gmx_lib_path}/libgmxana.a,[a="yes"],[a=""])
	if test -z "$a"; then
        	AC_MSG_ERROR([Could not find ${gromacs_lib_path}/libgmxana.a library file])
	fi
	
dnl       libmd.a
	AC_CHECK_FILE(${gmx_lib_path}/libmd.a,[a="yes"],[a=""])
	if test -z "$a"; then
        	AC_MSG_ERROR([Could not find ${gromacs_lib_path}/libmd.a library file])
	fi
else
	AC_MSG_RESULT([no])
	AC_MSG_NOTICE([WARNING: Path to GROMACS Libraries not provided])
	AC_MSG_NOTICE([checking for GROMACS Libraries in standard location])
	AC_CHECK_FILES([/usr/local/lib/libgmx.a
			/usr/local/lib/libgmxana.a
			/usr/local/lib/libmd.a ], [a="yes"],[a="no"])
	if test x"$a" = "xyes"; then
        	gmx_lib_path="/usr/local/lib"
	fi

	if test x"$a" = "xno" ; then
        	AC_CHECK_FILES([${HOME}/lib/libgmx.a
	                	${HOME}/lib/libgmxana.a
        	        	${HOME}/lib/libmd.a ], [b="yes"],[b="no"])
	fi

	if test x"$b" = "xyes"; then
        	gmx_lib_path="${HOME}/lib"
	fi

	if test [ x"$a" = "xno" ] -a  [ x"$b" = "xno" ] ; then
		AC_MSG_ERROR([GROMACS Libraries not found at the standard location. 
                Please provide the path e.g. --with-gmx-lib=/path/to/lib])
	fi
fi

dnl ---------------------------
dnl Checking for GROMACS version
dnl ----------------------------
AC_MSG_NOTICE([checking for GROMACS-4.5.x or GROMACS-4.6.x])
AC_MSG_CHECKING([for --enable-gmx46])
AC_ARG_ENABLE(
	[gmx46],
	AS_HELP_STRING([--enable-gmx46],[compile g_mmpbsa for GROMACS-4.6.x]),,
	[enable_gmx46=yes]
)
if test x"$enable_gmx46" == "xyes"; then
	AC_MSG_RESULT([yes])
	AC_MSG_NOTICE([g_mmpbsa will be compiled for GROMACS-4.6.x])
	
	AC_MSG_NOTICE([checking for fftw3f Libraries in standard location])
	AC_CHECK_FILE([/usr/local/lib/libfftw3f.a], [a="yes"],[a="no"])
	if test x"$a" = "xyes"; then
        	fftw_path="/usr/local/lib"
	fi

	if test x"$a" = "xno" ; then
        	AC_CHECK_FILE([${HOME}/lib/libfftw3f.a], [b="yes"],[b="no"])
	fi

	if test x"$b" = "xyes"; then
        	fftw_path="${HOME}/lib"
	fi
	
	if test x"$b" = "xno" ; then
        	AC_CHECK_FILE([/usr/lib/libfftw3f.a], [c="yes"],[c="no"])
	fi

	if test x"$c" = "xyes"; then
        	fftw_path="/usr/lib"
	fi

	if test [ x"$a" = "xno" ] -a  [ x"$b" = "xno" ] -a [ x"$c" = "xno" ] ; then
		AC_MSG_ERROR([WARNING: libfftw3f.a not found in standard location. 
                Please, provide path. e.g. export LDFLAGS=-L/path/to/lib])
		LIBS="$LIBS -pthread -lgmx -lgmxana -lmd -lfftw3f -ldl"
		LDFLAGS="$LDFLAGS -L${gmx_lib_path}"
	else
		LIBS="$LIBS -pthread -lgmx -lgmxana -lmd -lfftw3f -ldl"
		LDFLAGS="$LDFLAGS -L${fftw_path} -L${gmx_lib_path}"
	fi
	
else
	AC_MSG_NOTICE([g_mmpbsa will be compiled for GROMACS-4.5.x])
	LIBS="$LIBS -pthread -lgmx -lgmxana -lmd -ldl"
	LDFLAGS="$LDFLAGS -L${gmx_lib_path}"
fi

#lib_path_list="-L${gmx_lib_path}"
#lib_list="-lm -lgmx -lgmxana -lmd -ldl -lpthread"


# ---------------------------------------------------------------------------------
# To check ptherad library, require to make compaitable with gromacs installation
# ---------------------------------------------------------------------------------
#ACX_PTHREAD( , AC_MSG_ERROR([Cannot find POSIX threads - required for thread support using configure!]))
#ACX_PTHREAD([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#ACX_PTHREAD([[, b=""]])
#if test -z "$HAVE_PTHREAD"; then
#	LIBS="$PTHREAD_LIBS $LIBS"
#	CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
#	CC="$PTHREAD_CC "
#fi
	

#------------------------------------
# To check the Path of APBS Library
#-----------------------------------
AC_MSG_CHECKING([for --with-apbs-lib option])
AC_ARG_WITH(
		[apbs-lib],
		AS_HELP_STRING([--with-apbs-lib=PATH], [full path to the APBS Library Directory]),
		[apbs_lib_path=$withval]
	)
if test x"$apbs_lib_path" != "xno" -a -n "$apbs_lib_path"; then
	AC_MSG_RESULT([${apbs_lib_path}])
	
	#       libapbsmainroutines.a
        AC_CHECK_FILE(${apbs_lib_path}/libapbsmainroutines.a,[a="yes"],[a=""])
        if test -z "$a"; then
                AC_MSG_ERROR([Could not find ${apbs_lib_path}/libapbsmainroutines.a library file])
        fi

        #       libapbs.a
        AC_CHECK_FILE(${apbs_lib_path}/libapbs.a,[a="yes"],[a=""])
        if test -z "$a"; then
                AC_MSG_ERROR([Could not find ${apbs_lib_path}/libapbs.a library file])
        fi
else
	AC_MSG_RESULT([no])
	AC_MSG_NOTICE([WARNING: Path to APBS Libraries not provided])
	AC_MSG_NOTICE([checking for APBS Libraries in standard location])
	AC_CHECK_FILES([/usr/local/lib/libapbsmainroutines.a
                	/usr/local/lib/libapbs.a ], [a="yes"],[a="no"])
	if test x"$a" = "xyes"; then
        	apbs_lib_path="/usr/local/lib"
	fi

	if test x"$a" = "xno" ; then
        	AC_CHECK_FILES([${HOME}/lib/libapbsmainroutines.a
                        	${HOME}/lib/libapbs.a ], [b="yes"],[b="no"])
	fi

	if test x"$b" = "xyes"; then
        	apbs_lib_path="${HOME}/lib"
	fi

	if test [ x"$a" = "xno" ] -a [ x"$b" = "xno" ] ; then
        	AC_MSG_ERROR([APBS Libraries not found in standard location
		Please provide the path e.g. --with-apbs-lib=/path/to/lib
		which should contains libapbsmainroutines.a and libapbs.a files])
	fi
fi

LDFLAGS="$LDFLAGS -L${apbs_lib_path}"
LIBS="$LIBS -lapbsmainroutines -lapbs"

#lib_path_list="${lib_path_list} -L${apbs_lib_path}"
#lib_list="${lib_list} -lapbsmainroutines -lapbs"

#-------------------------------------------------
# To check the Path of APBS Source Code Directory
#-------------------------------------------------
AC_MSG_CHECKING([for --with-apbs option])
AC_ARG_WITH(
        [apbs],
        AS_HELP_STRING([--with-apbs=PATH], [full path to the APBS Source Code Directory]),
        [apbs_path=$withval]
)
if test x"$apbs_path" != "xno" -a -n "$apbs_path"; then
        AC_MSG_RESULT([${apbs_path}])
else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([Please Specify the path to APBS Source Code Directory])
fi

#	libapbsblas.a
AC_CHECK_FILE(${apbs_path}/contrib/lib/libapbsblas.a,[a="yes"],[a=""])
if test -z "$a"; then
        AC_MSG_ERROR([Could not find ${apbs_path}/contrib/lib/libapbsblas.a library file])
fi

#	libmaloc.a
AC_CHECK_FILE(${apbs_path}/contrib/lib/libmaloc.a,[a="yes"],[a=""])
if test -z "$a"; then
        AC_MSG_ERROR([Could not find ${apbs_path}/contrib/lib/libmaloc.a library file])
fi

#	libapbsgen.a
AC_CHECK_FILE(${apbs_path}/src/generic/.libs/libapbsgen.a,[a="yes"],[a=""])
if test -z "$a"; then
        AC_MSG_ERROR([Could not find ${apbs_path}/src/generic/.libs/libapbsgen.a library file])
fi

#	libz.a
AC_CHECK_FILE(${apbs_path}/contrib/zlib/.libs/libz.a,[a="yes"],[a=""])
if test -z "$a"; then
        AC_MSG_ERROR([Could not find ${apbs_path}/contrib/zlib/.libs/libz.a library file])
fi

LDFLAGS="$LDFLAGS -L${apbs_path}/contrib/lib -L${apbs_path}/src/generic/.libs -L${apbs_path}/contrib/zlib/.libs"
LIBS="$LIBS -lapbsblas -lmaloc -lapbsgen -lz"

#lib_path_list="${lib_path_list} -L${apbs_path}/contrib/lib -L${apbs_path}/src/generic/.libs -L${apbs_path}/contrib/zlib/.libs"
#lib_list="${lib_list} -lapbsblas -lmaloc -lapbsgen -lz"

#-------------------------------------------------
# To check the Path of APBS Source Code Directory
#-------------------------------------------------

AC_CHECK_FILES([${apbs_path}/src/aaa_inc/apbscfg.h
	        ${apbs_path}/src/aaa_inc/apbs/apbs.h
                ${apbs_path}/src/aaa_inc/apbs/vhal.h
	        ${apbs_path}/src/aaa_inc/apbs/nosh.h
	        ${apbs_path}/src/aaa_inc/apbs/vgrid.h
	        ${apbs_path}/src/aaa_inc/apbs/mgparm.h
	        ${apbs_path}/src/aaa_inc/apbs/pbeparm.h
	        ${apbs_path}/src/aaa_inc/apbs/femparm.h
		${apbs_path}/src/aaa_inc/apbs/vparam.h
                ${apbs_path}/src/aaa_inc/apbs/vacc.h],
              [a="yes"],[a=""])
if test -z "$a"; then
        AC_MSG_ERROR([Header files not found])
fi

AC_CHECK_FILE(${apbs_path}/contrib/include/maloc/maloc.h,[a="yes"],[a=""])
if test -z "$a"; then
        AC_MSG_ERROR([Could not find ${apbs_path}/contrib/include/maloc/maloc.h file])
fi

CPPFLAGS="$CPPFLAGS -I${apbs_path}/src/aaa_inc -I${apbs_path}/contrib/include"

dnl # -----------------------
dnl # CHECK WHETHER YOU WANT TO DISABLE OPENMP SUPPORT
dnl # -----------------------
AC_MSG_CHECKING([for --disable-openmp])
AC_ARG_ENABLE(
	[openmp],
	AS_HELP_STRING([--disable-openmp],[disables OpenMP compiler support]),,
	[disable_openmp=yes]
)
if test x"$disable_openmp" != "xyes"; then
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
	if [[ "${CC}" = "icc" ] && [ "${F77}" = "ifort" ]]; then
		FFLAGS="${FFLAGS} -openmp -assume cc_omp"
		LDFLAGS="${LDFLAGS} -static-intel"
	fi
	
	if [[ "${CC}" = "gcc" ] && [ "${F77}" = "gfortran" ]]; then
		gcc_version=`gcc --version`
		gcc_complete_version=`echo ${gcc_version} | awk '{printf("%s",$3)}' | awk -F'.' '{print $1$2}'`
		
		if test ["${gcc_complete_version}" -ge 42]; then
			FFLAGS="${FFLAGS} -fopenmp"
		else
		    AC_MSG_RESULT([WARNING])
			AC_MSG_RESULT([WARNING Version of gcc/gfortran does not support OpenMP. Disabling.])
			AC_MSG_RESULT([WARNING])
		fi
	fi
fi

##To add library path and libraies as a variable
#AC_SUBST(lib_path_list)
#AC_SUBST(lib_list)

## Extra CPP FLAGS
CPPFLAGS="$CPPFLAGS -c -std=gnu99 -fomit-frame-pointer -finline-functions -pthread -Wno-unused -msse2 -funroll-all-loops -fPIC -DPIC"


## MAKEFILES
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

## Configure output files
AC_OUTPUT

AC_MSG_NOTICE([
g_mmpbsa to be compiled with following settings:  ${CC} ${CPPFLAGS}
g_mmpbsa to be linked with following settings: ${F77} ${FFLAGS} ${LDFLAGS} ${LIBS}
])
