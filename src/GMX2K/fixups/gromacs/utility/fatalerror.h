/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \file
 * \brief
 * Declares fatal error handling and debugging routines for C code.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FATALERROR_H
#define GMX_UTILITY_FATALERROR_H

#include <stdarg.h>
#include <stdio.h>

/*! \brief
 * Debug log file.
 *
 * Functions can write to this file for debug info.
 * Before writing to it, it should be checked whether the file is not NULL:
 * \code
   if (debug)
   {
       fprintf(debug, "%s", "Debug text");
   }
   \endcode
 */
extern FILE* debug;

/*! \brief
 * Fatal error reporting routine for \Gromacs.
 *
 * This function prints a fatal error message with a header that contains the
 * source file and line number of the call, followed by the string specified by
 * \p fmt and supplied parameters.
 * If \p fatal_errno is 0, only the message and arguments are printed.
 * If \p fatal_errno is a legal system errno or -1, a perror()-like message is
 * printed after the first message; if fatal_errno is -1, the last system errno
 * will be used.
 * The format of \p fmt uses printf()-like formatting.
 *
 * In case all MPI processes want to stop with the same fatal error,
 * use gmx_fatal_collective(), declared in network.h,
 * to avoid having as many error messages as processes.
 *
 * The first three parameters can be provided through ::FARGS:
 * \code
   gmx_fatal(FARGS, fmt, ...);
   \endcode
 */
[[noreturn]] void gmx_fatal(int fatal_errno, const char* file, int line, /*gmx_fmtstr*/ const char* fmt, ...);
/** Helper macro to pass first three parameters to gmx_fatal(). */
#define FARGS 0, __FILE__, __LINE__

/*! \brief
 * Implementation for range_check() and range_check_mesg().
 *
 * \p warn_str can be NULL.
 */
void _range_check(int n, int n_min, int n_max, const char* warn_str, const char* var, const char* file, int line);

/*! \brief
 * Checks that a variable is within a range.
 *
 * If \p n is not in range [n_min, n_max), a fatal error is raised.
 * \p n_min is inclusive, but \p n_max is not.
 */
#define range_check_mesg(n, n_min, n_max, str) \
    _range_check(n, n_min, n_max, str, #n, __FILE__, __LINE__)

/*! \brief
 * Checks that a variable is within a range.
 *
 * This works as range_check_mesg(), but with a default error message.
 */
#define range_check(n, n_min, n_max) _range_check(n, n_min, n_max, NULL, #n, __FILE__, __LINE__)

#endif
