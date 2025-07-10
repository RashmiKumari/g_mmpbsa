/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
 * Declares common exception classes and macros for fatal error handling.
 *
 * The basic approach is the same as in boost::exception for storing additional
 * context information to exceptions, but since that functionality is a very
 * small and simple part of boost::exception, the code is duplicated here.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_EXCEPTIONS_H
#define GMX_UTILITY_EXCEPTIONS_H

#include <cstdio>
#include <cstdlib>

#include <exception>
#include <memory>
#include <string>
#include <type_traits>
#include <typeindex>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief
 * Helper function for terminating the program on an exception.
 *
 * \param[in] ex  Exception that is the cause for terminating the program.
 *
 * Does not throw, and does not return.
 */
[[noreturn]] void processExceptionAsFatalError(const std::exception& ex);

/*! \brief
 * Macro for catching exceptions at C++ -> C boundary.
 *
 * This macro is intended for uniform handling of exceptions when C++ code is
 * called from C code within Gromacs.  Since most existing code is written
 * using the assumption that fatal errors terminate the program, this macro
 * implements this behavior for exceptions.  It should only be used in cases
 * where the error cannot be propagated upwards using return values or such.
 *
 * Having this as a macro instead of having the same code in each place makes
 * it easy to 1) find all such locations in the code, and 2) change the exact
 * behavior if needed.
 *
 * Usage:
   \code
   try
   {
       // C++ code
   }
   GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
   \endcode
 */
#define GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR \
    catch (const std::exception& ex) { ::gmx::processExceptionAsFatalError(ex); }

//! \}

} // namespace gmx

#endif
