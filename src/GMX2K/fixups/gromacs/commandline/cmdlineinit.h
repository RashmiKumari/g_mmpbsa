/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2018,2019, by the GROMACS development team, led by
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
 * Declares functions for initializing the \Gromacs library for command line use.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEINIT_H
#define GMX_COMMANDLINE_CMDLINEINIT_H

/*! \brief
 * Implements a main() method that runs a given C main function.
 *
 * \param argc         \c argc passed to main().
 * \param argv         \c argv passed to main().
 * \param mainFunction The main()-like method to wrap.
 *
 * This method creates a dummy command line module that does its
 * processing by calling \p mainFunction.  It then runs this module as with
 * gmx::runCommandLineModule().
 * This allows the resulting executable to handle common options and do
 * other common actions (e.g., startup headers) without duplicate code
 * in the main methods.
 *
 * \p mainFunction should call parse_common_args() to process its command-line
 * arguments.
 *
 * Usage:
 * \code
   int my_main(int argc, char *argv[])
   {
       // <...>
   }

   int main(int argc, char *argv[])
   {
       return gmx_run_cmain(argc, argv, &my_main);
   }
   \endcode
 *
 * Does not throw.  All exceptions are caught and handled internally.
 */
int gmx_run_cmain(int argc, char* argv[], int (*mainFunction)(int, char*[]));

#endif
