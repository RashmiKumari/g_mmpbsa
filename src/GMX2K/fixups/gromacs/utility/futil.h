/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
 * Low-level wrappers for OS-specific file handling with some \Gromacs
 * customizations.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FUTIL_H
#define GMX_UTILITY_FUTIL_H

#include <cstdio>
#include <string>

/*! \brief
 * Opens a file, with \Gromacs-specific additions.
 *
 * If the file is in compressed format, opens a pipe which uncompresses the
 * file on the fly.  For this to work, gmx_ffclose() and frewind() should
 * always be used for files opened with gmx_ffopen() instead of fclose() and
 * rewind().  For compressed files, the \p file parameter should be passed
 * without the compressed extension; if that file is not found, then a few
 * compression extensions are tried.
 * Creates a backup if a file opened for writing already exists before
 * overwriting it.
 * A fatal error results if the file cannot be opened, for whatever reason.
 */
FILE* gmx_ffopen(const std::string& file, const char* mode);

/** Closes a file opened with gmx_ffopen(). */
int gmx_ffclose(FILE* fp);

/*! \brief
 * Creates unique name for temp file (wrapper around mkstemp).
 *
 * \p buf should be at least 7 bytes long
 */
void gmx_tmpnam(char* buf);

#endif
