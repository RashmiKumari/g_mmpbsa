/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_PBCUTIL_PBC_H
#define GMX_PBCUTIL_PBC_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

/* Maximum number of combinations of single triclinic box vectors
 * required to shift atoms that are within a brick of the size of
 * the diagonal of the box to within the maximum cut-off distance.
 */
#define MAX_NTRICVEC 12

/*! \brief Structure containing info on periodic boundary conditions */
typedef struct t_pbc
{
    //! The PBC type
    int ePBC;
    //! Number of dimensions in which PBC is exerted
    int ndim_ePBC;
    /*! \brief Determines how to compute distance vectors.
     *
     *  Indicator of how to compute distance vectors, depending
     *  on PBC type (depends on ePBC and dimensions with(out) DD)
     *  and the box angles.
     */
    int ePBCDX;
    /*! \brief Used for selecting which dimensions to use in PBC.
     *
     *  In case of 1-D PBC this indicates which dimension is used,
     *  in case of 2-D PBC this indicates the opposite
     */
    int dim;
    //! The simulation box
    matrix box;
    //! The lengths of the diagonal of the full box
    rvec fbox_diag;
    //! Halve of the above
    rvec hbox_diag;
    //! Negative of the above
    rvec mhbox_diag;
    //! Maximum allowed cutoff squared for the box and PBC used
    real max_cutoff2;
    /*! \brief Number of triclinic shift vectors.
     *
     *  Number of triclinic shift vectors depends on the skewedness
     *  of the box, that is mostly on the angles. For triclinic boxes
     *  we first take the closest image along each Cartesian dimension
     *  independently. When the resulting distance^2 is larger than
     *  max_cutoff2, up to ntric_vec triclinic shift vectors need to
     *  be tried. Because of the restrictions imposed on the unit-cell
     *  by GROMACS, ntric_vec <= MAX_NTRICVEC = 12.
     */
    int ntric_vec;
    //! The triclinic shift vectors in grid cells. Internal use only.
    ivec tric_shift[MAX_NTRICVEC];
    //!  The triclinic shift vectors in length units
    rvec tric_vec[MAX_NTRICVEC];
} t_pbc;

#define TRICLINIC(box) ((box)[YY][XX] != 0 || (box)[ZZ][XX] != 0 || (box)[ZZ][YY] != 0)

#define NTRICIMG 14
#define NCUCVERT 24
#define NCUCEDGE 36

enum
{
    ecenterTRIC, /* 0.5*(a+b+c)                  */
    ecenterRECT, /* (0.5*a[x],0.5*b[y],0.5*c[z]) */
    ecenterZERO, /* (0,0,0)                      */
    ecenterDEF = ecenterTRIC
};

/*! \brief Initiate the periodic boundary condition algorithms.
 *
 * pbc_dx will not use pbc and return the normal difference vector
 * when one or more of the diagonal elements of box are zero.
 * When ePBC=-1, the type of pbc is guessed from the box matrix.
 * \param[in,out] pbc The pbc information structure
 * \param[in] ePBC The PBC identifier
 * \param[in] box  The box tensor
 */
void set_pbc(t_pbc* pbc, int ePBC, const matrix box);

/*! \brief Compute distance with PBC
 *
 * Calculate the correct distance vector from x2 to x1 and put it in dx.
 * set_pbc must be called before ever calling this routine.
 *
 * Note that for triclinic boxes that do not obey the GROMACS unit-cell
 * restrictions, pbc_dx and pbc_dx_aiuc will not correct for PBC.
 * \param[in,out] pbc The pbc information structure
 * \param[in]    x1  Coordinates for particle 1
 * \param[in]    x2  Coordinates for particle 2
 * \param[out]   dx  Distance vector
 */
void pbc_dx(const t_pbc* pbc, const rvec x1, const rvec x2, rvec dx);

/*! \brief Put atoms inside triclinic box
 *
 * This puts ALL atoms in the triclinic unit cell, centered around the
 * box center as calculated by calc_box_center.
 * \param[in]    ecenter The pbc center type
 * \param[in]    box     The simulation box
 * \param[in,out] x       The coordinates of the atoms
 */
void put_atoms_in_triclinic_unitcell(int ecenter, const matrix box, gmx::ArrayRef<gmx::RVec> x);

#endif
