/**
 * @defgroup  MMEN Routines for Molecular-Mechanic Energy
 * @defgroup  RADIUS Routines for radius types
 * @defgroup  PBSA_PREP Routines for PBSA preparation
 * @defgroup PBSA_INPUT Routines to read input keywords file (mdp)
 * @defgroup PSIZE Routines to calculate grid dimensions
 */

/**
 * @file GMX46/g_mmpbsa.h
 * @brief Routines used in g_mmpbsa
 * @author Rashmi Kumari, Rajendra Kumar and Andrew Lynn
 * @ingroup MMEN RADIUS PBSA_PREP PBSA_INPUT PSIZE
 * @attention
 * @verbatim
 *
 *
 * This file is part of g_mmpbsa.
 *
 * Authors: Rashmi Kumari and Andrew Lynn
 * Contribution: Rajendra Kumar
 *
 * Copyright (C) 2013-2016 Rashmi Kumari and Andrew Lynn
 *
 * g_mmpbsa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_mmpbsa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 *
 */


#include <stdlib.h>
#include <math.h>
#include "gromacs/statutil.h"
#include "gromacs/typedefs.h"
#include "gromacs/smalloc.h"
#include "gromacs/vec.h"
#include "gromacs/tpxio.h"
#include "gromacs/rmpbc.h"
#include "gromacs/xvgr.h"


#ifndef G_MMPBSA_H_
#define G_MMPBSA_H_

/**
* @brief Data structure to hold information about non-bonded interactions
* @ingroup MMEN
* @author Rashmi Kumari
*/
typedef struct	{
	int nr_nb; /**< Number of non-bonded pairs excluding 1-2 and 1-3 pairs */
	int nr_14; /**<  Number of 1-4 pairs */
	int **pairNB;  /**< 2D array of of size `(nr_nb, 2)` containing all atom index of non-bonded pairs */
	int **pair14;  /**< 2D array of size `(nr_14, 2)` containing all 1-4 pairs atom index */
	int *pairtype; /**< array of pair-types for 1-4 Lenard-Jones parameters */
	gmx_bool *bItsA14; /**< To classify that first Atom of pairs is from subunit A in 1-4 pair list */
	gmx_bool *bItsA;   /**< To classify that first Atom of pairs is from subunit A in non-bonded pair list */
	gmx_bool *bItsB14; /**< To classify that second Atom of pairs is from subunit B in 1-4 pair list */
	gmx_bool *bItsB;   /**< To classify that second Atom of pairs is from subunit B in non-bonded pair list */

} t_non_bonded;

/**
* @brief Data structure to hold information about polar solvation energy calculations
* @ingroup PBSA_PREP
* @author Rashmi Kumari
*/
typedef struct {
	real cfac;            /**< Factor by which molecular dimensions should expand to get a coarse grid dimensions */
	real gridspace;       /**< Fine grid spacing in Angstrom */
	real gmemceil;        /**<  Maximum memory (MB) available*/
	real fadd;            /**<  The amount (in Å) to add to molecular dimensions to get a fine grid dimensions */
	real ofrac;           /**<  Used in mg-para: Overlap in mesh during parallel calculations */
	gmx_bool bParallel;   /**<  Whether parallel calculation is required */
	gmx_bool bFocus;      /**<  Whether focus type is enabled */

	int mg_type;          /**<  multi-grid calculation type */
	ivec dime;            /**<  grid dimenstions */
	ivec pdime;           /**<  grid dimenstions for parallel calculations */
	rvec cglen;           /**<  coarse grid lengths */
	rvec cgcent;          /**<  center of coarse-grid box */
	rvec fglen;           /**<  fine grid lengths */
	rvec fgcent;          /**<  center of fine grid box */
	int pbsolver;         /**<  PB Solver type: linear or non-linear */
	int bcfl;             /**<  Type of boundary conditions */
	real pcharge;         /**<  Magnitude of positive charge */
	real ncharge;         /**<  Magnitude of negative charge */
	real prad;            /**<  Radius of positive charged ion */
	real nrad;            /**<  Radius of negative charged ion */
	real pconc;           /**<  Concentration of positive charge */
	real nconc;           /**<  Concentration of negative charge */
	real pdie;            /**<  Solute dielectric constant */
	real sdie;            /**<  Solvent dielectric constant */
	real vdie;            /**<  Vacuum or reference dielectric constant*/
	int srfm;             /**<  To construct the dielectric and ion-accessibility coefficients */
	int chgm;             /**<  To map the biomolecular point charges to the grid */
	real sdens;           /**<  Number of grid points per Å^2 for constructing the molecular surface or solvent accessible surface */
	real srad;            /**<  Radius (in Å) of solvent molecules */
	real swin;            /**<  Value for cubic spline window for spline-based surface definitions */
	real temp;            /**<  Temperature in K */

} t_PolKey;

/**
* @brief Data structure to hold information about non-polar solvation energy calculations
* @ingroup PBSA_PREP
* @author Rashmi Kumari
*/
typedef struct	{
	gmx_bool    bWCA;   /**<  Whether to do WCA calculation*/
	real   bconc;       /**<  Bulk solvent density */
	real   wcarad;      /**<  Solvent radius */
	real   dpos;        /**<  Displacement in Å of the atomic positions for surface area  */
	real   sdens;       /**<  Number of grid points per Å^2 for constructing the molecular surface or solvent accessible surface */
	real   press;       /**<  Solvent pressure proportionality term of SAV model */
	real   savconst;    /**<  Offset or constant of SAV model */
	real   savrad;      /**<  Solvent radius for SAV */
	real   gamma;       /**<  Surface tension proportionality term of SASA model */
	real   sasaconst;   /**<  Offset or constant of SASA mode */
	real   sasrad;      /**<  Solvent radius for SASA */
	rvec   grid;        /**<  The quadrature grid spacing in Å for volume integral calculations */
	int    srfm;        /**<  To construct the solvent-related surface and volume */
	real   swin;        /**<  Size of spline window in Å for spline-based surface definitions */
	real   temp;        /**<  Temperature */
} t_APolKey;


/**
* @brief index for radius values in 2D array
* @ingroup PBSA_PREP
* @author Rashmi Kumari
*/
enum { eSASRAD, /**< For SASA model */
	     eSAVRAD, /**< For SAV model */
			 eWCARAD, /**< For WCA model */
			};

/**
* @brief index for PB solver
* @ingroup PBSA_INPUT
* @author Rashmi Kumari
*/
enum {eLPBE, /**< Linear solver */
	    eNPBE  /**< Non-linear solver */
		 };
static const char *PBsolver[] = { "lpbe", "npbe", NULL };

/**
* @brief index of keywords to boundary condition keywords
* @ingroup PBSA_INPUT
* @author Rashmi Kumari
*/
enum {ezero, esdh, emdh};
static const char *bcfl_words[] = { "zero", "sdh", "mdh", NULL };

/**
* @brief index of keywords to method for mapping the biomolecular point charges to the grid
* @ingroup PBSA_INPUT
* @author Rashmi Kumari
*/
enum { espl2, espl4, emol, esmol};
static const char *srfm_words[] = { "spl2", "spl4",  "mol", "smol", NULL};

enum {espl0 = 2};
static const char *chgm_words[] = {"spl2", "spl4", "spl0", NULL};

enum { sacc };
static const char *APsrfm_words[] = { "sacc", NULL};

enum { mg_auto, mg_para };
static const char *mg_words[] = { "mg-auto", "mg-para", NULL};


// energy_mm.c
/**
* @brief It generates all non-bonded atom-pair list from topology
* @ingroup MMEN
* @author Rashmi Kumari
* @param param Data structure \link t_non_bonded \endlink to hold information related to non-bonded pairs
* @param top Data structure t_topology containing all force-field parameters
* @param index A 1D array containing index of atoms (index group)
* @param isize Number of atoms in the input index
* @param splitIndex Atom index at which two groups are seperated
* @param bDiff TRUE or FALSE depends on user input for difference calculations
*/
int	energy_pair(t_non_bonded *param, t_topology *top, atom_id *index, int isize, int *splitIndex, gmx_bool bDiff);

/**
* @brief Data structure to hold information about non-bonded interactions
* @ingroup MMEN
* @author Rashmi Kumari
* @param[in] x 3D array containing atom coordinates
* @param[in] top Data structure t_topology containing all force-field parameters
* @param[in] param Data structure \link t_non_bonded \endlink to get information of non-bonded pairs
* @param[in] pdie Solute dielectric constant
* @param[in] bDiff TRUE or FALSE depends on user input for difference calculations
* @param[in] bDCOMP TRUE or FALSE depends on user input for energy decomposition calculations
* @param[out] EE 1D array containing electrostatic energy
* @param[out] Vdw 1D array containing van der Waals energy
*/
void Vac_MM(rvec *x, t_topology *top, t_non_bonded param, real pdie, gmx_bool bDiff, gmx_bool bDCOMP, double *EE, double *Vdw);

/**
* @brief Data structure to hold information about non-bonded interactions
* @ingroup MMEN
* @author Rashmi Kumari
* @param[in] x 3D array containing atom coordinates
* @param[in] top Data structure t_topology containing all force-field parameters
* @param[in] indexA 1D array containing atom index of subunit A
* @param[in] isizeA Size of `indexA` or number of atoms in subunit A
* @param[in] indexA 1D array containing atom index of subunit B
* @param[in] isizeB Size of `indexB` or number of atoms in subunit B
* @param[in] pdie Solute dielectric constant
* @param[in] bDCOMP TRUE or FALSE depends on user input for energy decomposition calculations
* @param[out] EE 1D array containing electrostatic energy
* @param[out] Vdw 1D array containing van der Waals energy
*/
void Vac_MM_without_14(rvec *x, t_topology *top, atom_id *indexA, int isizeA, atom_id *indexB, int isizeB, real pdie, gmx_bool bDCOMP, double *EE, double *Vdw);

// radius.c
/**
* @brief Data structure to contain information about atom properties: atomtype and its radius
* @ingroup RADIUS
* @author Rashmi Kumari
*/
typedef struct {
	int n;            /**< Number of atom-types */
	char **atomtype;  /**< List of atom-types */
	float *value;     /**< Radius values for respective atom-types */
} t_AtomProp;

/**
* @brief To contain radius-type options from which one type will be chosen by user
* @ingroup RADIUS
* @author Rashmi Kumari
*/
enum { eBondi,   /**< Bondi radius type */
	     eMbondi, /**< Modified Bondi radius type */
	     eMbondi2, /**< Modified-2 Bondi radius type */
			 eFF /**< Radius from force-field parameters */
			};

void substring(int start, int stop, const char *src, char *dest);

/**
* @brief To get radius of a particular atom-type
* @ingroup RADIUS
* @author Rashmi Kumari
* @param atomtype atom-type
* @param radtype Atom properties \link t_AtomProp \endlink for the selected radius type
* @param rvdw Default van der Waals radius for atom-types which are not present in \link t_AtomProp \endlink
* @returns Radius of the input atom-type
*/
real GetRad (char *atomtype, t_AtomProp *radtype, real rvdw);

/**
* @brief To store Bondi radii and atom-types in \link t_AtomProp \endlink
* @ingroup RADIUS
* @author Rashmi Kumari
* @param[in] radtype Atom properties \link t_AtomProp \endlink for the selected radius type
*/
void Bondi(t_AtomProp *bondi);

/**
* @brief To store Modified bondi radii and atom-types in \link t_AtomProp \endlink
* @ingroup RADIUS
* @author Rashmi Kumari
* @param[in] radtype Atom properties \link t_AtomProp \endlink for the selected radius type
*/
void mBondi(t_AtomProp *mbondi);

/**
* @brief To store Modified-2 Bondi radii and atom-types in \link t_AtomProp \endlink
* @ingroup RADIUS
* @author Rashmi Kumari
* @param[in] radtype Atom properties \link t_AtomProp \endlink for the selected radius type
*/
void mBondi2 (t_AtomProp *mbondi2);


// PbsaPrep.c
/**
* @brief To assign charge and radius of atoms and to store in t_topology data structure
* @ingroup PBSA_PREP
* @author Rashmi Kumari
* @param top Data structure t_topology containing all force-field parameters
* @param index A 1D array containing index of atoms (index group)
* @param isize Number of atoms in the input index
* @param eRadType \link eBondi \endlink, \link eMbondi \endlink, \link eMbondi2 \endlink, \link eFF \endlink
* @param[in] radtype Atom properties \link t_AtomProp \endlink for the selected radius type
* @param rvdw Default van der Waals radius for atom-types which are not present in \link t_AtomProp \endlink
* @return If successful 0
*/
int assignQR (t_topology *top, atom_id *index, int isize, int eRadType, t_AtomProp *radtype, real rvdw);

/**
* @brief To modify residue name. Residue may have same name in topology but have different structure.
* @ingroup PBSA_PREP
* @author Rashmi Kumari
* @param[in] top Data structure t_topology containing all force-field parameters
* @param[out] modresname List of unique residue name for all atoms
*/
void modResname (t_topology *top, char **modresname);

/**
* @brief To write a PQR file
* @ingroup PBSA_PREP
* @author Rashmi Kumari
* @param top Data structure t_topology containing all force-field parameters
* @param index A 1D array containing index of atoms (index group)
* @param isize Number of atoms in the input index
* @param ePBC PBC type
* @param matrix PBC Box dimension
* @param x 3D array containing atom coordinates
* @param modresname List of unique residue name for all atoms
* @param fnPQR Name of PQR file
* @return If successful 0
*/
int makePQR (t_topology *top, atom_id *index, int isize, int ePBC, matrix box, rvec *x, char **modresname, char *fnPQR);

/**
* @brief To assign radius of atoms for SASA and SAV calculation
* @ingroup PBSA_PREP
* @author Rashmi Kumari
* @param top Data structure t_topology containing all force-field parameters
* @param isize Number of atoms in the input index
* @param index A 1D array containing index of atoms (index group)
* @param[out] radius Pointer to a 2D array containing radius values. `radius[0]` for SASA, `radius[1]` for SAV and `radius[2]` for WCA calculations
* @param APolKey Data structure of t_APolKey type containing non-polar solvation input parameters
*/
void SasvRad(t_topology *top, int isize, atom_id *index, real ***radius, t_APolKey APolKey);

/**
* @brief To write a parameter file required by APBS for WCA calculation
* @ingroup PBSA_PREP
* @author Rashmi Kumari
* @param top Data structure t_topology containing all force-field parameters
* @param index A 1D array containing index of atoms (index group)
* @param isize Number of atoms in the input index
* @param modresname List of unique residue name for all atoms
* @param fnApbsParamAPol Name of output parameter file
*/
void ApbsParamAPol(t_topology *top, atom_id *index, int isize, char **modresname, char *fnApbsParamAPol);

/**
* @brief To write an input file required by APBS for polar solvation calculation
* @ingroup PBSA_PREP
* @author Rashmi Kumari
* @param param Data structure of t_PolKey type containing polar solvation input parameters
* @param fnPQR Name of PQR file
* @param fnPolAPBS Name of APBS input file
* @param bDECOMP TRUE or FALSE depends on user input for energy decomposition calculations
*/
int polarInAPBS(t_PolKey *param, char *fnPQR, char *fnPolAPBS, gmx_bool bDECOMP);

/**
* @brief To write an input file required by APBS for WCA non-polar solvation calculation
* @ingroup PBSA_PREP
* @author Rashmi Kumari
* @param param Data structure of t_APolKey type containing non-polar solvation input parameters
* @param fnPQR Name of PQR file
* @param fnAPolAPBS Name of APBS input file
* @param fnApbsParamAPol Name of input parameter file required by APBS
*/
int APolarInAPBS(t_APolKey *APolKey, char *fnPQR, char *fnAPolAPBS, char *fnApbsParamAPol);

/////////////////////////////////////////////////////////////////////////
////////This section is copied from nsc.h////////////////////////////////
////////Need for the surface area  and volume calculation////////////////
////////This section needs "gmxana" library /////////////////////////////
/////////////////////////////////////////////////////////////////////////
#define FLAG_DOTS       01
#define FLAG_VOLUME     02
#define FLAG_ATOM_AREA  04
extern int Mod_nsc_dclm_pbc(rvec *coords, real *radius, int nat,
			int  densit, int mode,
			real *value_of_area, double **at_area,
			real *value_of_vol, double ** at_vol,
			real **lidots, int *nu_dots,
			atom_id index[],int ePBC,matrix box);
////////////////////////////////////////////////////////////////////////


// InputPBSA.c
/**
* @brief To read input keyword file provided by user and store information in t_PolKey and t_APolKey
* @ingroup PBSA_INPUT
* @author Rashmi Kumari
* @param fnMdp Input mdp file name
* @param bPolar Whether polar solvation energy calculation will be done
* @param bAPolar Whether non-polar solvation energy calculation will be done
* @param fnApbsParamAPol Name of input parameter file required by APBS
* @param PolKey Data structure of t_PolKey type containing polar solvation input parameters
* @param ApolKey Data structure of t_APolKey type containing non-polar solvation input parameters
*/
int ReadInput	(char *fnMdp, gmx_bool *bPolar, gmx_bool *bAPolar, t_PolKey *PolKey, t_APolKey *ApolKey);

// psize.c
/**
* @brief Calculates fine and coarse grid dimensions with its center
* @ingroup PSIZE
* @author Rajendra Kumar
* @param top Data structure t_topology containing all force-field parameters
* @param index A 1D array containing index of atoms (index group)
* @param isize Number of atoms in the input index
* @param x 3D array containing atom coordinates
* @param PolKey Data structure of t_PolKey type containing polar solvation input parameters
* @param bCG [TRUE or FALSE]. If TRUE, new coarse-grid dimension and center will be copied in t_PolKey. If FALSE, previously assigned coarse-grid dimension and center will remain unchanged.
* @param bFocus If user has enabled the focusing on the specfic region of molecule.
* @return If successful 0
*/
int psize (t_topology *top, atom_id *index, int isize, rvec *x, t_PolKey *param, gmx_bool bCG, gmx_bool bFocus);

// apbs_main.c
#ifdef INT_APBS
int apbs( int argc,  char **argv, char *input_path, double *PolarEnergy, double *APolarEnergy, double *AtomEnergy);
#else
int ext_apbs(int isize, gmx_bool bVerbose, char *fnApbsOut, char *fnPolAPBS, double *PolarEnergy, double *APolarEnergy, double *AtomEnergy);
#endif

#endif /* G_MMPBSA_H_ */
