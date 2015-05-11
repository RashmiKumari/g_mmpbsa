/*
 * This file is part of g_mmpbsa.
 *
 * Authors: Rashmi Kumari and Andrew Lynn
 * Contribution: Rajendra Kumar
 *
 * Copyright (C) 2013, 2014, 2015 Rashmi Kumari and Andrew Lynn
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
 *
 */


#include <stdlib.h>
#include <math.h>
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
//#include "copyrite.h"
#include "vec.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "xvgr.h"


#ifndef G_MMPBSA_H_
#define G_MMPBSA_H_

//MM Energy
typedef struct
	{
    		int nr_nb, nr_14;
    		int **pairNB, **pair14;
    		int *pairtype;
    		gmx_bool *bItsA14, *bItsA, *bItsB14, *bItsB;
	} t_non_bonded;

//APBS Polar Input Keywords
typedef struct {
	real cfac;
	real gridspace;
	real gmemceil;
	real fadd;
	real ofrac;
	gmx_bool bParallel;
	gmx_bool bFocus;

	int mg_type;
	ivec dime;
	ivec pdime;
	rvec cglen;
	rvec cgcent;
	rvec fglen;
	rvec fgcent;
	int pbsolver;
	int bcfl;
	real pcharge;
	real ncharge;
	real prad;
	real nrad;
	real pconc;
	real nconc;
	real pdie;		//solute dielectric constant
	real sdie;		//Solvent dielectric constant
	real vdie;		//Vaccum dielectric constant
	int srfm;
	int chgm;
	real sdens;
	real srad;
	real swin;
	real temp;
} t_PolKey;

// Non-polar keywords
typedef struct          {
	gmx_bool    bWCA;
    real   bconc;
    real   wcarad;
    real   dpos;
    real   sdens;
    real   press;
    real   savconst;
    real   savrad;
    real   gamma;
    real   sasaconst;
    real   sasrad;
    rvec   grid;
    int    srfm;
    real   swin;
    real   temp;
} t_APolKey;

enum {eSASRAD, eSAVRAD, eWCARAD};

enum {eLPBE, eNPBE};
static const char *PBsolver[] = { "lpbe", "npbe", NULL };

enum {ezero, esdh, emdh};
static const char *bcfl_words[] = { "zero", "sdh", "mdh", NULL };

enum { espl2, espl4, emol, esmol};
static const char *srfm_words[] = { "spl2", "spl4",  "mol", "smol", NULL};

enum {espl0 = 2};
static const char *chgm_words[] = {"spl2", "spl4", "spl0", NULL};

enum { sacc };
static const char *APsrfm_words[] = { "sacc", NULL};

enum { mg_auto, mg_para };
static const char *mg_words[] = { "mg-auto", "mg-para", NULL};


// energy_mm.c
int	energy_pair(t_non_bonded *param, t_topology *top, atom_id *index, int isize, int *splitIndex, gmx_bool bDiff);

void Vac_MM(rvec *x, t_topology *top, t_non_bonded param, real pdie, gmx_bool bDiff, gmx_bool bDCOMP, double *EE, double *Vdw);

void Vac_MM_without_14(rvec *x, t_topology *top, atom_id *indexA, int isizeA, atom_id *indexB, int isizeB, real pdie, gmx_bool bDCOMP, double *EE, double *Vdw);

// radius.c
typedef struct {
	int n;
	char **atomtype;
	float *value;
} t_AtomProp;

enum { eBondi, eMbondi, eMbondi2, eFF };

void substring(int start, int stop, const char *src, char *dest);

real GetRad (char *atomtype, t_AtomProp *radtype, real rvdw);

void Bondi(t_AtomProp *bondi);

void mBondi(t_AtomProp *mbondi);

void mBondi2 (t_AtomProp *mbondi2);


// PbsaPrep.c
int assignQR (t_topology *top, atom_id *index, int isize, int eRadType, t_AtomProp *radtype, real rvdw);

void modResname (t_topology *top, char **modresname);

int makePQR (t_topology *top, atom_id *index, int isize, int ePBC, matrix box, rvec *x, char **modresname, char *fnPQR);

void SasvRad(t_topology *top, int isize, atom_id *index, real ***radius, t_APolKey APolKey);

void ApbsParamAPol(t_topology *top, atom_id *index, int isize, char **modresname, char *fnApbsParamAPol);

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
int ReadInput	(char *fnMdp, gmx_bool *bPolar, gmx_bool *bAPolar, t_PolKey *PolKey, t_APolKey *ApolKey);
int polarInAPBS(t_PolKey *param, char *fnPQR, char *fnPolAPBS, gmx_bool bDECOMP);
int APolarInAPBS(t_APolKey *APolKey, char *fnPQR, char *fnAPolAPBS, char *fnApbsParamAPol);

// psize.c
int psize (t_topology *top, atom_id *index, int isize, rvec *x, t_PolKey *param, gmx_bool bCG, gmx_bool bFocus);

// apbs_main.c
#ifdef INT_APBS
int apbs( int argc,  char **argv, char *input_path, double *PolarEnergy, double *APolarEnergy, double *AtomEnergy);
#else
int ext_apbs(int isize, gmx_bool bVerbose, char *fnApbsOut, char *fnPolAPBS, double *PolarEnergy, double *APolarEnergy, double *AtomEnergy);
#endif

#endif /* G_MMPBSA_H_ */
