/**
 * This file is part of g_mmpbsa.
 *
 * Authors: Rashmi Kumari and Andrew Lynn
 * Contribution: Rajendra Kumar
 *
 * Copyright (C) 2013-2015 Rashmi Kumari and Andrew Lynn
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
#include <string.h>

#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/rmpbc.h"
#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/xvgr.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/index.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/cmdlineinit.h"


#include "g_mmpbsa.h"

void CopyRightMsg()	{
	char *copyright[] =	{
			"                                                                        ",
			"                           :-)  g_mmpbsa (-:                            ",
			"                                                                        ",
			"               Authors: Rashmi Kumari and Andrew Lynn                   ",
			"               Contribution: Rajendra Kumar                             ",
			"                                                                        ",
			"       Copyright (C) 2013 - 2015 Rashmi Kumari and Andrew Lynn          ",
			"                                                                        ",
			"g_mmpbsa is free software: you can redistribute it and/or modify        ",
			"it under the terms of the GNU General Public License as published by    ",
			"the Free Software Foundation, either version 3 of the License, or       ",
			"(at your option) any later version.                                     ",
			"                                                                        ",
			"g_mmpbsa is distributed in the hope that it will be useful,             ",
			"but WITHOUT ANY WARRANTY; without even the implied warranty of          ",
			"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ",
			"GNU General Public License for more details.                            ",
			"                                                                        ",
			"You should have received a copy of the GNU General Public License       ",
			"along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.       ",
			"                                                                        ",
			"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     ",
			"\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ",
			"LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   ",
			"A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    ",
			"OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   ",
			"SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
			"TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  ",
			"PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ",
			"LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    ",
			"NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      ",
			"SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            ",
			"                                                                        ",
			"                           :-)  g_mmpbsa (-:                            ",
			"                                                                        ",
			"                                                                        "
	};
	int i = 0;
	char *str;
	for(i=0; i<34; i++)	{
		str = strdup(copyright[i]);
		fprintf(stderr,"%s\n", str);
	}
}


typedef struct {
  const char *key;
  const char *author;
  const char *title;
  const char *journal;
  int volume,year;
  const char *pages;
} t_citerec;

void show_citation(FILE *fp, const char *key)	{

	static const t_citerec citedb[] = {
	    { "Gromacs45",
	      "S Pronk, S Páll, R Schulz, P Larsson, P Bjelkmar, R Apostolov, M R Shirts, J C Smith, P M Kasson, D van der Spoel, B Hess, and E Lindahl",
	      "GROMACS 4.5: a high-throughput and highly parallel open source molecular simulation toolkit",
	      "Bioinformatics",
	      29, 2013, "845-854" },
	    { "APBS2001",
	      "N A Baker, D Sept, S Joseph, M J Holst, and J A McCammon",
	      "Electrostatics of nanosystems: Application to microtubules and the ribosome",
	      "Proc. Natl. Acad. Sci. USA",
	      98, 2001, "10037-10041" },
	    { "APBS2006",
	      "J A Wagoner and N A Baker",
	      "Assessing implicit models for nonpolar mean solvation forces: The importance of dispersion and volume terms",
	      "Proc. Natl. Acad. Sci. USA",
	      103, 2006, "8331-8336" },
	    { "Eisenhaber95",
	      "F Eisenhaber and P Lijnzaad and P Argos and C Sander and M Scharf",
	      "The Double Cube Lattice Method: Efficient Approaches to Numerical Integration of Surface Area and Volume and to Dot Surface Contouring of Molecular Assemblies",
	      "J. Comp. Chem.",
	      16, 1995, "273-284" },
	};

	#define NSTR (int)asize(citedb)

	int  j,index;
	char *author;
	char *title;
#define LINE_WIDTH 79

	if (fp == NULL)
    return;

	for(index=0; (index<NSTR) && (strcmp(citedb[index].key,key) != 0); index++)
    ;


  if (index < NSTR) {
    /* Insert newlines */
    author = wrap_lines(citedb[index].author,LINE_WIDTH,0,FALSE);
    title  = wrap_lines(citedb[index].title,LINE_WIDTH,0,FALSE);
    fprintf(fp,"%s\n%s\n%s %d (%d) pp. %s\n",
	    author,title,citedb[index].journal,
	    citedb[index].volume,citedb[index].year,
	    citedb[index].pages);
    sfree(author);
    sfree(title);
  	  }
  	  else {
  		  fprintf(fp,"Entry %s not found in citation database\n",key);
  	  }

  fflush(fp);
}


void cite(gmx_bool bPolar, gmx_bool bAPolar, t_APolKey APolarKeyWords){

	fprintf(stderr,"\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n\n");
	fprintf(stderr,"-------- -------- ------------------- -------- --------\n");
	show_citation(stderr, "Gromacs45");
	fprintf(stderr,"-------- -------- ------------------- -------- --------\n\n");
	if(bPolar)	{
		fprintf(stderr,"-------- -------- ------------------- -------- --------\n");
		show_citation(stderr, "APBS2001");
		fprintf(stderr,"-------- -------- ------------------- -------- --------\n\n");
	}

	if(bAPolar)	{
		if ((APolarKeyWords.gamma != 0) || (APolarKeyWords.gamma != 0))	{
			fprintf(stderr,"-------- -------- ------------------- -------- --------\n");
			show_citation(stderr, "Eisenhaber95");
			fprintf(stderr,"-------- -------- ------------------- -------- --------\n\n");
		}

		if (APolarKeyWords.bWCA)	{
			fprintf(stderr,"-------- -------- ------------------- -------- --------\n");
			show_citation(stderr, "APBS2006");
			fprintf(stderr,"-------- -------- ------------------- -------- --------\n\n");
		}
	}
	fprintf(stderr,"-------- -------- ------------------- -------- --------\n");
	fprintf(stderr,"g_mmpbsa—A GROMACS Tool for High-Throughput MM-PBSA Calculations.\n");
	fprintf(stderr,"Kumari R. et al. (2014)\n");
	fprintf(stderr,"J. Chem. Inf. Model., 54 (7), 1951–1962\n");
	fprintf(stderr,"URL: http://pubs.acs.org/doi/abs/10.1021/ci500020m \n");
	fprintf(stderr,"-------- -------- ------------------- -------- --------\n\n");

	fprintf(stderr,"-------- -------- --- Thank You --- -------- --------\n\n");
}

void decomp_calc_energy (t_topology top, int sizeA, atom_id *indexA, int sizeB, atom_id *indexB, int sizeAB, atom_id *indexAB,
		                   double *EnergyAtomA, double *EnergyAtomB, double *EnergyAtomAB, gmx_bool *bResA, gmx_bool *bResB, double *ResEnergy)	{

  int i = 0;
  double *ResOnly, *ResComplex;
  snew(ResOnly,top.atoms.nres);
  snew(ResComplex,top.atoms.nres);

  for(i=0;i<top.atoms.nres;i++){
	  ResEnergy[i] =0;
	  ResOnly[i]=0;
	  ResComplex[i]=0;
  }

   //Residues in Complex
  for(i=0;i<sizeAB;i++)
	  ResComplex[top.atoms.atom[indexAB[i]].resind] += EnergyAtomAB[i];

  //Residues in A only
  for(i=0;i<sizeA;i++)
	  ResOnly[top.atoms.atom[indexA[i]].resind] += EnergyAtomA[i];

  //Residues in B only
  for(i=0;i<sizeB;i++)
	  ResOnly[top.atoms.atom[indexB[i]].resind] += EnergyAtomB[i];

  for(i=0;i<top.atoms.nres;i++)
	  ResEnergy[i] = (ResComplex[i] - ResOnly[i]);

  sfree(ResOnly);
  sfree(ResComplex);
}

int gmx_do_mmpbsa(int argc, char *argv[]) {
  const char *desc[] =
  {	"g_mmpbsa calculates relative binding free energy using the MM-PBSA method for "
    "bio-molecular associations such as protein-protein, protein-ligand, protein-DNA etc. "
    "It calculates three components of the binding energy in separate files, so that",
    "user will have choice to calculate MM, PB and SA energy values according to their ",
	"objective. It also calculates contribution of each residue to the net binding energy",
	"and provides information about important contributing residues to the molecular",
	"association.\n\n",
	"For more detail, see please visit <http://rashmikumari.github.io/g_mmpbsa>"
  };

  /* Command-line arguments */
  static real rvdw = 0.1, pdie = 1 ;
  static int ndots = 24;
  static gmx_bool bDIFF = TRUE, bMM = TRUE, bDCOMP = FALSE, bPBSA = FALSE, bFocus = FALSE;
  gmx_bool bIncl14 = FALSE, bVerbose=FALSE;
  static const char *rtype[] = { NULL, "bondi", "mbondi", "mbondi2", "force-field", NULL };

  t_pargs pa[] =
        {
          { "-silent", FALSE,  etBOOL, {&bVerbose}, "Display messages, output and errors from external APBS program. Only works with external APBS program" },
          { "-rad",    FALSE,  etENUM, { rtype   }, "van der Waal radius type" },
          { "-rvdw",   FALSE,  etREAL, { &rvdw   }, "Default van der Waal radius (in nm) if radius not found for any atom-types" },
          { "-mme",    TRUE,   etBOOL, { &bMM    }, "To calculate vacuum molecular mechanics energy" },
          { "-pdie",   TRUE,   etREAL, { &pdie   }, "Dielectric constant of solute. Should be same as of polar solvation" },
          { "-incl_14",FALSE,  etBOOL, { &bIncl14}, "Include 1-4 atom-pairs, exclude 1-2 and 1-3 atom pairs during MM calculation. Should be \"yes\" when groups are bonded with each other." },
          { "-focus",  TRUE,   etBOOL, { &bFocus }, "To enable focusing on the specfic region of molecule, group of atoms must be provided in index file" },
          { "-pbsa",   FALSE,  etBOOL, { &bPBSA  }, "To calculate polar and/or non-polar solvation energy" },
          { "-ndots",  FALSE,  etINT,  { &ndots  }, "Number of dots per sphere in the calculation of SASA, more dots means more accuracy" },
          { "-diff",   TRUE,   etBOOL, { &bDIFF  }, "Calculate the energy difference between two group otherwise only calculates for one group" },
          { "-decomp", FALSE,  etBOOL, { &bDCOMP }, "Decomposition of energy for each residue" }
        };

  t_filenm fnm[] =
    {
      { efTRX, "-f", NULL, ffREAD },
      { efTPX, "-s", NULL, ffREAD },
      { efMDP, "-i", NULL, ffOPTRD },
      { efNDX, "-n", "index", ffOPTRD },
      { efXVG, "-mm", "energy_MM", ffOPTWR },
      { efXVG, "-pol", "polar", ffOPTWR },
      { efXVG, "-apol", "apolar", ffOPTWR },
      { efDAT, "-mmcon", "contrib_MM", ffOPTWR },
      { efDAT, "-pcon", "contrib_pol", ffOPTWR },
      { efDAT, "-apcon", "contrib_apol", ffOPTWR }
    };

  output_env_t oenv;
  int ngrps, nframe = 0;
  int g;
  int rc, i, j, k;
  gmx_bool *bAindex, *bBindex, *bABindex;

  //Variable for topology (read_tps_conf)
  char title[STRLEN], buf[256];
  t_topology top;
  int ePBC;
  rvec *xtop;
  matrix box;
  t_atoms *atoms;
  char *atomtype, *resname, *atomname, **modresname = NULL;
  double *tmpEE=NULL, *tmpVdw=NULL;
  real tmpArea = 0, tmpVolume = 0, Area = 0, Volume = 0, **radiusA = NULL, **radiusB = NULL, **radiusAB = NULL;
  double *TempAtomAPolA=NULL, *TempAtomAPolB=NULL, *TempAtomAPolAB=NULL;
  real **atomAreaA=NULL, **atomAreaB=NULL,**atomAreaAB=NULL;
  t_AtomProp radtype, gamma;
  int eRadType = 0;
  gmx_bool *bAtomA, *bAtomB, *bAtomAB;
  gmx_bool bPolar = FALSE, bAPolar = FALSE;
  double TmpApolarEnergy;

  //Variable for index file
  int *isize, foc_isize;	//Number of index group
  atom_id **index, *foc_index;
  char **grpnm;

  //Variable for trajectory
  int natoms; //Number of atoms in the trajectory
  real t; //Time
  rvec *x; //Coordinate
  t_trxstatus *status;

  //Variable related to output files
  const char *fnVacMM = NULL, *fnDecompMM = NULL, *fnDecompPol = NULL, *fnDecompAPol = NULL, *fnPolar = NULL, *fnAPolar = NULL;
  FILE *fVacMM = NULL, *fDecompMM = NULL, *fDecompPol = NULL,*fDecompAPol = NULL, *fPolar = NULL, *fAPolar = NULL;
  char **legVacMM = NULL, **legPolar = NULL, **legAPolar = NULL;
  char inMdp[256];

  char fnPQR[256], fnPolAPBS[256], fnApbsParamAPol[256], fnAPolAPBS[256], fnApbsOut[256];
  t_PolKey PolarKeyWords;
  t_APolKey APolarKeyWords;


#define NFILE asize(fnm)

  //Copyright message to the screen
  CopyRightMsg();

  //To show the option on the screen and to take the all option
  parse_common_args(&argc, argv,
      PCA_CAN_TIME | PCA_CAN_VIEW | PCA_TIME_UNIT | PCA_BE_NICE, NFILE, fnm,
      asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

  // change rvdw to Angstrom
  rvdw = rvdw * 10;

  /* Reverting silent to verbosity for easy use */
  if (bVerbose==FALSE)
	  bVerbose = TRUE;
  else
	  bVerbose = FALSE;

  if(!fn2bTPX(ftp2fn(efTPX, NFILE, fnm)))
	  gmx_fatal(FARGS, "tpr file is necessary\n");

  if ((!bPBSA) && (!bMM))
    gmx_fatal(FARGS, "No calculation opted. Use either \"-pbsa\" or \"-mm\" option.\n");

  if (bPBSA)
    {
      if (ftp2fn_null(efMDP, NFILE, fnm) == NULL )
        gmx_fatal(FARGS, "Input parameter file for the PBSA calculation is missing, use \"-i\" option\n");
      else
        {
          sprintf(inMdp, "%s", ftp2fn(efMDP, NFILE, fnm));
          //Reading input file and creating structure for all parameters
          ReadInput(inMdp, &bPolar, &bAPolar, &PolarKeyWords,&APolarKeyWords);
        }

      //Creating name for temporary data file
      strcpy(buf, "pbXXXXX");
      gmx_tmpnam(buf);
      sprintf(fnPQR, "%sA.pqr", buf);
      sprintf(fnPolAPBS,"%sA.in" , buf);
      sprintf(fnApbsParamAPol,"%s_param.dat", buf);
      sprintf(fnAPolAPBS,"AP%sA.in" , buf);
      sprintf(fnApbsOut,"%s.out" , buf);
      remove(buf);
    }


  read_tps_conf(ftp2fn(efTPX, NFILE, fnm), title, &top, &ePBC, &xtop, NULL, box, FALSE);
  atoms = &(top.atoms);

  if ((bMM) && (!bDIFF) && (!bIncl14))
	{
	  	  printf("\n\nWARNING: For single group calculations, 1-4 interactions are also included.\n\n");
	  	  bIncl14 = TRUE;
	}

  ////Modify residue name having same name but different structure////
  ////Important for APBS Parameter file prepared for non-polar solvation energy////
  snew(modresname, top.atoms.nr);
  for (i = 0; i < top.atoms.nr; i++)
    snew(modresname[i], 4);
  modResname(&top, modresname);
  ////END////


////////To change radius according to the given option//
  if (strcmp(rtype[0], "bondi") == 0)
  {
	  Bondi(&radtype);
	  eRadType = eBondi;
  }

  if (strcmp(rtype[0], "mbondi") == 0)
  {
    mBondi(&radtype);
    eRadType = eMbondi;
  }

  if (strcmp(rtype[0], "mbondi2") == 0)
  {
    mBondi2(&radtype);
    eRadType = eMbondi2;
  }

  if (strcmp(rtype[0], "force-field") == 0)
	  eRadType = eFF;

///////////END/////


  if (!bDIFF)
    {
//////////////TO SELECT SPECIFIC INEDEX GROUP for single group calculation////
      snew(isize, 1);
      snew(index, 1);
      snew(grpnm, 1);
      printf("\nEnter a group number for energy calculation:\n");
      get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize[0], &index[0], grpnm);
//////////END////

////////Assigning CHARGE and RADIUS for single group calculation////
      if((bPBSA) && (!bDIFF))
		assignQR(&top, index[0], isize[0], eRadType,  &radtype, rvdw);

      //snew(radiusA, top.atoms.nr);
        if(bAPolar)
      	  SasvRad(&top, isize[0], index[0], &radiusA, APolarKeyWords);
///////////END////

/////////Focusing index if enabled////////
      if((bPolar) && (bFocus)){
        	  PolarKeyWords.bFocus = TRUE;
        	  printf("\n\n\nEnter the group number to select the region for focusing \n\n");
        	  get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &foc_isize, &foc_index, grpnm);
      }
/////////Focusing index if enabled////////



/////////LEGENDS for XVG output files////
      if (bMM)
        {
          snew(legVacMM, 3);
          sprintf(buf, "VdW Energy");
          legVacMM[0] = strdup(buf);
          sprintf(buf, "Electrostatic Energy");
          legVacMM[1] = strdup(buf);
          sprintf(buf, "Total Energy");
          legVacMM[2] = strdup(buf);
        }
      if (bPolar)
        {
          snew(legPolar, 1);
          sprintf(buf, "PB Energy");
          legPolar[0] = strdup(buf);
        }
      if (bAPolar)
      {
    	  snew(legAPolar, 4);
	  	  sprintf(buf, "Surf-ten energy");
	  	  legAPolar[0] = strdup(buf);
	  	  sprintf(buf, "Press-Vol energy");
	  	  legAPolar[1] = strdup(buf);
	  	  sprintf(buf, "WCA energy");
	  	  legAPolar[2] = strdup(buf);
	  	  sprintf(buf, "Total energy");
	  	  legAPolar[3] = strdup(buf);

          if(APolarKeyWords.bWCA)	{
 ////////////////Generation of APBS Parameter file for non-polar calculation
        	  ApbsParamAPol(&top, index[0], isize[0], modresname, fnApbsParamAPol);
///////////////END//
          }
        }
    }
  else
    {

////////To select index group for difference calculation///
      snew(isize, 2);
      snew(index, 2);
      snew(grpnm, 2);
      printf(
          "\n\n\nEnter the group number for Protein or first Protein or first group:\n");
      get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize[0], &index[0],
          &grpnm[0]);
      printf(
          "\n\n\nEnter the group number of Ligand or second Protein or second group:\n");
      get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize[1], &index[1],
          &grpnm[1]);
/////////END/////

//////////Merging above selected two index group//
      srenew(index, 3);
      srenew(isize, 3);
      isize[2] = isize[0] + isize[1];
      snew(index[2], isize[0]+isize[1]);
      for (i = 0; i < isize[0]; i++)
        index[2][i] = index[0][i];
      for (i = 0; i < isize[1]; i++)
        index[2][i + isize[0]] = index[1][i];
//////////END//


/////////Assigning CHARGE and RADIUS for the Difference calculation//
      if((bPBSA) && (bDIFF))
    	  assignQR(&top, index[2], isize[2], eRadType, &radtype, rvdw);
      if(bAPolar)	{
    	  SasvRad(&top, isize[0], index[0], &radiusA, APolarKeyWords);
    	  SasvRad(&top, isize[1], index[1], &radiusB, APolarKeyWords);
    	  SasvRad(&top, isize[2], index[2], &radiusAB, APolarKeyWords);
      }
/////////END//

/////////Focusing index if enabled////////
      if ((bPolar) && (bFocus))	{
        	  PolarKeyWords.bFocus = TRUE;
        	  printf("\n\n\nEnter the group number to select the region for focusing \n\n");
        	  get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &foc_isize, &foc_index, grpnm);
      }
/////////END//


      if (bMM)
        {
          snew(legVacMM, 7);
          sprintf(buf, "%s VdW Energy", grpnm[0]);
          legVacMM[0] = strdup(buf);
          sprintf(buf, "%s Elec. Energy", grpnm[0]);
          legVacMM[1] = strdup(buf);
          sprintf(buf, "%s VdW Energy", grpnm[1]);
          legVacMM[2] = strdup(buf);
          sprintf(buf, "%s Elec. Energy", grpnm[1]);
          legVacMM[3] = strdup(buf);
          if(bIncl14)	{
        	  sprintf(buf, "%s+%s VdW Energy", grpnm[0], grpnm[1]);
        	  legVacMM[4] = strdup(buf);
        	  sprintf(buf, "%s+%s Elec. Energy", grpnm[0], grpnm[1]);
        	  legVacMM[5] = strdup(buf);
          }
          else	{
              sprintf(buf, "%s-%s VdW Energy", grpnm[0], grpnm[1]);
              legVacMM[4] = strdup(buf);
              sprintf(buf, "%s-%s Elec. Energy", grpnm[0], grpnm[1]);
              legVacMM[5] = strdup(buf);
          }
          sprintf(buf, "%s-%s Total Energy", grpnm[0], grpnm[1]);
          legVacMM[6] = strdup(buf);
        }

      if (bPolar)
        {
          snew(legPolar, 3);
          sprintf(buf, "%s PB energy", grpnm[0]);
          legPolar[0] = strdup(buf);
          sprintf(buf, "%s PB energy", grpnm[1]);
          legPolar[1] = strdup(buf);
          sprintf(buf, "%s=%s PB energy", grpnm[0], grpnm[1]);
          legPolar[2] = strdup(buf);

        }

      if(bAPolar)
      {
    	  if(APolarKeyWords.bWCA)
////////////////Generation of APBS parameters file for non-polar solvation energy calculation//
    		  ApbsParamAPol(&top, index[2], isize[2], modresname, fnApbsParamAPol);
///////////////END////

		  snew(legAPolar, 9);
		  sprintf(buf, "%s-Surf-ten energy",grpnm[0]);
		  legAPolar[0] = strdup(buf);
		  sprintf(buf, "%s-Surf-ten energy",grpnm[1]);
		  legAPolar[1] = strdup(buf);
		  sprintf(buf, "%s+%s-Surf-ten energy",grpnm[0], grpnm[1]);
		  legAPolar[2] = strdup(buf);


		  sprintf(buf, "%s-Press-Vol energy",grpnm[0]);
		  legAPolar[3] = strdup(buf);
		  sprintf(buf, "%s-Press-Vol energy",grpnm[1]);
		  legAPolar[4] = strdup(buf);
		  sprintf(buf, "%s+%s-Press-Vol energy",grpnm[0], grpnm[1]);
		  legAPolar[5] = strdup(buf);


		  sprintf(buf, "%s-WCA energy",grpnm[0]);
		  legAPolar[6] = strdup(buf);
		  sprintf(buf, "%s-WCA energy",grpnm[1]);
		  legAPolar[7] = strdup(buf);
		  sprintf(buf, "%s+%s-WCA energy",grpnm[0], grpnm[1]);
		  legAPolar[8] = strdup(buf);
        }
    }
//////END////

  //BUILDING INDEX OF EACH RESIDUE FOR DECOMPOSITION CALCULATION
  int nres = 0, prev_res, curr_res; //nres = total number of residue
  int *resnmr = NULL, *ResIsize = NULL, *ResIstart = NULL; //resnmr = residue number, ResIsize = No. of atoms in each residue; Array indices according to nres
  gmx_bool *bResA=NULL, *bResB=NULL;

  if (bDCOMP)
    {
      snew(resnmr, 1);
      snew(ResIstart,1);
      snew(ResIsize, 1);

      snew(bResA,top.atoms.nres);
      snew(bResB,top.atoms.nres);

      for(i = 0; i<top.atoms.nres;i++)
      {
    	  bResA[i] = FALSE;
    	  bResB[i] = FALSE;
      }

      resnmr[0] = atoms->resinfo[atoms->atom[index[2][0]].resind].nr;
      nres = 1;
      ResIsize[0] = 1;
      ResIstart[0] = index[2][0];
      prev_res = atoms->atom[index[2][0]].resind;
      bResA[atoms->atom[index[2][0]].resind] = TRUE;
      for(i=0;i<isize[2];i++)
      {
    	  curr_res = atoms->atom[index[2][i]].resind;
    	  if(curr_res != prev_res)
    	  {
    		  nres++;
    		  srenew(resnmr,nres);
    		  srenew(ResIsize,nres);
    		  srenew(ResIstart,nres);

    		  resnmr[nres-1] = atoms->resinfo[curr_res].nr;
    		  ResIsize[nres-1]= 1;
    		  ResIstart[nres-1] = index[2][i];

    		  if(i<isize[0])
    			  bResA[atoms->atom[index[2][i]].resind] = TRUE;
    		  else
    			  bResB[atoms->atom[index[2][i]].resind] = TRUE;

    		  prev_res = curr_res;
    	  }
      }
      if(bMM)
    	  fnDecompMM = opt2fn("-mmcon", NFILE, fnm);
      if(bPolar)
    	  fnDecompPol = opt2fn("-pcon", NFILE, fnm);
	  if(bAPolar)	{
	  snew(TempAtomAPolA,isize[0]);
	  snew(TempAtomAPolB,isize[1]);
	  snew(TempAtomAPolAB,isize[2]);
	  for(i=0;i<isize[0];i++)
		  TempAtomAPolA[i] =0;
	  for(i=0;i<isize[1];i++)
		  TempAtomAPolB[i] =0;
	  for(i=0;i<isize[2];i++)
		  TempAtomAPolAB[i] =0;
	  fnDecompAPol = opt2fn("-apcon", NFILE, fnm);
	  }
    }

  t_non_bonded paramNonBond;
  if ((bMM) && (bIncl14))
    printf("\nGenerating non-bonded pair and 1-4 pair list...\n");
  if (!bDIFF)
    {
      if (bMM)
        {
///////////GENERATION OF NON BONDED PAIRs and 1-4 PAIRs List for single calculation
          energy_pair(&paramNonBond, &top, index[0], isize[0], NULL, bDIFF);
///////////END////
          fnVacMM = opt2fn("-mm", NFILE, fnm);
          snew(tmpEE,1);
          snew(tmpVdw,1);
        }
      if (bPolar)
        fnPolar = opt2fn("-pol", NFILE, fnm);

      if (bAPolar)
        fnAPolar = opt2fn("-apol", NFILE, fnm);
    }
  else
    {
      if (bMM)
        {
////////////Generation of NON BONDED PAIRs and 1-4 PAIRs List for difference calculation
    	  if (bIncl14)
    		  energy_pair(&paramNonBond, &top, index[2], isize[2], &isize[0], bDIFF);
///////////END///
          fnVacMM = opt2fn("-mm", NFILE, fnm);
          snew(tmpEE,3);
          snew(tmpVdw,3);
        }
      if (bPolar)
        fnPolar = opt2fn("-pol", NFILE, fnm);
      if (bAPolar)
        fnAPolar = opt2fn("-apol", NFILE, fnm);
    }

////Removing PBC from first frame
  natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
  gmx_rmpbc_t gpbc = NULL;
  gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms);
////END////

////OPENING of OUTPUT XVG FILES and Writing Legends
  if (!bDIFF)
    {
      if (bMM)
        {
          fVacMM = xvgropen(fnVacMM, "Vacuum MM Energy", output_env_get_time_label(oenv), "Energy (kJ/mol)", oenv);
          xvgr_legend(fVacMM, 3, (const char**) legVacMM, oenv);
        }
      if (bPolar)
        {
          fPolar = xvgropen(fnPolar, "Polar solvation energy", output_env_get_time_label(oenv), "Energy (kJ/mol)", oenv);
          xvgr_legend(fPolar, 1, (const char**) legPolar, oenv);
        }
      if (bAPolar)
        {
          fAPolar = xvgropen(fnAPolar, "Apolar solvation energy", output_env_get_time_label(oenv), " Energy (kJ/mol) ", oenv);
          xvgr_legend(fAPolar, 4, (const char**) legAPolar, oenv);
        }

    }
  else
    {
      if (bMM)
        {
          fVacMM = xvgropen(fnVacMM, "Vacuum MM Energy", output_env_get_time_label(oenv), "Energy (kJ/mol)", oenv);
          xvgr_legend(fVacMM, 7, (const char**) legVacMM, oenv);
        }
      if (bPolar)
        {
          fPolar = xvgropen(fnPolar, "Polar solvation energy", output_env_get_time_label(oenv), "Energy (kJ/mol)", oenv);
          xvgr_legend(fPolar, 3, (const char**) legPolar, oenv);
        }
      if (bAPolar)
        {
          fAPolar = xvgropen(fnAPolar, "APolar solvation energy", output_env_get_time_label(oenv), "Energy (kJ/mol)", oenv);
          xvgr_legend(fAPolar,9, (const char**) legAPolar, oenv);

        }
    }
///END/////

  if (bDCOMP)
    {
	  if(bMM)	{
		  fDecompMM = gmx_ffopen(fnDecompMM,"w");
		  fprintf(fDecompMM,"# Time\t");
		  for(i=0;i<top.atoms.nres;i++)
		  {
			  if((bResA[i]) || (bResB[i]))	{
				  resname= *(atoms->resinfo[i].name);
				  fprintf(fDecompMM,"%s-%d\t",resname,top.atoms.resinfo[i].nr);
			  }

		  }
		  fprintf(fDecompMM,"\n");
		  srenew(tmpEE,top.atoms.nres+3);
		  srenew(tmpVdw,top.atoms.nres+3);
	  }

	  if(bPolar)	{
		  fDecompPol = gmx_ffopen(fnDecompPol,"w");
		  fprintf(fDecompPol,"# Time\t");
		  for(i=0;i<top.atoms.nres;i++)
		  {
			  if((bResA[i]) || (bResB[i]))	{
				  resname= *(atoms->resinfo[i].name);
				  fprintf(fDecompPol,"%s-%d\t",resname,top.atoms.resinfo[i].nr);
			  }
		  }
		  fprintf(fDecompPol,"\n");
	  }

	  if(bAPolar)	{
		  fDecompAPol = gmx_ffopen(fnDecompAPol,"w");
		  fprintf(fDecompAPol,"# Time\t");
		  for(i=0;i<top.atoms.nres;i++)
		  {
			  if((bResA[i]) || (bResB[i]))	{
				  resname= *(atoms->resinfo[i].name);
				  fprintf(fDecompAPol,"%s-%d\t",resname,top.atoms.resinfo[i].nr);
			  }
		  }
		  fprintf(fDecompAPol,"\n");
	  }
    }
  //Reading trajectory file
  do
    {
	  double *TempAtomAreaA=NULL, *TempAtomAreaB=NULL, *TempAtomAreaAB=NULL;
	  double *TempAtomVolA=NULL, *TempAtomVolB=NULL, *TempAtomVolAB=NULL;
	  double *TempAtomWcaA=NULL, *TempAtomWcaB=NULL, *TempAtomWcaAB=NULL;
	  double *TempAtomPolA=NULL, *TempAtomPolB=NULL, *TempAtomPolAB=NULL;
	  double *ResEnergyPol=NULL, *ResEnergyAPol=NULL;

      gmx_rmpbc(gpbc, natoms, box, x);

      if (!bDIFF) //Single calculation start
        {
          if (bMM)
            {
              Vac_MM(x, &top, paramNonBond, pdie, bDIFF, bDCOMP, tmpEE, tmpVdw);
              fprintf(fVacMM, "%15.3lf%15.3lf%15.3lf%15.3lf\n", t, tmpVdw[0],
                  tmpEE[0], tmpEE[0] + tmpVdw[0]);
              fflush(fVacMM);
            } // bMM End
          if (bPolar)
            {
              double tempPolar;

              // Setting fine grids for the focusing calculation
              if(bFocus)
            	  psize(&top, foc_index, foc_isize, x, &PolarKeyWords,FALSE,FALSE);

              psize(&top, index[0], isize[0], x, &PolarKeyWords,TRUE,bFocus);
              makePQR(&top, index[0], isize[0], ePBC, box, x, modresname, fnPQR);
              polarInAPBS(&PolarKeyWords,fnPQR,fnPolAPBS, bDCOMP);
#ifdef INT_APBS
              apbs(argc,argv,fnPolAPBS,&tempPolar,NULL,NULL);
#else
              ext_apbs(isize[0], bVerbose, fnApbsOut, fnPolAPBS, &tempPolar, NULL, NULL);
#endif
              fprintf(fPolar,"%15.3lf%15.3lf\n",t, tempPolar);
              remove(fnPQR);remove(fnPolAPBS);
              fflush(fPolar);
            } //bPolar end
          if (bAPolar)
            {
        	  if(APolarKeyWords.gamma != 0)
        	  {
        		  if (Mod_nsc_dclm_pbc(x, radiusA[eSASRAD], isize[0], ndots, FLAG_ATOM_AREA , &tmpArea, NULL, &tmpVolume, NULL, NULL, NULL, index[0], ePBC, NULL ) != 0)
        			  gmx_fatal(FARGS, "Something wrong in the nsc_dclm_pbc\n");
        		  fprintf(fAPolar, "%15.3lf%15.3lf", t, (tmpArea*APolarKeyWords.gamma*100)+APolarKeyWords.sasaconst);
        	  }
        	  else
        		  fprintf(fAPolar, "%15.3lf%15.3lf", t, 0.0000);

        	  if (APolarKeyWords.press != 0)	{
        		  if (Mod_nsc_dclm_pbc(x, radiusA[eSAVRAD], isize[0], ndots, FLAG_VOLUME, &tmpArea, NULL, &tmpVolume, NULL, NULL, NULL, index[0], ePBC, NULL ) != 0)
        			  gmx_fatal(FARGS, "Something wrong in the nsc_dclm_pbc\n");
        		  fprintf(fAPolar, "%15.3lf",(tmpVolume*APolarKeyWords.press*1000)+APolarKeyWords.savconst);
        	  }
        	  else
        		  fprintf(fAPolar, "%15.3lf", 0.000);

        	  if(APolarKeyWords.bWCA)
        	  {
                  makePQR(&top, index[0], isize[0], ePBC, box, x, modresname, fnPQR);
                  APolarInAPBS(&APolarKeyWords, fnPQR, fnAPolAPBS,fnApbsParamAPol);
#ifdef INT_APBS
        		  apbs(argc,argv,fnAPolAPBS,NULL,&TmpApolarEnergy, NULL);
#else
        		  ext_apbs(isize[0], bVerbose,fnApbsOut, fnAPolAPBS,NULL,&TmpApolarEnergy, NULL);
#endif
        		  fprintf(fAPolar, "%15.3lf\n", TmpApolarEnergy);
        	  }
        	  else
        		  fprintf(fAPolar, "%15.3lf\n", 0.0000);
        	  fflush(fAPolar);
            } //bApolar End
        } //Single calculation end

      else //Difference calculation start
        {
    	  snew(TempAtomWcaA,isize[0]); snew(TempAtomWcaB,isize[1]); snew(TempAtomWcaAB,isize[2]);
    	  snew(TempAtomPolA,isize[0]); snew(TempAtomPolB,isize[1]); snew(TempAtomPolAB,isize[2]);

          if (bMM)
            {
        	  if(bIncl14)
        		  Vac_MM(x, &top, paramNonBond, pdie, bDIFF, bDCOMP, tmpEE, tmpVdw);
        	  else
        		  Vac_MM_without_14(x, &top, index[0], isize[0], index[1], isize[1], pdie, bDCOMP, tmpEE, tmpVdw);
              fprintf(fVacMM, "%15.3lf%15.3lf%15.3lf", t, tmpVdw[0], tmpEE[0]);
              fprintf(fVacMM, "%15.3lf%15.3lf", tmpVdw[1], tmpEE[1]);
              fprintf(fVacMM, "%15.3lf%15.3lf%15.3lf\n", tmpVdw[2], tmpEE[2], (tmpVdw[2] + tmpEE[2]) - (tmpEE[0]+tmpVdw[0]+tmpEE[1]+tmpVdw[1]));
              fflush(fVacMM);
              if(bDCOMP)
              {
            	  fprintf(fDecompMM, "%15.3lf", t);
				  for(i=0; i<top.atoms.nres; i++)
					  if((bResA[i]) || (bResB[i]))
						  fprintf(fDecompMM, "%15.3lf", (tmpEE[i+3]+tmpVdw[i+3])/2);
				  fprintf(fDecompMM, "\n");
				  fflush(fDecompMM);
              }
            }// bMM End

          if (bPolar)
            {
        	  if(bDCOMP)	{
        	  snew(ResEnergyPol, top.atoms.nres);
        	  for(i=0;i<top.atoms.nres;i++)
        		  ResEnergyPol[i] = 0;
        	  }

              double tempPolar;
              //Setting cglen and cgcent according to the Complex
              psize(&top, index[2], isize[2], x, &PolarKeyWords,TRUE,bFocus);

              // Setting fine grids for the focusing calculation
              if(bFocus)
            	  psize(&top, foc_index, foc_isize, x, &PolarKeyWords,FALSE, FALSE);

              psize(&top, index[0], isize[0], x, &PolarKeyWords,FALSE,bFocus);
              makePQR(&top, index[0], isize[0], ePBC, box, x, modresname,fnPQR);
              polarInAPBS(&PolarKeyWords, fnPQR, fnPolAPBS, bDCOMP);
#ifdef INT_APBS
              apbs(argc, argv, fnPolAPBS, &tempPolar, NULL, TempAtomPolA);
#else
              ext_apbs(isize[0], bVerbose, fnApbsOut, fnPolAPBS, &tempPolar, NULL, TempAtomPolA);
#endif
              fprintf(fPolar, "%15.3lf%15.3lf\t", t, tempPolar);
              remove(fnPQR);
              remove(fnPolAPBS);

              psize(&top, index[1], isize[1], x, &PolarKeyWords,FALSE,bFocus);
              makePQR(&top, index[1], isize[1], ePBC, box, x, modresname,fnPQR);
              polarInAPBS(&PolarKeyWords, fnPQR, fnPolAPBS, bDCOMP);
#ifdef INT_APBS
              apbs(argc, argv, fnPolAPBS, &tempPolar, NULL, TempAtomPolB);
#else
              ext_apbs(isize[1], bVerbose, fnApbsOut, fnPolAPBS, &tempPolar, NULL, TempAtomPolB);
#endif
              fprintf(fPolar, "%15.3lf\t", tempPolar);
              remove(fnPQR);
              remove(fnPolAPBS);

              psize(&top, index[2], isize[2], x, &PolarKeyWords,FALSE,bFocus);
              makePQR(&top, index[2], isize[2], ePBC, box, x, modresname,fnPQR);
              polarInAPBS(&PolarKeyWords, fnPQR, fnPolAPBS, bDCOMP);
#ifdef INT_APBS
              apbs(argc, argv, fnPolAPBS, &tempPolar, NULL, TempAtomPolAB);
#else
              ext_apbs(isize[2], bVerbose, fnApbsOut, fnPolAPBS, &tempPolar, NULL, TempAtomPolAB);
#endif
              fprintf(fPolar, "%15.3lf\n", tempPolar);
              remove(fnPQR);
              remove(fnPolAPBS);

              fflush(fPolar);

            }
          if (bAPolar)
            {
        	  if(bDCOMP)	{
        		  snew(ResEnergyAPol, top.atoms.nres);
        		  	  for(i=0;i<top.atoms.nres;i++)
        		  		  ResEnergyAPol[i] = 0;
        	  }

        	  if(APolarKeyWords.gamma != 0)
        	  {
				  if (Mod_nsc_dclm_pbc(x, radiusA[eSASRAD], isize[0], ndots, FLAG_ATOM_AREA, &tmpArea, &TempAtomAreaA,  &tmpVolume, NULL,  NULL, NULL, index[0], ePBC, NULL ) != 0)
					gmx_fatal(FARGS, "Something wrong in the nsc_dclm_pbc\n");
				  fprintf(fAPolar, "%15.3lf%15.3lf", t, (tmpArea*APolarKeyWords.gamma*100) + APolarKeyWords.sasaconst);

				  if (Mod_nsc_dclm_pbc(x, radiusB[eSASRAD], isize[1], ndots, FLAG_ATOM_AREA, &tmpArea, &TempAtomAreaB,  &tmpVolume, NULL, NULL, NULL, index[1], ePBC, NULL ) != 0)
					gmx_fatal(FARGS, "Something wrong in the nsc_dclm_pbc\n");
				  fprintf(fAPolar, "%15.3lf", (tmpArea*APolarKeyWords.gamma*100) + APolarKeyWords.sasaconst);

				  if (Mod_nsc_dclm_pbc(x, radiusAB[eSASRAD], isize[2], ndots, FLAG_ATOM_AREA, &tmpArea, &TempAtomAreaAB, &tmpVolume, NULL, NULL, NULL, index[2], ePBC, NULL ) != 0)
					gmx_fatal(FARGS, "Something wrong in the nsc_dclm_pbc\n");
				  fprintf(fAPolar, "%15.3lf", (tmpArea*APolarKeyWords.gamma*100) + APolarKeyWords.sasaconst);

				  if(bDCOMP)	{
					  for(i=0;i<isize[0];i++)
						  TempAtomAPolA[i] +=  (TempAtomAreaA[i] * APolarKeyWords.gamma*100) + APolarKeyWords.sasaconst;
					  for(i=0;i<isize[1];i++)
						  TempAtomAPolB[i] +=  (TempAtomAreaB[i] * APolarKeyWords.gamma*100) + APolarKeyWords.sasaconst;
					  for(i=0;i<isize[2];i++)
						  TempAtomAPolAB[i] +=  (TempAtomAreaAB[i] * APolarKeyWords.gamma * 100) + APolarKeyWords.sasaconst;
				  }

            }
        	  else	{
        		  fprintf(fAPolar, "%15.3lf%15.3lf%15.3lf%15.3lf", t, 0.0000,0.0000,0.0000);
        	  }

        	  if(APolarKeyWords.press != 0)
        	  {
				  if (Mod_nsc_dclm_pbc(x, radiusA[eSAVRAD], isize[0], ndots, FLAG_VOLUME, &tmpArea, NULL, &tmpVolume, &TempAtomVolA,  NULL, NULL, index[0], ePBC, NULL ) != 0)
					gmx_fatal(FARGS, "Something wrong in the nsc_dclm_pbc\n");
				  fprintf(fAPolar, "%15.3lf", (tmpVolume*APolarKeyWords.press*1000)+APolarKeyWords.savconst);

				  if (Mod_nsc_dclm_pbc(x, radiusB[eSAVRAD], isize[1], ndots, FLAG_VOLUME, &tmpArea, NULL, &tmpVolume, &TempAtomVolB, NULL, NULL, index[1], ePBC, NULL ) != 0)
					gmx_fatal(FARGS, "Something wrong in the nsc_dclm_pbc\n");
				  fprintf(fAPolar, "%15.3lf", (tmpVolume*APolarKeyWords.press*1000)+APolarKeyWords.savconst);

				  if (Mod_nsc_dclm_pbc(x, radiusAB[eSAVRAD], isize[2], ndots, FLAG_VOLUME, &tmpArea, NULL, &tmpVolume, &TempAtomVolAB, NULL, NULL, index[2], ePBC, NULL ) != 0)
					gmx_fatal(FARGS, "Something wrong in the nsc_dclm_pbc\n");
				  fprintf(fAPolar, "%15.3lf", (tmpVolume*APolarKeyWords.press*1000)+APolarKeyWords.savconst);

				  if(bDCOMP)	{
					  for(i=0;i<isize[0];i++)
						  TempAtomAPolA[i] +=  (TempAtomVolA[i] * APolarKeyWords.press * 1000)+APolarKeyWords.savconst;
					  for(i=0;i<isize[1];i++)
						  TempAtomAPolB[i] +=  (TempAtomVolB[i] * APolarKeyWords.press * 1000)+APolarKeyWords.savconst;
					  for(i=0;i<isize[2];i++)
						  TempAtomAPolAB[i] +=  (TempAtomVolAB[i] * APolarKeyWords.press * 1000)+APolarKeyWords.savconst;
				  }

            }
        	  else	{
        		  fprintf(fAPolar, "%15.3lf%15.3lf%15.3lf", 0.0000,0.0000,0.0000);
        	  }

        	  if(APolarKeyWords.bWCA)
        	  {
                  makePQR(&top, index[0], isize[0], ePBC, box, x, modresname, fnPQR);
                  APolarInAPBS(&APolarKeyWords, fnPQR, fnAPolAPBS,fnApbsParamAPol);
#ifdef INT_APBS
                  apbs(argc,argv,fnAPolAPBS,NULL,&TmpApolarEnergy, TempAtomWcaA);
#else
                  ext_apbs(isize[0], bVerbose, fnApbsOut,fnAPolAPBS,NULL,&TmpApolarEnergy, TempAtomWcaA);
#endif
           		  fprintf(fAPolar, "%15.3lf",TmpApolarEnergy);
           		  remove(fnPQR);
           		  remove(fnAPolAPBS);

                  makePQR(&top, index[1], isize[1], ePBC, box, x, modresname, fnPQR);
                  APolarInAPBS(&APolarKeyWords, fnPQR, fnAPolAPBS,fnApbsParamAPol);
#ifdef INT_APBS
                  apbs(argc,argv,fnAPolAPBS,NULL,&TmpApolarEnergy, TempAtomWcaB);
#else
                  ext_apbs(isize[1], bVerbose,fnApbsOut,fnAPolAPBS,NULL,&TmpApolarEnergy, TempAtomWcaB);
#endif
           		  fprintf(fAPolar, "%15.3lf", TmpApolarEnergy);
          		  remove(fnPQR);
           		  remove(fnAPolAPBS);

                  makePQR(&top, index[2], isize[2], ePBC, box, x, modresname, fnPQR);
                  APolarInAPBS(&APolarKeyWords, fnPQR, fnAPolAPBS,fnApbsParamAPol);
#ifdef INT_APBS
                  apbs(argc,argv,fnAPolAPBS,NULL,&TmpApolarEnergy, TempAtomWcaAB);
#else
                  ext_apbs(isize[2], bVerbose,fnApbsOut,fnAPolAPBS,NULL,&TmpApolarEnergy, TempAtomWcaAB);
#endif
           		  fprintf(fAPolar, "%15.3lf\n", TmpApolarEnergy);

          		  remove(fnPQR);
           		  remove(fnAPolAPBS);

				  if(bDCOMP)	{
					  for(i=0;i<isize[0];i++)
						  TempAtomAPolA[i] +=  TempAtomWcaA[i];
					  for(i=0;i<isize[1];i++)
						  TempAtomAPolB[i] +=  TempAtomWcaB[i];
					  for(i=0;i<isize[2];i++)
						  TempAtomAPolAB[i] +=  TempAtomWcaAB[i];
				  }
        	  }
        	  else	{
        		  fprintf(fAPolar, "%15.3lf%15.3lf%15.3lf\n", 0.0000,0.0000,0.0000);
        	  }

        	  fflush(fAPolar);
          }
        }

      if (bDCOMP)
        {
          if(bPolar)	{
              fprintf(fDecompPol, "%15.3lf", t);
        	  decomp_calc_energy(top,isize[0], index[0], isize[1], index[1], isize[2], index[2],	\
          		             TempAtomPolA, TempAtomPolB,TempAtomPolAB, bResA, bResB, ResEnergyPol);
              for(i=0; i<top.atoms.nres; i++)
    			  if((bResA[i]) || (bResB[i]))
    				  fprintf(fDecompPol, "%15.3lf", ResEnergyPol[i]);
              fprintf(fDecompPol, "\n");
        	  sfree(ResEnergyPol);
        	  sfree(TempAtomPolA); sfree(TempAtomPolB); sfree(TempAtomPolAB);
        	  fflush(fDecompPol);
          }

          if(bAPolar)	{
        	  fprintf(fDecompAPol, "%15.3lf", t);
        	  decomp_calc_energy(top,isize[0], index[0], isize[1], index[1], isize[2], index[2],	\
          		             TempAtomAPolA, TempAtomAPolB,TempAtomAPolAB, bResA, bResB, ResEnergyAPol);
              for(i=0; i<top.atoms.nres; i++)
    			  if((bResA[i]) || (bResB[i]))
    				  fprintf(fDecompAPol, "%15.3lf", ResEnergyAPol[i]);
              fprintf(fDecompAPol, "\n");

              sfree(ResEnergyAPol);
        	  if(APolarKeyWords.gamma != 0)	{
        		  sfree(TempAtomAreaA); sfree(TempAtomAreaB); sfree(TempAtomAreaAB);
        	  }
        	  if(APolarKeyWords.press != 0)	{
        		  sfree(TempAtomVolA); sfree(TempAtomVolB); sfree(TempAtomVolAB);
        	  }
        	  if(APolarKeyWords.bWCA)	{
        		  sfree(TempAtomWcaA); sfree(TempAtomWcaB); sfree(TempAtomWcaAB);
        	  }
			  for(i=0;i<isize[0];i++)
				  TempAtomAPolA[i]  = 0;
			  for(i=0;i<isize[1];i++)
				  TempAtomAPolB[i] =  0;
			  for(i=0;i<isize[2];i++)
				  TempAtomAPolAB[i] = 0;
			  fflush(fDecompAPol);
          }
        }
      nframe++;
    }
  while (read_next_x(oenv, status, &t, x, box));

  cite(bPolar, bAPolar, APolarKeyWords);
  remove("io.mc");
  fprintf(stderr,"\n\nThanks for using g_mmpbsa.\n\n");
  return 0;
}


int main(int argc, char *argv[])
{

# ifndef INT_APBS
#ifdef GMX_NO_SYSTEM
	  gmx_fatal(FARGS,"No calls to system(3) supported on this platform.");
#endif
#endif

  gmx_run_cmain(argc, argv, &gmx_do_mmpbsa);

  return 0;
}
