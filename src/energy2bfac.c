/*
 * This file is part of g_mmpbsa.
 *
 * Authors: Rashmi Kumari and Andrew Lynn
 * Contribution: Rajendra Kumar
 *
 * Copyright (C) 2013 Rashmi Kumari and Andrew Lynn
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
#include "copyrite.h"
#include "vec.h"
#include "tpxio.h"

void CopyRightMsg()	{
	char *copyright[] =	{
			"                                                                        ",
			"                           :-)  g_mmpbsa (-:                            ",
			"                                                                        ",
			"               Authors: Rashmi Kumari and Andrew Lynn                   ",
			"               Contribution: Rajendra Kumar                             ",
			"                                                                        ",
			"         Copyright (C) 2013  Rashmi Kumari and Andrew Lynn              ",
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

int get_energy( char *fnEnergy, int nres, real *energy)	{
	int i =0, count=-1, dum =0;
	FILE *fin;
	real *tmpData;

	fin = ffopen(fnEnergy,"r");
	snew(tmpData, 1);

	while (fgetc(fin) != EOF)	{
		srenew(tmpData,count+2);
		fscanf(fin,"%d", &dum);
		fscanf(fin,"%g", &tmpData[count+1]);
		count++;
	}

	if(nres==count)
		printf("\nTotal number of residues matched in the energy and chosen index groups.\n");
	else
		gmx_fatal(FARGS, "Mismatched total number of Residues between energy and chosen index groups.\n");

	for (i=0; i<nres; i++)
		energy[i] = tmpData[i];

	sfree(tmpData);
	return 0;
}


int	gmx_energy2bfac (int argc, char *argv[])		{
  const char *desc[] =
    { "This tool maps the binding energy contribution of each residue on the structure.",
      "The energy will be written in the B-factor field of output PDB file/s",
      "This PDB file can be used with any molecular visualizer and ",
      "residues can be viewed in color according to the energy of respective residues.",
      "The molecular visualizer should support method to color residues by the B-factor values."
	};


  t_filenm fnm[] =
  {
    { efTPS, "-s", NULL, ffREAD },
    { efDAT, "-i", "decomp_energy", ffREAD },
    { efNDX, "-n", "index", ffOPTRD },
    { efPDB, "-c", "complex.pdb", ffOPTWR },
    { efPDB, "-s1", "subunit_1.pdb", ffOPTWR },
    { efPDB, "-s2", "subunit_2.pdb", ffOPTWR },
  };

  char title[STRLEN], buf[256];
  t_topology top;
  int ePBC;
  rvec *xtop;
  matrix box;
  t_atoms *atoms;
  char *atomtype, *resname, *atomname, **modresname = NULL;
  output_env_t oenv;

  //Variable for index file
  int *isize;	//Number of index group
  atom_id **index;
  char **grpnm;

  int i=0,j=0,k=0;
  gmx_bool *bAtomA, *bAtomB, *bAtomAB;
  int nres = 0, prev_res, curr_res;
  gmx_bool *bResA, *bResB;
  real *energy;
  char fnEnergy[256];
  FILE *fComplex, *fS1, *fS2;

  #define NFILE asize(fnm)
  CopyRightMsg();
  //To show the option on the screen and to take the all option
  parse_common_args(&argc, argv,
      PCA_CAN_TIME | PCA_CAN_VIEW | PCA_TIME_UNIT | PCA_BE_NICE, NFILE, fnm,
      0, NULL, asize(desc), desc, 0, NULL, &oenv);

  read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &xtop, NULL, box, FALSE);
  atoms = &(top.atoms);

  if (opt2fn_null("-i", NFILE, fnm) == NULL )
    gmx_fatal(FARGS, "Residue energy file is missing, use \"-i\" option\n");
  else
      sprintf(fnEnergy, "%s", opt2fn("-i", NFILE, fnm));

  snew(isize, 2);
  snew(index, 2);
  snew(grpnm, 2);
  printf(
      "\n\n\nEnter the group number for Protein or first Protein or first group\n\n");
  get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize[0], &index[0],
      &grpnm[0]);
  printf(
      "\n\n\nEnter the group number of Ligand or second Protein or second group\n\n");
  get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &isize[1], &index[1],
      &grpnm[1]);

  srenew(index, 3);
  srenew(isize, 3);
  isize[2] = isize[0] + isize[1];
  snew(index[2], isize[0]+isize[1]);
  for (i = 0; i < isize[0]; i++)
	  index[2][i] = index[0][i];
  for (i = 0; i < isize[1]; i++)
	  index[2][i + isize[0]] = index[1][i];

  snew(bResA,top.atoms.nres);
  snew(bResB,top.atoms.nres);

  for(i = 0; i<top.atoms.nres;i++)
  {
	  bResA[i] = FALSE;
	  bResB[i] = FALSE;
  }

  nres = 1;
  prev_res = atoms->atom[index[2][0]].resind;
  bResA[atoms->atom[index[2][0]].resind] = TRUE;

  for(i=0;i<isize[2];i++)
  {
	  curr_res = atoms->atom[index[2][i]].resind;
	  if(curr_res != prev_res)
	  {
		  nres++;
		  if(i<isize[0])
			  bResA[atoms->atom[index[2][i]].resind] = TRUE;
		  else
			  bResB[atoms->atom[index[2][i]].resind] = TRUE;

		  prev_res = curr_res;
	  }
  }

  printf("\nTotal number of residues in topology: %d \n", nres);

  snew(energy,nres);
  get_energy( fnEnergy, nres, energy);

  snew(atoms->pdbinfo, atoms->nr);

  for(i = 0; i<isize[2];i++)
  {
	  if((bResA[atoms->atom[index[2][i]].resind]) || (bResB[atoms->atom[index[2][i]].resind]))	{
		 atoms->pdbinfo[index[2][i]].bfac = energy[atoms->atom[index[2][i]].resind];
		 atoms->pdbinfo[index[2][i]].occup = 1.00;
	  }
	  else{
		  atoms->pdbinfo[index[2][i]].bfac = 0.00;
		  atoms->pdbinfo[index[2][i]].occup = 1.00;
	  }
  }

  fComplex = ffopen(opt2fn("-c", NFILE, fnm),"w");
  write_pdbfile_indexed(fComplex,NULL,atoms,xtop,ePBC,box,' ',-1,isize[2],index[2],NULL,TRUE);

  fS1 = ffopen(opt2fn("-s1", NFILE, fnm),"w");
  write_pdbfile_indexed(fS1,NULL,atoms,xtop,ePBC,box,' ',-1,isize[0],index[0],NULL,TRUE);

  fS2 = ffopen(opt2fn("-s2", NFILE, fnm),"w");
  write_pdbfile_indexed(fS2,NULL,atoms,xtop,ePBC,box,' ',-1,isize[1],index[1],NULL,TRUE);

return 0;
}


int main(int argc, char *argv[])
{
  gmx_energy2bfac(argc, argv);
  return 0;
}

