/**
 * This file is part of g_mmpbsa.
 *
 * Authors: Rajendra Kumar, Rashmi Kumari and Andrew Lynn
 *
 * Copyright (C) 2013-2021 Rashmi Kumari and Andrew Lynn
 * Copyright (C) 2022- Rajendra Kumar and Rashmi Kumari
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
#include <string.h>
#include <math.h>

#include "gromacs/utility/fatalerror.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/cmdlineinit.h"

/**
 * Copyright message to dispaly at the start of execution
 */
void CopyRightMsg()
{
    std::string copyright =
        "                                                                        \n"
        "                        :-)  g_mmpbsa (-:                               \n"
        "                                                                        \n"
        "                                                                        \n"
        "       Copyright (C) 2013 - 2021 Rashmi Kumari and Andrew Lynn          \n"
        "       Copyright (C) 2022- Rajendra Kumar and Rashmi Kumari             \n"
        "                                                                        \n"
        "g_mmpbsa is free software: you can redistribute it and/or modify        \n"
        "it under the terms of the GNU General Public License as published by    \n"
        "the Free Software Foundation, either version 3 of the License, or       \n"
        "(at your option) any later version.                                     \n"
        "                                                                        \n"
        "g_mmpbsa is distributed in the hope that it will be useful,             \n"
        "but WITHOUT ANY WARRANTY; without even the implied warranty of          \n"
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           \n"
        "GNU General Public License for more details.                            \n"
        "                                                                        \n"
        "You should have received a copy of the GNU General Public License       \n"
        "along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.       \n"
        "                                                                        \n"
        "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     \n"
        "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     \n"
        "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   \n"
        "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    \n"
        "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   \n"
        "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED\n"
        "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  \n"
        "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  \n"
        "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    \n"
        "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      \n"
        "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            \n"
        "                                                                        \n"
        "                                                                        \n"
        "                                                                        \n";
    fprintf ( stderr,"%s\n", copyright.c_str() );
}

int get_energy( char *fnEnergy, int nres, real *energy)	{
	int i =0, count=-1, dum =0;
    double tempDouble;
	FILE *fin;
	real *tmpData;
    const char *informat = sizeof(tmpData) == sizeof(tempDouble) ? "%lf" : "%f";

	fin = gmx_ffopen(fnEnergy,"r");
	snew(tmpData, 1);

	while (fgetc(fin) != EOF)	{
		srenew(tmpData,count+2);
		fscanf(fin,"%d", &dum);
		fscanf(fin,informat, &tmpData[count+1]);
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
    { "It maps the binding energy contribution of each residue on the structure.",
      "The energy will be written in the B-factor field of the output PDB file/s",
      "This PDB file can be used with any molecular visualizer and ",
      "residues can be viewed in color according to the energy of respective residue.",
      "The molecular visualizer should support method to color residues by B-factor values."
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

  char title[256], buf[256];
  t_topology top;
  PbcType ePBC;
  rvec *xtop;
  matrix box;
  t_atoms *atoms;
  gmx_output_env_t *oenv;

  //Variable for index file
  int *isize;	//Number of index group
  int **index;
  char **grpnm;

  int i=0,j=0,k=0;
  int nres = 0, prev_res, curr_res;
  gmx_bool *bResA, *bResB;
  real *energy;
  char fnEnergy[256];
  FILE *fComplex, *fS1, *fS2;

  #define NFILE asize(fnm)
  CopyRightMsg();
  //To show the option on the screen and to take the all option
  parse_common_args(&argc, argv,
      PCA_CAN_TIME | PCA_CAN_VIEW | PCA_TIME_UNIT , NFILE, fnm,
      0, NULL, asize(desc), desc, 0, NULL, &oenv);

  read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &ePBC, &xtop, NULL, box, FALSE);
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

  fComplex = gmx_ffopen(opt2fn("-c", NFILE, fnm),"w");
  write_pdbfile_indexed(fComplex,NULL,atoms,xtop,ePBC,box,' ',-1,isize[2],index[2],NULL,TRUE);

  fS1 = gmx_ffopen(opt2fn("-s1", NFILE, fnm),"w");
  write_pdbfile_indexed(fS1,NULL,atoms,xtop,ePBC,box,' ',-1,isize[0],index[0],NULL,TRUE);

  fS2 = gmx_ffopen(opt2fn("-s2", NFILE, fnm),"w");
  write_pdbfile_indexed(fS2,NULL,atoms,xtop,ePBC,box,' ',-1,isize[1],index[1],NULL,TRUE);

return 0;
}


int main(int argc, char *argv[])
{
  gmx_run_cmain(argc, argv, &gmx_energy2bfac);
  return 0;
}
