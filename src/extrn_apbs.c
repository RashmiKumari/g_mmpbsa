/*
 * This file is part of g_mmpbsa.
 *
 * Authors: Rashmi Kumari and Rajendra Kumar
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "smalloc.h"
#include "gmx_fatal.h"
#include "ExtractData.h"

double get_totEnergy(char *line)	{
	double energy;
	char **split_data=NULL;
	int nwords, i;

	remove_leading_white_space(line);
	split_data = split_by_space(line, &nwords);
	energy = strtod(split_data[1], NULL);

	free(split_data);
	return energy;
}

double get_AtomEnergy(char *line)	{
	double energy;
	char **split_data=NULL;
	int nwords, i;

	remove_leading_white_space(line);
	split_data = split_by_space(line, &nwords);
	energy = strtod(split_data[2], NULL);

	free(split_data);
	return energy;
}

double get_WCAEnergy(char *line)	{
	double energy;
	char **split_data=NULL;
	int nwords, i;

	remove_leading_white_space(line);
	split_data = split_by_space(line, &nwords);
	energy = strtod(split_data[3], NULL);

	free(split_data);
	return energy;
}

double get_AtomWCAEnergy(char *line)	{
	double energy;
	char **split_data=NULL;
	int nwords, i;

	remove_leading_white_space(line);
	split_data = split_by_space(line, &nwords);
	energy = strtod(split_data[5], NULL);

	free(split_data);
	return energy;
}

int have_mpirun(const char *line)	{
	int i = 0;
	char *mpi_dict[] = {
			"mpirun",
			"mpiexec",
			"orterun"
	};
	bool bMPI=FALSE;

	for(i=0; i<3; i++)	{
		if (strstr(line, mpi_dict[i]) != NULL)	{
			bMPI = TRUE;
			break;
		}
	}

	if (bMPI)
		return TRUE;
	else
		return FALSE;
}

int ext_apbs(int isize, gmx_bool bVerbose, char *fnApbsOut, char *fnPolAPBS, double *PolarEnergy, double *APolarEnergy, double *AtomEnergy)	{
	char *apbs_cmd = NULL;
	const char *apbs_env = NULL;
	FILE *fApbsOut;
	char **data=NULL;
	bool bID_A=FALSE, bID_B =FALSE, bWCA_end=FALSE;
	int nlines, i;
	int mol1_start = 0, mol2_start = 0;
	int mol1_lastID = 0, mol2_lastID = 0;

	double totEnergy1, totEnergy2;
	double *atEnergy1=NULL, *atEnergy2=NULL;
	int at_count = 0;

	/* Getting path from $APBS environment */
	apbs_env = getenv("APBS");
	snew(apbs_cmd, strlen(apbs_env)+64);

	/* Constructing APBS command for polar solvation energy*/
	if (PolarEnergy != NULL)
		sprintf(apbs_cmd, "$APBS %s --output-file=%s %s", fnPolAPBS, fnApbsOut, bVerbose?"":" >/dev/null 2>&1");

	/* Constructing APBS command for WCA energy*/
	if (APolarEnergy != NULL)	{
		if (have_mpirun(apbs_env))
			gmx_fatal(FARGS,"Do not use MPI for WCA...");
		sprintf(apbs_cmd, "$APBS %s > %s", fnPolAPBS, fnApbsOut);
	}

	/* Executing APBS command */
	if(0 != system(apbs_cmd))
		  gmx_fatal(FARGS,"Failed to execute command: %s", apbs_cmd);

	/* Reading APBS output file */
	fApbsOut = fopen(fnApbsOut, "r");
	data = get_all_lines(fApbsOut, &nlines);
	fclose(fApbsOut);

	/* Memory allocation for atom energy */
	if (AtomEnergy != NULL)	{
		snew(atEnergy1, 1);
		snew(atEnergy2, 1);
		at_count = 1;
	}

	if (PolarEnergy != NULL)	{

		/* Identifying indices (line number) of mol1 and mol2 */
		for(i=0; i<nlines; i++)	{
			if (strstr(data[i],"elec name mol1") !=NULL)	{
				mol1_start = i;
			}
			if (strstr(data[i],"elec name mol2") !=NULL)	{
				mol2_start = i;
			}
		}

		/* Identifying last energy calculation index */
		for(i=0; i<nlines; i++)	{
			if ( (strstr(data[i],"id") !=NULL) && (i < mol2_start) )	{
				mol1_lastID = i;
			}
			if (strstr(data[i],"id") !=NULL)	{
				mol2_lastID = i;
			}
		}

		/* Extracting energy values */
		for(i=0; i<nlines; i++)	{

			if ((strstr(data[i],"id")!=NULL) && (i == mol1_lastID) )
				bID_A = TRUE;

			if( (strstr(data[i],"id")!=NULL) && (i == mol2_lastID) )
				bID_B = TRUE;

			if ( (strstr(data[i],"end")!=NULL) && (bID_A == TRUE) )	{
				bID_A =FALSE;
				if ((AtomEnergy != NULL) && (at_count-1 != isize) )	{
					gmx_fatal(FARGS,"Number of atoms in selected index (%d) does not match with number of atoms (%d) in APBS output. \n", isize, at_count);
				}
				at_count = 1;
			}

			if ( (strstr(data[i],"end")!=NULL) && (bID_B == TRUE) )	{
				bID_B =FALSE;
				if ((AtomEnergy != NULL) && (at_count-1 != isize) )	{
					gmx_fatal(FARGS,"Number of atoms in selected index (%d) does not match with number of atoms (%d) in APBS output. \n", isize, at_count);
				}
				at_count = 1;
				break;
			}

			if ( (strstr(data[i], "totEnergy") != NULL) && (bID_A == TRUE) )	{
				totEnergy1 = get_totEnergy(data[i]);
			}

			if ( (strstr(data[i], "totEnergy") != NULL) && (bID_B == TRUE) )	{
				totEnergy2 = get_totEnergy(data[i]);
			}

			if ((AtomEnergy != NULL) && (strstr(data[i], "atom") != NULL) && (bID_A == TRUE)) {
				if (at_count > 1)
					srenew(atEnergy1, at_count);
				atEnergy1[at_count-1] = get_AtomEnergy(data[i]);
				at_count += 1;
			}

			if ((AtomEnergy != NULL) && (strstr(data[i], "atom") != NULL) && (bID_B == TRUE)) {
				if (at_count > 1)
					srenew(atEnergy2, at_count);
				atEnergy2[at_count-1] = get_AtomEnergy(data[i]);
				at_count += 1;
			}
		}

		*PolarEnergy = totEnergy1 - totEnergy2;

		if (AtomEnergy != NULL)
			for(i=0; i<isize; i++)	{
				AtomEnergy[i] = atEnergy1[i] - atEnergy2[i];
			}
	}

	if (APolarEnergy != NULL)	{

		for(i=0; i<nlines; i++)	{

			if (strstr(data[i], "Total WCA energy") != NULL)	{
				*APolarEnergy = get_WCAEnergy(data[i]);
				bWCA_end = TRUE;
			}

			if (bWCA_end)	{
				if ((AtomEnergy != NULL) && (at_count-1 != isize) )	{
					gmx_fatal(FARGS,"Number of atoms in selected index (%d) does not match with number of atoms (%d) in APBS output. \n", isize, at_count);
				}
				at_count = 1;
				break;
			}

			if (strstr(data[i], "WCA energy for atom") != NULL)	{
				if (at_count > 1)
					srenew(atEnergy1, at_count);
				atEnergy1[at_count-1] = get_AtomWCAEnergy(data[i]);
				at_count += 1;
			}
		}

		if (AtomEnergy != NULL)
			for(i=0; i<isize; i++)	{
				AtomEnergy[i] = atEnergy1[i];
			}
	}

	/* memory cleanup */
	for(i=0;i<nlines;i++)
		free(data[i]);
	free(data);
	sfree(apbs_cmd);
	if (AtomEnergy != NULL)	{
		sfree(atEnergy1);
		sfree(atEnergy2);
	}

	return 0;
}
