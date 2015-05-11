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
#include <stdio.h>
#include <math.h>
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "copyrite.h"
#include "vec.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "xvgr.h"

#include "g_mmpbsa.h"


int energy_pair(t_non_bonded *param, t_topology *top, atom_id *index, int isize, int *splitIndex, gmx_bool bDiff){

	int i, j=0, k,l,n;
	int itypeA, itypeB, ntype = top->idef.atnr;
	gmx_bool bExclude=FALSE, bLJ14=FALSE;

	param->nr_nb = 1;
	param->nr_14 = 1;
	param->pairNB = malloc(sizeof(int*));
	param->pair14 = malloc(sizeof(int*));
	param->pairtype = malloc(sizeof(int));
	if(bDiff)
	{
		snew(param->bItsA,1);
		snew(param->bItsB,1);
		snew(param->bItsA14,1);
		snew(param->bItsB14,1);
	}

	float progress=0;
	for(i=0;i<isize;i++)
	{
		for(j=0;j<i;j++)	{
			bExclude=FALSE;
			bLJ14=FALSE;
			//To loop over 4th atoms for every pair, time consuming
			/*
			for(k=0;k<top->excls.nr;k++)
			{
				n = top->excls.index[k+1]-top->excls.index[k];
				for(l=0;l<n;l++){
					if(((index[i]== k) && (index[j] == top->excls.a[top->excls.index[k]+l])) || \
							((index[j]== k) && (index[i] == top->excls.a[top->excls.index[k]+l])))
					{
						//printf("===%d=%d=%d=%d==\n",i,j,k,top->excls.a[top->excls.index[k]+l]);
						bExclude=TRUE;
						break;
					}
				}
			}*/

			n = top->excls.index[index[i]+1] - top->excls.index[index[i]];
			for(l=0;l<n;l++){
				if(index[j] == top->excls.a[top->excls.index[index[i]]+l])
				{
					bExclude=TRUE;
					break;
				}
			}
			k = 0;
			while(k<top->idef.il[F_LJ14].nr)
			{
				if(((index[i] == top->idef.il[F_LJ14].iatoms[k+1]) && (index[j] == top->idef.il[F_LJ14].iatoms[k+2])) || \
											((index[j]== top->idef.il[F_LJ14].iatoms[k+1]) && (index[i] == top->idef.il[F_LJ14].iatoms[k+2])))
				{
					param->pairtype = realloc(param->pairtype,(param->nr_14*sizeof(int)));
					param->pairtype[param->nr_14-1] = k;

					param->pair14 = realloc(param->pair14, (param->nr_14*sizeof(int*)));
					param->pair14[param->nr_14-1] = malloc(2*sizeof(int));
					param->pair14[param->nr_14-1][0] = index[i];
					param->pair14[param->nr_14-1][1] = index[j];
					bLJ14=TRUE;

					if(bDiff)	{
						srenew(param->bItsA14,param->nr_14);
						srenew(param->bItsB14,param->nr_14);
						param->bItsA14[param->nr_14-1] = FALSE;
						param->bItsB14[param->nr_14-1] = FALSE;
						if(( i < *splitIndex) && (  j  < *splitIndex ))
						{
							param->bItsA14[param->nr_14-1] = TRUE;
							param->bItsB14[param->nr_14-1] = FALSE;
						}

						if(( i >= *splitIndex) && ( j >= *splitIndex))
						{
							param->bItsB14[param->nr_14-1] = TRUE;
							param->bItsA14[param->nr_14-1] = FALSE;
						}
					}
					param->nr_14++;
					break;
				}
				k= k+3;
			}
			k=0;

			if(!bExclude)	{
				param->pairNB = realloc(param->pairNB, (param->nr_nb*sizeof(int*)));
				param->pairNB[param->nr_nb-1] = malloc(2*sizeof(int));
				param->pairNB[param->nr_nb-1][0] = index[i];
				param->pairNB[param->nr_nb-1][1] = index[j];

				if(bDiff)	{
					srenew(param->bItsA,param->nr_nb);
					srenew(param->bItsB,param->nr_nb);
					param->bItsA[param->nr_nb-1] = FALSE;
					param->bItsB[param->nr_nb-1] = FALSE;
					if(( i <  *splitIndex) && ( j  <  *splitIndex))
					{
						param->bItsA[param->nr_nb-1] = TRUE;
						param->bItsB[param->nr_nb-1] = FALSE;
					}
					if(( i  >= *splitIndex) && ( j  >= *splitIndex))
					{
						param->bItsB[param->nr_nb-1] = TRUE;
						param->bItsA[param->nr_nb-1] = FALSE;
					}
				}
				param->nr_nb++;
			}
		}
		progress = ((float)i/(float)isize) * 100;
		fprintf(stderr,"\r %5.0f %% completed...",progress);
		fflush(stdout);
	}
	printf("\n Finished pair generation....\nTotal %d 1-4 pairs and %d non-bonded pairs generated.\n\n",param->nr_14-1,param->nr_nb-1);
	param->nr_14 = param->nr_14-1;
	param->nr_nb = param->nr_nb-1;
	return 0;
}


void Vac_MM(rvec *x, t_topology *top, t_non_bonded param, real pdie, gmx_bool bDiff, gmx_bool bDCOMP, double *EE, double *Vdw)
{
	rvec dx;
	real colmb_factor = 138.935485;
	double qi, qj, c6, c12, c6j, c12j,rij, c6ij, c12ij;
	int itypeA, itypeB, ntype = top->idef.atnr;
	int i, j,k,l,n;
	int atomA, atomB, resA, resB;
	real TempEE, TempVdw;
	int nres = top->atoms.nres;

	real *EERes, *VdwRes;

	if(bDiff)
	{
		EE[0] = 0; Vdw[0] = 0;
		EE[1] = 0; Vdw[1] = 0;
		EE[2] = 0; Vdw[2] = 0;
	}
	else
	{
		EE[0] = 0; Vdw[0] = 0;
	}

	snew(EERes,nres);
	snew(VdwRes,nres);
	if(bDCOMP)
		for (i=0; i<nres;i++)
		{
			EE[i+3] = 0;
			Vdw[i+3] = 0;
			EERes[i] = 0;
			VdwRes[i] = 0;
		}

	for(i=0;i<param.nr_nb;i++)	{
		atomA = param.pairNB[i][0];
		atomB = param.pairNB[i][1];
		resA = top->atoms.atom[atomA].resind;
		resB = top->atoms.atom[atomB].resind;

		rij = sqrt(pow((x[atomA][0]-x[atomB][0]),2)+pow((x[atomA][1]-x[atomB][1]),2) + pow((x[atomA][2]-x[atomB][2]),2));

		itypeA = top->atoms.atom[atomA].type;
		itypeB = top->atoms.atom[atomB].type;

		if(itypeA<=itypeB)	{
			c6 = top->idef.iparams[itypeA*ntype+itypeB].lj.c6;
			c12 = top->idef.iparams[itypeA*ntype+itypeB].lj.c12;
		}
		else	{
			c6 = top->idef.iparams[itypeB*ntype+itypeA].lj.c6;
			c12 = top->idef.iparams[itypeB*ntype+itypeA].lj.c12;
		}
		TempEE = ((colmb_factor/pdie) * (top->atoms.atom[atomA].q*top->atoms.atom[atomB].q))/rij;
		TempVdw = (c12/pow(rij,12)) - (c6/pow(rij,6));

		if(bDiff)	{
			if(param.bItsA[i])
			{
				EE[0] += TempEE;
				Vdw[0] += TempVdw;
			}
			if(param.bItsB[i])
			{
				EE[1] += TempEE;
				Vdw[1] += TempVdw;
			}
			EE[2] += TempEE;
			Vdw[2] += TempVdw;
		}
		else
		{
			EE[0] += TempEE;
			Vdw[0] += TempVdw;
		}

		if(bDCOMP)
			if((!param.bItsA[i]) && (!param.bItsB[i]))
			{
				EERes[resA] += TempEE;
				EERes[resB] += TempEE;
				VdwRes[resA] += TempVdw;
				VdwRes[resB] += TempVdw;
			}
	}
	for(i=0;i<param.nr_14;i++)
	{
		atomA = param.pair14[i][0];
		atomB = param.pair14[i][1];
		resA = top->atoms.atom[atomA].resind;
		resB = top->atoms.atom[atomB].resind;

		rij = sqrt(pow((x[atomA][0]-x[atomB][0]),2)+pow((x[atomA][1]-x[atomB][1]),2) + pow((x[atomA][2]-x[atomB][2]),2));

		c6 = top->idef.iparams[top->idef.il[F_LJ14].iatoms[param.pairtype[i]]].lj14.c6A;
		c12 = top->idef.iparams[top->idef.il[F_LJ14].iatoms[param.pairtype[i]]].lj14.c12A;

		TempEE = top->idef.fudgeQQ * ((colmb_factor/pdie) * (top->atoms.atom[atomA].q*top->atoms.atom[atomB].q))/rij;
		TempVdw = (c12/pow(rij,12)) - (c6/pow(rij,6));

		if(bDiff)
		{
			if(param.bItsA14[i])
			{
				EE[0] += TempEE;
				Vdw[0] += TempVdw;
			}
			if(param.bItsB14[i])
			{
				EE[1] += TempEE;
				Vdw[1] += TempVdw;
			}
			EE[2] += TempEE;
			Vdw[2] += TempVdw;
		}
		else
		{
			EE[0] += TempEE;
			Vdw[0] += TempVdw;
		}
		if (bDCOMP)
			if((!param.bItsA14[i]) && (!param.bItsB14[i]))
			{
				EERes[resA] += TempEE;
				EERes[resB] += TempEE;
				VdwRes[resA] += TempVdw;
				VdwRes[resB] += TempVdw;
			}
	}

	if(bDCOMP)
		for (i=0; i<nres; i++)
		{
			EE[i+3] = EERes[i];
			Vdw[i+3] = VdwRes[i];
		}
}

void Vac_MM_without_14(rvec *x, t_topology *top, atom_id *indexA, int isizeA, atom_id *indexB, int isizeB, real pdie, gmx_bool bDCOMP, double *EE, double *Vdw){

	int i, j=0, k,l,n;
	rvec dx;
	real colmb_factor = 138.935485;
	double qi, qj, c6, c12, c6j, c12j,rij, c6ij, c12ij;
	int itypeA, itypeB, ntype = top->idef.atnr;
	int resA, resB;
	real TempEE, TempVdw;
	int nres = top->atoms.nres;

	real *EERes, *VdwRes;

	EE[0] = 0; Vdw[0] = 0;
	EE[1] = 0; Vdw[1] = 0;
	EE[2] = 0; Vdw[2] = 0;


	snew(EERes,nres);
	snew(VdwRes,nres);
	if(bDCOMP)
		for (i=0; i<nres;i++)
		{
			EE[i+3] = 0;
			Vdw[i+3] = 0;
			EERes[i] = 0;
			VdwRes[i] = 0;
		}

	for(i=0;i<isizeA;i++)
	{
		for(j=0;j<isizeB;j++)	{

			resA = top->atoms.atom[indexA[i]].resind;
			resB = top->atoms.atom[indexB[j]].resind;

			rij = sqrt(pow((x[indexA[i]][0]-x[indexB[j]][0]),2)+pow((x[indexA[i]][1]-x[indexB[j]][1]),2) + pow((x[indexA[i]][2]-x[indexB[j]][2]),2));

			itypeA = top->atoms.atom[indexA[i]].type;
			itypeB = top->atoms.atom[indexB[j]].type;

			if(itypeA<=itypeB)	{
				c6 = top->idef.iparams[itypeA*ntype+itypeB].lj.c6;
				c12 = top->idef.iparams[itypeA*ntype+itypeB].lj.c12;
			}
			else	{
				c6 = top->idef.iparams[itypeB*ntype+itypeA].lj.c6;
				c12 = top->idef.iparams[itypeB*ntype+itypeA].lj.c12;
			}
			TempEE = ((colmb_factor/pdie) * (top->atoms.atom[indexA[i]].q*top->atoms.atom[indexB[j]].q))/rij;
			TempVdw = (c12/pow(rij,12)) - (c6/pow(rij,6));

			EE[2] += TempEE;
			Vdw[2] += TempVdw;

			if(bDCOMP)
			{
				EERes[resA] += TempEE;
				EERes[resB] += TempEE;
				VdwRes[resA] += TempVdw;
				VdwRes[resB] += TempVdw;
			}
		}
	}

	if(bDCOMP)
		for (i=0; i<nres; i++)
		{
			EE[i+3] = EERes[i];
			Vdw[i+3] = VdwRes[i];
		}
}
