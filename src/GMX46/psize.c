/**
 * @file GMX46/psize.c
 * @brief Definition of routines for grid dimension calculation
 * @ingroup PSIZE
 * @author Rajendra Kumar
 * @attention
 * @verbatim
 *
 * Originally written by Dave Sept.
 * Additional APBS-specific features added by Nathan Baker.
 * Ported to Python/Psize class by Todd Dolinsky and subsequently hacked by Nathan Baker.
 *
 * ------------------------------------------------------
 * Ported to C by Rajendra Kumar for g_mmpbsa
 * ------------------------------------------------------
 * @endverbatim
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gromacs/statutil.h"
#include "gromacs/typedefs.h"
#include "gromacs/smalloc.h"
#include "gromacs/vec.h"
#include "gromacs/tpxio.h"
#include "gromacs/rmpbc.h"
#include "gromacs/xvgr.h"

#include "g_mmpbsa.h"

int maxIndex (ivec n){
	int i, max=0, imax=0;
	for(i=0;i<DIM;i++)	{
		if(n[i]>max)	{
			max = n[i];
			imax = i;
		}
	}
	return imax;
}

int psize (t_topology *top, atom_id *index, int isize, rvec *x, t_PolKey *param, gmx_bool bCG, gmx_bool bFocus)	{
	int i,j;
	rvec minlen, maxlen;
	rvec olen, clen, flen, cen;
	ivec n, np, tn, nsmall;
	real nsmem, gmem, zofac;
	real r, np_float;

	minlen[XX] = 9999; minlen[YY] = 9999; minlen[ZZ] = 9999;
	maxlen[XX] = -9999; maxlen[YY] = -9999; maxlen[ZZ] = -9999;

	for (i=0;i<isize;i++)	{
		r = top->atoms.pdbinfo[index[i]].bfac;
		for(j=0;j<DIM;j++)	{
			if((x[index[i]][j]*10)-r < minlen[j])
				minlen[j] = (x[index[i]][j]*10) - r;

			if((x[index[i]][j]*10)+r > maxlen[j])
				maxlen[j] = (x[index[i]][j]*10) + r;

		}
	}

	for (i=0;i<DIM;i++)		{
		olen[i] = maxlen[i] - minlen[i];
		clen[i] = param->cfac * olen[i];
		flen[i] = param->fadd + olen[i];
		if(flen[i]>clen[i])
			flen[i] = clen[i];
		cen[i] = (maxlen[i] + minlen[i])/2;

		tn[i] = (int) flen[i]/param->gridspace + 0.5;
		n[i] = 32*((int)((tn[i] - 1) / 32.0 + 0.5)) + 1;
		nsmall[i] = 32*((int)((tn[i] - 1) / 32.0 + 0.5)) + 1;

		if (nsmall[i] < 33)
			nsmall[i] = 33;
	}

	//To Check the available memory
	gmem = 200.0 * n[XX] * n[YY] * n[ZZ] / 1024 / 1024;
	while(1)	{
		nsmem = 200.0 * nsmall[XX] * nsmall[YY] * nsmall[ZZ] / 1024 / 1024;
		if(nsmem<param->gmemceil)
			break;
		else	{
			i = maxIndex(nsmall);
			nsmall[i] = 32 * ((nsmall[i] - 1)/32 - 1) + 1;
			if (nsmall <= 0)	{
				gmx_fatal(FARGS, "You picked a memory ceiling that is too small\n");
			}


		}
	}

	// Calculating pdime => np
	if (gmem >= param->gmemceil)	{
		zofac = 1 + 2 * param->ofrac;
		for (i=0;i<DIM;i++)	{
			np_float = n[i]/(float)nsmall[i];
			if (np_float > 1)
				np[i] = (int)(zofac*n[1]/nsmall[i] + 1.0);
		}
	}


	if(gmem >= param->gmemceil)
		param->mg_type = mg_para;
	else
		param->mg_type = mg_auto;

	if (bCG)	{
		copy_rvec(clen,param->cglen);
		copy_rvec(cen,param->cgcent);
	}

	if(!bFocus)	{
		copy_rvec(flen,param->fglen);
		copy_rvec(cen,param->fgcent);
	}
	copy_ivec(nsmall,param->dime);

	if(param->mg_type == mg_para)
		copy_ivec(np, param->pdime);

	return 0;
}
