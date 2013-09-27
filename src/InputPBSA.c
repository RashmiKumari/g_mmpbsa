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


#include "typedefs.h"
#include "readinp.h"
#include "warninp.h"
#include "string2.h"
#include "g_mmpbsa.h"

#define MAXPTR 254
static char grid[STRLEN], polar[STRLEN], apolar[STRLEN], wca[STRLEN];


//Copied from the GROMACS/src/kernel,readir.h"
static int str_nelem(const char *str,int maxptr,char *ptr[])
{
  int  np=0;
  char *copy0,*copy;

  copy0=strdup(str);
  copy=copy0;
  ltrim(copy);
  while (*copy != '\0') {
    if (np >= maxptr)
      gmx_fatal(FARGS,"Too many groups on line: '%s' (max is %d)", str,maxptr);
    if (ptr)
      ptr[np]=copy;
    np++;
    while ((*copy != '\0') && !isspace(*copy))
      copy++;
    if (*copy != '\0') {
      *copy='\0';
      copy++;
    }
    ltrim(copy);
  }
  if (ptr == NULL)
    sfree(copy0);

  return np;
}



int ReadInput	(char *fnMdp, gmx_bool *bPolar, gmx_bool *bAPolar, t_PolKey *PolKey, t_APolKey *APolKey)	{
	warninp_t wi;
	gmx_bool bAllowWarnings=FALSE;
	int maxwarning = 99;
	t_inpfile *inp;
	int ninp;
	const char *tmp;


	wi = init_warning(bAllowWarnings,maxwarning);
	inp = read_inpfile(fnMdp, &ninp, NULL, wi);

	char *ptr[MAXPTR];
	int n;

	//To check for polar solvation calculation
	STYPE ("polar",                 polar,                  NULL);
	n= str_nelem(polar,MAXPTR,ptr);
	if ( (n>1) || ((strcmp(ptr[0],"yes")!=0) && (strcmp(ptr[0],"no")!=0)) )
		gmx_fatal(FARGS, "Polar should be yes or no\n");
	if(strcmp(ptr[0],"yes")==0)
		*bPolar = TRUE;

	//To check for APolar calculation
	STYPE ("apolar",               apolar,           NULL);
	n= str_nelem(apolar,MAXPTR,ptr);
	if ( (n>1) || ((strcmp(ptr[0],"yes")!=0) && (strcmp(ptr[0],"no")!=0)) )
		gmx_fatal(FARGS, "Polar should be yes or no\n");
	if(strcmp(ptr[0],"yes")==0)
		*bAPolar = TRUE;


	if(*bPolar)
	{
		//Psize keywords
		RTYPE ("cfac", 			PolKey->cfac, 		2);
		RTYPE ("gridspace", 	PolKey->gridspace, 	0.5);
		RTYPE ("gmemceil", 		PolKey->gmemceil, 	400);
		RTYPE ("fadd",			PolKey->fadd, 		20);
		RTYPE ("ofrac",			PolKey->ofrac, 		0.1);

		//Polar Keywords
		RTYPE ("pcharge",       PolKey->pcharge,    	0);
		RTYPE ("ncharge",   	PolKey->ncharge,    	0);
		RTYPE ("prad", 	 	 	PolKey->prad,    	0);
		RTYPE ("nrad",   		PolKey->nrad,    	0);
		RTYPE ("pconc",   		PolKey->pconc,    	0);
		RTYPE ("nconc",   		PolKey->nconc,    	0);
		RTYPE ("pdie",   		PolKey->pdie,    	4);
		RTYPE ("sdie",   		PolKey->sdie,    	78.4);
		RTYPE ("vdie",   		PolKey->vdie,    	1);
		RTYPE ("srad",   		PolKey->srad,    	1.4);
		RTYPE ("swin",   		PolKey->swin,    	0.30);
		RTYPE ("sdens",   		PolKey->sdens,    	10);
		RTYPE ("temp",   		PolKey->temp,    	300);
		EETYPE("srfm",			PolKey->srfm, 		srfm_words);
		EETYPE("chgm",			PolKey->chgm,    	chgm_words);
		EETYPE("bcfl",			PolKey->bcfl,    	bcfl_words);
		EETYPE("PBsolver",		PolKey->pbsolver,   	PBsolver);
	}

	if (*bAPolar)
	{
		//APolar Keywords
		RTYPE ("gamma",         APolKey->gamma,           0.030096);
		RTYPE ("sasconst",     APolKey->sasaconst,           0.0);
		RTYPE ("sasrad",        APolKey->sasrad,          1.4);

		RTYPE ("press",         APolKey->press,           0.0);
		RTYPE ("savconst",      APolKey->savconst,        0.0);
		RTYPE ("savrad",        APolKey->savrad,          1.4);

		APolKey->bWCA = FALSE;
		STYPE("WCA",               wca,           NULL);
		n= str_nelem(wca,MAXPTR,ptr);
		if ( (n>1) || ((strcmp(ptr[0],"yes")!=0) && (strcmp(ptr[0],"no")!=0)) )
			gmx_fatal(FARGS, "Polar should be yes or no\n");
		if(strcmp(ptr[0],"yes")==0)
			APolKey->bWCA = TRUE;

		RTYPE ("wcarad",        APolKey->wcarad,          1.4);
		RTYPE ("bconc",         APolKey->bconc,          0.033428);
		RTYPE ("dpos",          APolKey->dpos,           0.05);
		RTYPE ("APsdens",       APolKey->sdens,          200);
		STYPE ("grid",          grid,                    NULL);
		EETYPE("APsrfm",        APolKey->srfm,          APsrfm_words);
		RTYPE ("APswin",        APolKey->swin,          0.2);
		RTYPE ("APtemp",        APolKey->temp,          300);

		n= str_nelem(grid,MAXPTR,ptr);
		if ((n<3) || (n>3))
			gmx_fatal(FARGS, "Number of grid point should be three\n");
		APolKey->grid[0] = strtod(ptr[0],NULL);
		APolKey->grid[1] = strtod(ptr[1],NULL);
		APolKey->grid[2] = strtod(ptr[2],NULL);
	}
	return 0;
}
