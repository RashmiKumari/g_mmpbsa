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

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#include "g_mmpbsa.h"

void substring(int start, int stop, const char *src, char *dest)
{
   sprintf(dest,"%.*s", stop - start, &src[start]);
}

void Bondi(t_AtomProp *bondi)	{
	int i =0;
	char *atomtype[] = {"O", "o", "S","s", "N", "n", "C", "c", "H", "h", "P", "p", "F", "f",
			            "I", "i", "CL", "cl", "Cl", "Br", "BR", "br", "CA", "CB", "CC", "CN",
			            "CR", "CV", "CW", "C*", "cd", "cc", "ca", "HA", "H4", "H5", "ha", "h4", "h5"};
    real radius[] = {1.520, 1.520, 1.830, 1.830, 1.550, 1.550, 1.700, 1.700, 1.200, 1.200, 1.800, 1.800, 1.470, 1.470,
    	  	        2.060, 2.060, 1.770, 1.770, 1.770, 1.920, 1.920, 1.20, 1.770, 1.770, 1.770, 1.770,
    	  	        1.770, 1.770, 1.770, 1.770, 1.770, 1.770, 1.770, 1.000, 1.000, 1.000, 1.000, 1.000, 1.00};
	bondi->n = 38;
	bondi->atomtype = (char**) malloc(sizeof(char*)*(bondi->n+1));
	bondi->value = (real*) malloc(sizeof(real)*(bondi->n+1));

	for(i=0;i<=bondi->n;i++)	{
		 bondi->atomtype[i] = (char*) malloc(sizeof(atomtype[i]));
		 sprintf(bondi->atomtype[i],"%s",atomtype[i]);
		 bondi->value[i] = radius[i];
	}
}

void mBondi(t_AtomProp *mbondi)	{
	int i=0;
	char *atomtype[] ={ "O", "o", "S", "s", "N", "n", "C", "c", "H", "h", "P", "p", "F", "f", "I", "i",
						"CL", "Cl", "cl", "BR", "Br", "br", "HC", "CA", "CB", "CC", "CN", "CR", "CV",
						"CW", "C*", "cd", "cc", "ca", "HA", "H4", "H5", "ha", "h4", "h5", "hc", "HN",
						"hn", "HO", "ho", "HS", "hs", "HP", "hp"};
	real radius[] = {1.500,1.5000,1.800,1.800,1.550,1.550,1.700,1.700,1.200,1.200,1.850,1.850,1.470,1.470,1.980,1.980,
			          1.770,1.770,1.770,1.850,1.850,1.850,1.300,1.770,1.770,1.770,1.770,1.770,1.770,
			          1.770,1.770,1.770,1.770,1.770,1.000,1.000,1.000,1.000,1.000,1.000,1.30,1.30,
			          1.300,0.80,0.80,0.80,0.80,1.30,1.30};
	mbondi->n = 48;
	mbondi->atomtype = (char**) malloc(sizeof(char*)*(mbondi->n+1));
	mbondi->value = (real*) malloc(sizeof(real)*(mbondi->n+1));

	for(i=0;i<=mbondi->n;i++)	{
		 mbondi->atomtype[i] = (char*) malloc(sizeof(atomtype[i]));
		 sprintf(mbondi->atomtype[i],"%s",atomtype[i]);
		 mbondi->value[i] = radius[i];
	}
}

void mBondi2 (t_AtomProp *mbondi2)	{
	int i =0;
	char *atomtype[] = {"O", "o", "S", "s", "N", "n", "C", "c", "H", "h", "P", "p", "F", "I", "f", "i",
	                    "CL", "cl", "Cl", "BR", "br", "Br", "HC", "hc", "HN", "hn", "HO", "ho", "HS",
	                    "hs", "HP", "hp"};
	real radius[] = {1.500, 1.500, 1.800, 1.800, 1.550, 1.550, 1.700, 1.700, 1.200, 1.200, 1.850, 1.850, 1.470, 1.980, 1.470, 1.980,
	 	              1.770, 1.770, 1.770, 1.850, 1.850, 1.850, 1.300, 1.300, 1.300, 1.300, 0.800, 0.800, 0.800,
	 	              0.800, 1.300, 1.300};
	mbondi2->n = 31;
	mbondi2->atomtype = (char**) malloc(sizeof(char*)*(mbondi2->n+1));
	mbondi2->value = (real*) malloc(sizeof(real)*(mbondi2->n+1));

	for(i=0;i<=mbondi2->n;i++)	{
		 mbondi2->atomtype[i] = (char*) malloc(sizeof(atomtype[i]));
		 sprintf(mbondi2->atomtype[i],"%s",atomtype[i]);
		 mbondi2->value[i] = radius[i];
	}
}

real GetRad (char *atomtype, t_AtomProp *radtype, real rvdw)	{
	char *a, *b;
	int i =0, done=0;
	real r =0;
	a = (char*) malloc(sizeof(char)*1);
	b = (char*) malloc(sizeof(char)*2);

	substring(0,1, atomtype,a);
	for(i=0;i<=radtype->n;i++)	{
		if(strcmp(a,radtype->atomtype[i])==0)	{
			r = radtype->value[i];
			done=1;
		}
	}

	if(atomtype[1])	{
		substring(0,2, atomtype,b);
		for(i=0;i<=radtype->n;i++)	{
			if(strcmp(b,radtype->atomtype[i])==0)	{
				r = radtype->value[i];
				done = 1;
			}
		}
	}

	if(done==0)	{
		printf("\nWARNING: Radius for atomtype \"%s\" not found. Assigning radius to %3.2lf \n",atomtype,rvdw);
		return rvdw;
	}
	else
		return r;
}
