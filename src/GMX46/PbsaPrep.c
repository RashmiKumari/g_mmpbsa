/**
 * @file GMX46/PbsaPrep.c
 * @brief Definition of routines to setup PBSA calculations
 * @ingroup PBSA_PREP
 * @author Rashmi Kumari, Rajendra Kumar and Andrew Lynn
 * @attention
 * @verbatim
 *
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
 * @endverbatim
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gromacs/statutil.h"
#include "gromacs/typedefs.h"
#include "gromacs/smalloc.h"
#include "gromacs/vec.h"
#include "gromacs/tpxio.h"

#include "g_mmpbsa.h"

//To modify residue name as many time residue have same name in topology
//but have different structure.
void modResname(t_topology *top, char **modresname)	{

  int n = top->atoms.nr;
  int i = 0;
  int resnmrCur, resnmrTer;
  char *atomname, *tmpresname, *resname, *atomtype;

  for (i = 0; i < n; i++)
    {
      resname = *(top->atoms.resinfo[top->atoms.atom[i].resind].name);
      modresname[i] = strdup(resname);
    }

  for (i = 0; i < n; i++)
    {
      atomname = *(top->atoms.atomname[i]);
      atomtype = *(top->atoms.atomtype[i]);
      resname = *(top->atoms.resinfo[top->atoms.atom[i].resind].name);
      resnmrCur = top->atoms.resinfo[top->atoms.atom[i].resind].nr;
      //To modify name of two different HISA/HID and HISB/HIE

      if (((strcmp(atomname, "HD1") == 0) || (strcmp(atomname, "HE2") == 0)) && (strcmp(resname, "HIS") == 0))
        {
          int j = 0, k = 0;
          char tmpresname[16];

          if (strcmp(atomname, "HD1") == 0)
            sprintf(tmpresname, "HID");
          if (strcmp(atomname, "HE2") == 0)
            sprintf(tmpresname, "HIE");

          resnmrTer = top->atoms.resinfo[top->atoms.atom[i].resind].nr;
          while (1)
            {
              if (i - j < 0)
                break;
              resnmrCur = top->atoms.resinfo[top->atoms.atom[i - j].resind].nr;
              if (resnmrTer != resnmrCur)
                break;
              j++;
            }
          j = j - 1;
          resnmrCur = top->atoms.resinfo[top->atoms.atom[i - j].resind].nr;

          while (1)
            {
              if (i - j + k >= n)
                break;
              resnmrCur = top->atoms.resinfo[top->atoms.atom[i - j + k].resind].nr;
              if (resnmrTer != resnmrCur)
                break;
              modresname[i - j + k] = strdup(tmpresname);
              k++;
            }
          i = i - j + k - 1;
        }
      // END //

      //To modify name of CYS linked through di-sulfide linkage (CYS -CYS)
      if (strcmp(resname, "CYS") == 0 && strcmp(atomtype, "S") == 0)
        {
          int j = 0, k = 0;
          char tmpresname[16];

          sprintf(tmpresname, "CYX");
          resnmrTer = top->atoms.resinfo[top->atoms.atom[i].resind].nr;
          while (1)
            {
              if (i - j < 0)
                break;
              resnmrCur = top->atoms.resinfo[top->atoms.atom[i - j].resind].nr;
              if (resnmrTer != resnmrCur)
                break;
              j++;
            }
          j = j - 1;
          resnmrCur = top->atoms.resinfo[top->atoms.atom[i - j].resind].nr;

          while (1)
            {
              if (i - j + k >= n)
                break;
              resnmrCur = top->atoms.resinfo[top->atoms.atom[i - j + k].resind].nr;
              if (resnmrTer != resnmrCur)
                break;
              modresname[i - j + k] = strdup(tmpresname);
              k++;
            }
          i = i - j + k - 1;
        }
      // END //
    }

  for (i = 0; i < n; i++)
    {
      atomname = *(top->atoms.atomname[i]);
      atomtype = *(top->atoms.atomtype[i]);
      resname = *(top->atoms.resinfo[top->atoms.atom[i].resind].name);
      resnmrCur = top->atoms.resinfo[top->atoms.atom[i].resind].nr;

      //To modify terminal residue name
      if (   ((strcmp(atomname, "OC1") == 0) || (strcmp(atomname, "H1") == 0))
          && (
            		 (strcmp(resname, "ARG") == 0) || (strcmp(resname, "HIS") == 0)
                  || (strcmp(resname, "LYS") == 0) || (strcmp(resname, "ASP") == 0)
                  || (strcmp(resname, "GLU") == 0) || (strcmp(resname, "SER") == 0)
                  || (strcmp(resname, "THR") == 0) || (strcmp(resname, "ASN") == 0)
                  || (strcmp(resname, "GLN") == 0) || (strcmp(resname, "CYS") == 0)
                  || (strcmp(resname, "GLY") == 0) || (strcmp(resname, "PRO") == 0)
                  || (strcmp(resname, "ALA") == 0) || (strcmp(resname, "VAL") == 0)
                  || (strcmp(resname, "ILE") == 0) || (strcmp(resname, "LEU") == 0)
                  || (strcmp(resname, "MET") == 0) || (strcmp(resname, "PHE") == 0)
                  || (strcmp(resname, "TYR") == 0) || (strcmp(resname, "TRP") == 0)
                  || (strcmp(resname, "ASH") == 0) || (strcmp(resname, "GLH") == 0)
                  || (strcmp(resname, "LYN") == 0) || (strcmp(resname, "CYM") == 0)
                  || (strcmp(resname, "HIP") == 0) || (strcmp(resname, "CYX") == 0)
                  || (strcmp(resname, "HID") == 0) || (strcmp(resname, "HIE") == 0)
             )
         )
        {
          int j = 0, k = 0;
          char tmpresname[16];

          if (strcmp(atomname, "OC1") == 0)
            sprintf(tmpresname, "C%s", modresname[i]);
          if (strcmp(atomname, "H1") == 0)
            sprintf(tmpresname, "N%s", modresname[i]);

          resnmrTer = top->atoms.resinfo[top->atoms.atom[i].resind].nr;
          while (1)
            {
              if (i - j < 0)
                break;
              resnmrCur = top->atoms.resinfo[top->atoms.atom[i - j].resind].nr;
              if (resnmrTer != resnmrCur)
                break;
              j++;
            }
          j = j - 1;
          resnmrCur = top->atoms.resinfo[top->atoms.atom[i - j].resind].nr;

          while (1)
            {
              if (i - j + k >= n)
                break;
              resnmrCur =
                  top->atoms.resinfo[top->atoms.atom[i - j + k].resind].nr;
              if (resnmrTer != resnmrCur)
               break;
              modresname[i - j + k] = strdup(tmpresname);
              k++;
            }
          i = i - j + k - 1;
        }
    }
  // END //
}

int assignQR(t_topology *top, atom_id *index, int isize, int eRadType,  t_AtomProp *radtype, real rvdw) {

  int i, itype;
  double c6, c12;
  int ntype = top->idef.atnr;
  real sig6, rad;
  char *atomtype, *resname, *atomname;

  snew(top->atoms.pdbinfo, top->atoms.nr);

  for (i = 0; (i < isize); i++)
    {
      itype = top->atoms.atom[index[i]].type;
      c12 = top->idef.iparams[itype * ntype + itype].lj.c12;
      c6 = top->idef.iparams[itype * ntype + itype].lj.c6;
      if ((c6 != 0) && (c12 != 0))
        {
          sig6 = c12 / c6;
          rad = 0.5 * pow(sig6, 1.0 / 6.0);
        }
      else
        rad = 0.000;
      //printf("%f\n", atoms->pdbinfo[i].occup);
      rad *= 10; //Conversion of nano meter to angstroms

      top->atoms.pdbinfo[index[i]].occup = top->atoms.atom[index[i]].q;
      top->atoms.pdbinfo[index[i]].bfac = rad;

      atomtype = *(top->atoms.atomtype[index[i]]);
      atomname = *(top->atoms.atomname[index[i]]);
      resname = *(top->atoms.resinfo[top->atoms.atom[index[i]].resind].name);

      if ((eRadType == eBondi) || ( eRadType == eMbondi) || (eRadType == eMbondi2) )
        top->atoms.pdbinfo[index[i]].bfac = GetRad(atomtype, radtype, rvdw);
      if (eRadType == eFF)
        top->atoms.pdbinfo[index[i]].bfac = rad;

      //printf("%s\t%s\t%s\t",atomname,atomtype,resname);
      //printf("%f\t%f\n",top->atoms.pdbinfo[i].occup,top->atoms.pdbinfo[i].bfac);
    }
  return 0;
}

int makePQR(t_topology *top, atom_id *index, int isize, int ePBC, matrix box, rvec x[], char **modresname, char *fnPQR)	{

  FILE *fPQR;
  int i;
  char *resname, *atomname;
  int resnmr;
  real q, r;
  static const char *pdb = "%-6s%5u  %-5.4s%4.4s  %4d    %8.3f%8.3f%8.3f";
  //set_pdb_wide_format(TRUE);
  fPQR = ffopen(fnPQR, "w");
  //write_pdbfile_indexed(fPQR,NULL,&(top->atoms),x,ePBC,box,' ',-1,isize,index,NULL,TRUE);
  for (i = 0; i < isize; i++)
    {
      atomname = *(top->atoms.atomname[index[i]]);
      resname = *(top->atoms.resinfo[top->atoms.atom[index[i]].resind].name);
      resnmr = top->atoms.resinfo[top->atoms.atom[index[i]].resind].nr;
      q = top->atoms.pdbinfo[index[i]].occup;
      r = top->atoms.pdbinfo[index[i]].bfac;
                         //%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f";
      fprintf(fPQR, pdb, "ATOM", i + 1, atomname, modresname[index[i]], resnmr, 10 * x[index[i]][XX], 10 * x[index[i]][YY], 10 * x[index[i]][ZZ]);
      fprintf(fPQR, "%8.3f%8.3f\n", q, r);
    }
  ffclose(fPQR);
  return 0;
}

void SasvRad(t_topology *top, int isize, atom_id *index, real ***radius, t_APolKey APolKey)
{
  gmx_bool *bRad;
  int i = 0;
  int natoms;
  natoms = top->atoms.nr;
  real **localRadius;

  snew(localRadius,3);
  for(i=0;i<3;i++)
	  snew(localRadius[i],natoms);

  snew(bRad, natoms);
  for (i = 0; i < isize; i++)
    {
      bRad[index[i]] = TRUE;
    }


  for (i = 0; i < natoms; i++)
    {
      if (bRad[i])	{
        localRadius[eSASRAD][i] = (top->atoms.pdbinfo[i].bfac/10) + (APolKey.sasrad/10);
        localRadius[eSAVRAD][i] = (top->atoms.pdbinfo[i].bfac/10) + (APolKey.savrad/10);
        localRadius[eWCARAD][i] = (top->atoms.pdbinfo[i].bfac/10) + (APolKey.wcarad/10);
      }
      else	{
    	  localRadius[eSASRAD][i] = APolKey.sasrad;
    	  localRadius[eSAVRAD][i] = APolKey.savrad;
    	  localRadius[eWCARAD][i] = APolKey.wcarad;
      }
    }
  *radius = localRadius;
}

void ApbsParamAPol(t_topology *top, atom_id *index, int isize, char **modresname, char *fnApbsParamAPol)
{
  char *atomname, *resname, *atomtype;
  real q, m;
  int n, resid, itype, ntype = top->idef.atnr;
  real sig6, vdw, sigma, epsilon;
  double c6, c12;
  FILE *fOut;
  fOut = ffopen(fnApbsParamAPol, "w");

  for (n = 0; n < isize; n++)
    {
      itype = top->atoms.atom[index[n]].type;
      c12 = top->idef.iparams[itype * ntype + itype].lj.c12;
      c6 = top->idef.iparams[itype * ntype + itype].lj.c6;

      //c6 = 4 * epsilon * sigma^6
      //c12 = 4 * epsilon * sigma^12
      //C12/c6  = sigma^2

      sig6 = c12 / c6;
      if (c6 != 0)
        vdw = 10 * 0.5 * pow(sig6, 1.0 / 6.0);
      else
        vdw = 0;

      sigma = pow(sig6, 1.000 / 6.000);
      epsilon = c6 / (4 * sig6);
      if (vdw == 0)
        epsilon = 0.000;

      atomname = *(top->atoms.atomname[index[n]]);
      atomtype = *(top->atoms.atomtype[index[n]]);
      resname = *(top->atoms.resinfo[top->atoms.atom[index[n]].resind].name);
      resid = top->atoms.resinfo[top->atoms.atom[index[n]].resind].nr;
      q = top->atoms.atom[index[n]].q;
      m = top->atoms.atom[index[n]].m;
      //printf("%10d %6s %6s %6s %10.3f %10.3f %10.3f %10.3f %12.5g %12.5g\n",n+1, atomname, atomtype, resname,10*x[index[n]][0],10*x[index[n]][1],10*x[index[n]][2],q,c6,c12);
      fprintf(fOut, "%-6s     %-6s     %5.3lf     %5.3lf     %10.5lf\n",
          modresname[index[n]], atomname, q, vdw, epsilon);
      //printf("%4.3f",q);
      //printf("%4.3f",m);
      //printf("%4.3f%4.3f4.3%4.3f",x[index[n]][0],x[index[n]][1],x[index[n]][2]);
      //printf("%4.5g",c6);
      //printf("%4.5g\n",c12);
    }
  fprintf(fOut, "WAT     OW      -0.834        1.575305  6.36386e-01\n");
  fprintf(fOut, "WAT     HW1      0.417        0.000000  0.000000000\n");
  fprintf(fOut, "WAT     HW2      0.417        0.000000  0.000000000");
  ffclose(fOut);

}

int polarInAPBS(t_PolKey *param, char *fnPQR, char *fnPolAPBS, gmx_bool bDECOMP)	{

  FILE *fIn;
  fIn = ffopen(fnPolAPBS, "w");

  fprintf(fIn, "read\n    mol pqr %s\nend\n", fnPQR);

  if(param->mg_type == mg_para)	{
	  fprintf(fIn, "\nelec name mol1\n    mg-para\n");
	  fprintf(fIn, "    dime  %d %d %d\n", param->dime[XX], param->dime[YY], param->dime[ZZ]);
	  fprintf(fIn, "    pdime  %d %d %d\n", param->pdime[XX], param->pdime[YY], param->pdime[ZZ]);
	  fprintf(fIn, "    ofrac %g\n", param->ofrac);
  }
  else	{
	  fprintf(fIn, "\nelec name mol1\n    mg-auto\n");
	  fprintf(fIn, "    dime  %d %d %d\n", param->dime[XX], param->dime[YY], param->dime[ZZ]);
  }
  fprintf(fIn, "    cglen %6.3lf %6.3lf %6.3lf\n", param->cglen[XX], param->cglen[YY], param->cglen[ZZ]);
  fprintf(fIn, "    fglen %6.3lf %6.3lf %6.3lf\n", param->fglen[XX], param->fglen[YY], param->fglen[ZZ]);
  fprintf(fIn, "    cgcent %6.3lf %6.3lf %6.3lf\n", param->cgcent[XX], param->cgcent[YY], param->cgcent[ZZ]);
  fprintf(fIn, "    fgcent %6.3lf %6.3lf %6.3lf\n", param->fgcent[XX], param->fgcent[YY], param->fgcent[ZZ]);
  fprintf(fIn, "    mol 1\n");
  fprintf(fIn, "    %s\n", PBsolver[param->pbsolver]);
  fprintf(fIn, "    bcfl %s\n", bcfl_words[param->bcfl]);
  fprintf(fIn, "    ion %.1g %.3g %.3g\n", param->pcharge, param->pconc, param->prad);
  fprintf(fIn, "    ion %.1g %.3g %.3g\n", param->ncharge, param->nconc, param->nrad);
  fprintf(fIn, "    pdie %g\n", param->pdie);
  fprintf(fIn, "    sdie %g\n", param->sdie);
  fprintf(fIn, "    srfm %s\n", srfm_words[param->srfm]);
  fprintf(fIn, "    chgm %s\n", chgm_words[param->chgm]);
  fprintf(fIn, "    sdens %g\n", param->sdens);
  fprintf(fIn, "    srad %g\n", param->srad);
  fprintf(fIn, "    swin %g\n", param->swin);
  fprintf(fIn, "    temp %g\n", param->temp);
  if(bDECOMP)
	  fprintf(fIn, "    calcenergy comps\n");
  else
	  fprintf(fIn, "    calcenergy total\n");
  fprintf(fIn, "end\n");

  if(param->mg_type == mg_para)	{
	  fprintf(fIn, "\nelec name mol2\n    mg-para\n");
	  fprintf(fIn, "    dime  %d %d %d\n", param->dime[XX], param->dime[YY], param->dime[ZZ]);
	  fprintf(fIn, "    pdime  %d %d %d\n", param->pdime[XX], param->pdime[YY], param->pdime[ZZ]);
	  fprintf(fIn, "    ofrac %g\n", param->ofrac);
  }
  else	{
	  fprintf(fIn, "\nelec name mol2\n    mg-auto\n");
	  fprintf(fIn, "    dime  %d %d %d\n", param->dime[XX], param->dime[YY], param->dime[ZZ]);
  }

  fprintf(fIn, "    cglen %6.3lf %6.3lf %6.3lf\n", param->cglen[XX], param->cglen[YY], param->cglen[ZZ]);
  fprintf(fIn, "    fglen %6.3lf %6.3lf %6.3lf\n", param->fglen[XX], param->fglen[YY], param->fglen[ZZ]);
  fprintf(fIn, "    cgcent %6.3lf %6.3lf %6.3lf\n", param->cgcent[XX], param->cgcent[YY], param->cgcent[ZZ]);
  fprintf(fIn, "    fgcent %6.3lf %6.3lf %6.3lf\n", param->fgcent[XX], param->fgcent[YY], param->fgcent[ZZ]);
  fprintf(fIn, "    mol 1\n");
  fprintf(fIn, "    %s\n", PBsolver[param->pbsolver]);
  fprintf(fIn, "    bcfl %s\n", bcfl_words[param->bcfl]);
  fprintf(fIn, "    ion %.1g %.3g %.3g\n", param->pcharge, param->pconc, param->prad);
  fprintf(fIn, "    ion %.1g %.3g %.3g\n", param->ncharge, param->nconc, param->nrad);
  fprintf(fIn, "    pdie %g\n", param->pdie);
  fprintf(fIn, "    sdie %g\n", param->vdie);
  fprintf(fIn, "    srfm %s\n", srfm_words[param->srfm]);
  fprintf(fIn, "    chgm %s\n", chgm_words[param->chgm]);
  fprintf(fIn, "    sdens %g\n", param->sdens);
  fprintf(fIn, "    srad %g\n", param->srad);
  fprintf(fIn, "    swin %g\n", param->swin);
  fprintf(fIn, "    temp %g\n", param->temp);
  if(bDECOMP)
	  fprintf(fIn, "    calcenergy comps\n");
  else
	  fprintf(fIn, "    calcenergy total\n");
  fprintf(fIn, "end\n");
  fprintf(fIn, "print elecEnergy mol1 - mol2 end\n");
  fprintf(fIn, "quit\n");
  ffclose(fIn);

  return 0;
}

int APolarInAPBS(t_APolKey *APolKey, char *fnPQR, char *fnAPolAPBS, char *fnApbsParamAPol)	{
  FILE *fIn;
  fIn = ffopen(fnAPolAPBS, "w");

  fprintf(fIn, "read\n    mol pqr %s\n    parm flat %s\nend\n", fnPQR, fnApbsParamAPol);
  fprintf(fIn, "\nAPOLAR name mol1\n");
  fprintf(fIn, "    mol 1\n");
  fprintf(fIn, "    bconc %g\n", APolKey->bconc);
  fprintf(fIn, "    dpos %g\n", APolKey->dpos);
  fprintf(fIn, "    sdens %g\n", APolKey->sdens);
  fprintf(fIn, "    press %g\n", APolKey->press);
  fprintf(fIn, "    gamma %g\n", APolKey->gamma);
  fprintf(fIn, "    grid %g %g %g\n", APolKey->grid[0],APolKey->grid[1],APolKey->grid[2]);
  fprintf(fIn, "    srad %g\n", APolKey->wcarad);
  fprintf(fIn, "    srfm %s\n", APsrfm_words[APolKey->srfm]);
  fprintf(fIn, "    swin %g\n", APolKey->swin);
  fprintf(fIn, "    temp %g\n", APolKey->temp);
  fprintf(fIn, "    calcenergy total\n");
  fprintf(fIn, "end\n");
  fprintf(fIn, "\nprint apolEnergy mol1 end\n");
  fprintf(fIn, "quit\n");
  ffclose(fIn);

  return 0;
}
