---
layout: page
---

### Input Keywords for Polar Solvation Energy


##### polar

    polar          = yes

This will allow the calculation of polar solvation energy. The value can change to “no” if one doesn’t want to do polar calculation.

***
    
##### cfac
 
    cfac           = 2

The factor by which molecular dimensions should expand to get a coarse grid dimensions. For detail see help of a APBS script `apbs-1.3-source/tools/manip/psize.py` (psize.py -h).

***

##### fadd

    fadd           = 20

The amount (in Å) to add to molecular dimensions to get a fine grid dimensions. For detail see help of a APBS script `apbs-1.3-source/tools/manip/psize.py` (psize.py -h).

***

##### gridspace

    gridspace      = 0.2
    
It specifies the value (in Å)  for  fine grid spacing. For detail see help of a APBS script `apbs-1.3-source/tools/manip/psize.py` (psize.py -h).

***
    
##### gmemceil

    gmemceil       = 4000
    
Usage: Sets memory (in MB) which will be used per-processor for a calculation. For detail see help of a APBS script `apbs-1.3-source/tools/manip/psize.py` (psize.py -h).

***

##### PBsolver

    PBsolver        = npbe

This specifies whether linear or nonlinear Poisson Boltzmann equation should be solved. The accepted keywords are `lpbe` and `npbe` for linear and traditional non-linear PB equation, respectively. The effects of different value on the polar calculation is checked during this implementation.

* * *

##### mg-type

    mg-type        = mg-auto

How multigrid PB calculation should be performed? **Accepted keywords:** `mg-auto` and `mg-para`.

[**mg-auto**](http://www.poissonboltzmann.org/docs/elec-calcs/#mgauto){:target="_blank"} will perform automatically-configured sequential focusing multigrid PB calculation.


[**mg-para**](http://www.poissonboltzmann.org/docs/elec-calcs/#mgpara){:target="_blank"} will perform automatically-configured parallel focusing multigrid PB calculation. **Note:** This keyword only works with external APBS executable and mpirun.

***


##### pcharge

    pcharge        = 1

The charge of positive ions in bulk solution.

***

##### prad

    prad           = 0.95
    
Radius of positive ions.

***

##### pconc

    pconc          = 0.150

Concentration of positive ion.

***

##### ncharge

    ncharge        = -1

The charge of negative ions in bulk solution.

***

##### nrad

    nrad           = 1.81

Radius of negative ion.

***

##### nconc

    nconc          = 0.150

Concentration of negative ion.

***

##### pdie

    pdie           = 4
    
The value of solute dielectric constant. This can be change depending on the solute used for calculation. For highly charged solute high dielectric value will produce more accurate polar solvation energy.

***

##### sdie

    sdie           = 80

The value of solvent dielectric constant. 

***

##### vdie

    vdie           = 1

The value of vacuum dielectric constant.

***

##### srad

    srad           = 1.4 
    
This specify the radius (in Å)  of solvent molecules. This is used in case of probe-based surface definition. **For more details, see Elec keyword `srad` [here][elec-keywords]{:target="_blank"}**.

***

##### swin

    swin           = 0.30

This specify the value for cubic spline window for spline-based surface definitions. Not used when probe-based surface are used in calculation. **For more details, see Elec keyword `swin` [here][elec-keywords]{:target="_blank"}**.

***

##### srfm 

    srfm           = smol

This specify the model used to construct the dielectric and ion-accessibility coefficients. The accepted keywords are `mol`, `smol`, `spl2` and `spl4` and it may affect the polar energy calculation. **For more details, see Elec keyword `srfm` [here][elec-keywords]{:target="_blank"}**.

***

##### sdens

    sdens          = 10

Specify the number of grid points per Å<sup>2</sup> for constructing the molecular surface or solvent accessible surface. Not taken in consideration when `srad = 0.0` or `srfm = spl2`. **For more details, see Elec keyword `sdens` [here][elec-keywords]{:target="_blank"}**.

***

##### temp

    temp           = 300

This specify the temperature used for Poisson-Boltzmann calculation. **For more details, see Elec keyword `temp` [here][elec-keywords]{:target="_blank"}**.

***

##### chgm

    chgm            = spl4

This specify the method used to map the biomolecular point charges to the grid for a multigrid Poisson-Boltzmann calculation. The accepted keywords are `spl0`, `spl2` and `spl4`. The effects of these keywords on energy are not tested in this implementation. **For more details, see Elec keyword `chgm` [here][elec-keywords]{:target="_blank"}**.

***

##### bcfl

    bcfl            = mdh

It specifies the type of boundary conditions used to solve the Poisson-Boltzmann equation. The accepted keywords are `zero`, `sdh`, `mdh`, `focus`, and `map`. However, use of `focus`, and `map` will terminate g_mmpbsa with error. The change in bcfl keywords may affect the polar energy calculation. The effects of these keywords on energy are not tested in this implementation. **For more details, see Elec keyword `bdfl` [here][elec-keywords]{:target="_blank"}**.

***


[elec-keywords]: http://www.poissonboltzmann.org/docs/apbs-overview/#elec