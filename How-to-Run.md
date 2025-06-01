---
layout: page
---

## This is now deprecated. [Latest version of g_mmpbsa is now HERE](https://g-mmpbsa.readthedocs.io/). 

### Executing g_mmpbsa

g_mmpbsa is a console application which is executed from terminal/console by command options similar to other GROMACS module. All the input options should be provided on command line depending on the type of calculation. To calculate for single atomic group (e.g. only protein, only ligand, only specific region of protein etc.) one group is required to choose from index file, whereas for complex (e.g. protein-ligand, protein-protein etc.) two separate group is required to choose from index file.

<p>To get extensive details of usage with options, please click on <a href="Usage.html">Usage</a>.</p>

**Only molecular mechanics (vdw and electrostatic) vacuum energy with energy decomposition**

    g_mmpbsa -f traj.xtc -s topol.tpr -n index.ndx -mme -mm energy_MM.xvg -decomp -mmcon contrib_MM.dat
     
     
**Only polar solvation energy with energy decomposition**

    g_mmpbsa -f traj.xtc -s topol.tpr -i mmpbsa.mdp -n index.ndx -nomme -pbsa -decomp -pol polar.xvg -pcon contrib_pol.dat

An example [mmpbsa.mdp](https://github.com/RashmiKumari/g_mmpbsa/blob/master/test/polar_orig/mmpbsa.mdp){:target="_blank"} is provided in `g_mmpbsa/test/polar_orig`.

**Only non-polar solvation energy with energy decomposition**

    g_mmpbsa -f traj.xtc -s topol.tpr -i mmpbsa.mdp -n index.ndx -nomme -pbsa -decomp -apol apolar.xvg -apcon contrib_apol.dat

*   An example [mmpbsa.mdp](https://github.com/RashmiKumari/g_mmpbsa/blob/master/test/sasa_orig/mmpbsa.mdp){:target="_blank"} for SASA model is provided in `g_mmpbsa/test/sasa_orig`.
*   An example [mmpbsa.mdp](https://github.com/RashmiKumari/g_mmpbsa/blob/master/test/sav_orig/mmpbsa.mdp){:target="_blank"} for SAV model is provided in `g_mmpbsa/test/sav_orig`.
*   An example [mmpbsa.mdp](https://github.com/RashmiKumari/g_mmpbsa/blob/master/test/wca_orig/mmpbsa.mdp){:target="_blank"} for WCA model is provided in `g_mmpbsa/test/wca_orig`.

**All energetic term with energy decomposition**

    g_mmpbsa -f traj.xtc            -s topol.tpr \
             -i mmpbsa.mdp          -n index.ndx -pbsa \
             -mm energy_MM.xvg      -pol polar.xvg \
             -apol apolar.xvg       -decomp \
             -mmcon contrib_MM.dat  -pcon contrib_pol.dat \
             -apcon contrib_apol.dat 
              

**NOTE:** Please monitor the RAM because combined Molecular-Mechanics and Polar-Solvation energy calculation may require GBs of memory.

***

### <a name="extrnAPBS"></a>Using g_mmpbsa with external APBS

When g_mmpbsa is compiled and installed without linking with APBS libraries, an external APBS executable could be used. An environment variable `$APBS` should be defined for path to APBS executable before using g_mmpbsa.

    export APBS=/usr/local/bin/apbs
    
***


### <a name="openmp"></a> Parallel computation using `OpenMP`

g_mmpbsa inherits OpenMP parallel comuputation implemented in APBS 1.2.x, 1.3.x and 1.4.x. This parallel support is enabled by default and all processors/core will be used during runtime.

To control the usage of processors/cores, number of threads can be changed by defining environment variable. For example in bash, one can write following command:

    export OMP_NUM_THREADS=X

where `X` is number of core/processors.

***

### <a name="mpirun"></a> Parallel computation using `mpirun`

Although g_mmpbsa does not support `mpirun`, it can use external APBS with mpirun. This external APBS should be compiled and linked against `MPI`. Environment variable `$APBS` should be defined before executing g_mmpbsa. For example:

    export APBS="mpirun -np 8 apbs"
    
Afterwards, g_mmpbsa uses this `apbs` with `mpirun`, and calculates polar-solvation energy using 8 processors.

**Note:**{: style="color: red"} Presently, MPI support in APBS-1.4.1 is broken, therefore, only use APBS-1.3.x versions. 

***
