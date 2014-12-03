---
layout: page
---


<a name="top"></a>

### [g_mmpbsa](#g_mmpbsa)

### [MmPbSaStat.py](#mmpbsastatpy)

### [MmPbSaDecomp.py](#mmpbsadecomppy)

### [energy2bfac](#energy2bfac)

### <a name="g_mmpbsa"></a> g_mmpbsa

g_mmpbsa calculates binding energy of biomolecular associations like protein-protein, protein-ligand protein-DNA etc using MM-PBSA. It gives the different component of energy term in separate file so that user will have choice to have either MM, PB and SA energy values or all energies according to their objective. The tool also gives residue wise contribution to total binding energy which will provide information about important contributing residues to the molecular association.

    Option     Filename  Type         Description
    ------------------------------------------------------------
      -f       traj.xtc  Input        Trajectory: xtc trr trj gro g96 pdb cpt
      -s      topol.tpr  Input        Run input file: tpr tpb tpa
      -i     grompp.mdp  Input, Opt.  grompp input file with MD parameters
      -n      index.ndx  Input, Opt.  Index file
     -mm  energy_MM.xvg  Output, Opt. xvgr/xmgr file
    -pol      polar.xvg  Output, Opt. xvgr/xmgr file
    -apol    apolar.xvg  Output, Opt. xvgr/xmgr file
    -mmcon contrib_MM.dat  Output, Opt. Generic data file
    -pcon contrib_pol.dat  Output, Opt. Generic data file
    -apcon contrib_apol.dat  Output, Opt. Generic data file
    

    Option       Type   Value   Description
    ------------------------------------------------------
    -[no]h       bool   yes     Print help info and quit
    -[no]version bool   no      Print version info and quit
    -nice        int    19      Set the nicelevel
    -b           time   0       First frame (ps) to read from trajectory
    -e           time   0       Last frame (ps) to read from trajectory
    -dt          time   0       Only use frame when t MOD dt = first time (ps)
    -tu          enum   ps      Time unit: fs, ps, ns, us, ms or s
    -[no]w       bool   no      View output .xvg, .xpm, .eps and .pdb files
    -xvg         enum   xmgrace  xvg plot formatting: xmgrace, xmgr or none
    -[no]silent  bool   no      Display messages, output and errors from external
                                APBS program. Only works with external APBS
                                program
    -rad         enum   bondi   van der Waal radius type: bondi, mbondi, mbondi2
                                or amber
    -rvdw        real   1       Default van der Waal radius (in nm) if not found
    -[no]mme     bool   yes     To calculate vacuum molecular mechanics energy
    -pdie        real   1       Dielectric constant of solute. Should be same as
                                of polar solvation
    -[no]incl_14 bool   no      Include 1-4 atom-pairs, exclude 1-2 and 1-3 atom
                                pairs during MM calculation. Should be "yes" when
                                groups are bonded with each other.
    -[no]focus   bool   no      To enable focusing on the specfic region of
                                molecule, group of atoms must be provided in
                                index file
    -[no]pbsa    bool   no      To calculate polar and/or non-polar solvation
                                energy
    -ndots       int    24      Number of dots per sphere in the calculation of
                                SASA, more dots means more accuracy
    -[no]diff    bool   yes     Calculate the energy difference between two group
                                otherwise only calculates for one group
    -[no]decomp  bool   no      Decomposition of energy for each residue


#### File Options

    -s topol.tpr

Input tpr/tpx file of molecule.

***
    -f traj.xtc

Input trajectory xtc/trr format file. 

**WARNING:** Trajectory should be PBC corrected and molecule should not be PBC broken. To make molecule whole in trajectory, please follow these links: [PBC](http://www.gromacs.org/Documentation/Terminology/Periodic_Boundary_Conditions) and [trjconv](http://manual.gromacs.org/current/online/trjconv.html).

***
    -n index.ndx

Input atomic index file. User will get choice to select atomic groups.

***
    -i mmpbsa.mdp

Input parameter file for polar and non-polar solvation energy. For more details about accepted keywords and options, follow these two links: [Polar-Solvation Keywords](Polar-solvation-input-parameters.html) and [Non-polar Solvation Keywords](Non-polar-solvation-input-parameters.html).

***
    -mm energy_MM.xvg

van der Waal and electrostatic energy of the selected atom group/s.  

* With `-nodiff` option, only one index group can be selected. In this case, this file contains vacuum MM energy of this selected group. Always **USE**{: style="color: red"} `-incl_14` option for single group calculations.

* By default, two groups can be selected, and this file contains only interaction energy between two groups. Energy of each group and thier complex is not calculated. 

* However, with `-incl_14` option, vacuum MM energy componenets for each group and their complex is calculated. Interaction energy can be calculated later using the provided Python scripts.  

***
    -pol polar.xvg

Polar solvation energy of the selected atoms group/s.

* With `-nodiff` option, only one index group can be selected. In this case, this file contains energy of this selected group. 

* By default, two groups can be selected, and this file contains energy of each group and thier complex. 

***
    -apol apolar.xvg

Non-polar solvation energy of the selected atoms group/s. 

* With `-nodiff` option, only one index group can be selected. In this case, this file contains energy of this selected group. 

* By default, two groups can be selected, and this file contains energy of each group and thier complex. 

***
    -mmcon contrib_MM.dat

Vacuum MM van der Waals and electrostatic energy contribution per residue per frame/snaspshot.

***
    -pcon contrib_pol.dat

Polar solvation energy contribution per residue frame wise.

***
    -apcon contrib_apol.dat

Non-polar solvation energy contribution per residue frame wise.

***

#### Other options

    -diff or -nodiff

`Default: yes` 

* By default, selection of two atom groups will be prompted. For example, atom group A and B is selected by user. Then, third AB group will be automatically created by combining A and B. Subsequently, all energy calculation will be performed for these three atom groups A, B and AB.

* If this option is switched off by `-nodiff`, selection of only one atom group will be prompted and all energy calculation will be performed on this selected atom group.

***
    -mme or -nomme

`Default: yes` 

* By default, van der Waals and electrostatic energy of the selected group/s will be calculated. 

* To prevent calculation of the vacuum MM energy, this option can be switched off using `-nomme`.

***
    -pdie 1

Value of solute dielectric constant in the vacuum electrostatic calculation. It should be similar to that of the polar-solvation energy calculation.

***
    -pbsa or -nopbsa

`Default: no`

To calculate polar or non-polar solvation enerby, use of `-pbsa` option is required, and additionally an input paramaeter is required with `-i` option.

***
    -rad bondi

Three keywords are accepted `bondi`, `mbondi` and `mbondi2` which corresponds to three type of radius discussed in the publication of g_mmpbsa.

***
    -ndots 100

Number of dots per sphere used in the calculation of solvent accessible surface area and volume. Higher will be the number, more will be the accuracy.

***
    -decomp  or  -nodecomp

`Default: no` 

To calculate energetic contribution of each residue to total binding energy, use of `-decomp` option is required.

***
    -silent  or  -nosilent

`Default: Yes` 

When external APBS is used, `-silent` option can be used to suppress all the messages from APBS program. This option will not work when g_mmpbsa is compiled with APBS libraries.  

***

[ ↑ top ↑](#top)

***

### <a name="mmpbsastatpy"> </a> MmPbSaStat.py

This script calculates the average binding energy and its standard deviation/error from the output files, which are obtained from g_mmpbsa.

#### Usage:

    python MmPbSaStat.py [-h] [-mt] [-mf metafile.dat] [-m energy_MM.xvg]
                         [-p polar.xvg] [-a apolar.xvg] [-bs] [-nbs 500]
                         [-of full_energy.dat] [-os summary_energy.dat]
                         [-om meta_energy.dat]

***

#### Description

    -h, --help

show this help message and exit

***
    
    -mt, --multiple

If given, calculates for multiple complexes. Need metafile containing path of energy files

***
    -mf metafile.dat, --metafile metafile.dat
    
Metafile containing path to energy files of each complex in a row obtained from g_mmpbsa in following order: `[MM file] [Polar file] [ Non-polar file]`

***
    -m energy_MM.xvg, --molmech energy_MM.xvg

Vacuum Molecular Mechanics energy file obtained from g_mmpbsa.

***
    -p polar.xvg, --polar polar.xvg

Polar solvation energy file obtained from g_mmpbsa.

***
    -a apolar.xvg, --apolar apolar.xvg

Non-Polar solvation energy file obtained from g_mmpbsa.

***
    -bs, --bootstrap

If given, Enable Boot Strap analysis to calculate standard error.

***
    -nbs 500, --nbstep 500

Number of boot strap steps for average energy and standard error calculation.

***
    -of full_energy.dat, --outfr full_energy.dat

Energy File: All energy components in function of time.

***
    -os summary_energy.dat, --outsum summary_energy.dat

Final Energy File: Summary of energy components.

***
    -om meta_energy.dat, --outmeta meta_energy.dat

Final Energy File for Multiple Complexes: Complex wise net binding energy.

***

[ ↑ top ↑](#top)


***

### <a name="mmpbsadecomppy"> </a> MmPbSaDecomp.py

This scripts calculate final contribution energy of each residue from individual energetic terms obtained from the g_mmpbsa

#### Usage
    python MmPbSaDecomp.py [-h] [-m contrib_MM.dat] [-p contrib_pol.dat]
                            [-a contrib_apol.dat] [-bs] [-nbs 500] [-ct 999]
                            [-o final_contrib_energy.dat] [-om energyMapIn.dat]

***

#### Description

    -h, --help

show this help message and exit

***
    -m contrib_MM.dat, --molmech contrib_MM.dat

Molecular Mechanics energy contribution file obtained from g_mmpbsa.

***
    -p contrib_pol.dat, --polar contrib_pol.dat

Polar solvation energy contribution file obtained from g_mmpbsa.

***
    -a contrib_apol.dat, --apolar contrib_apol.dat

Non-Polar solvation energy contribution file obtained from g_mmpbsa.

***
    -bs, --bootstrap

If given, Enable Boot Strap analysis to calculate standard error.

***
    -nbs 500, --nbstep 500

Number of boot strap steps for average energy and standard error calculation.

***
    -ct 999, --cutoff 999

Absolute Cutoff in kJ/mol: Energy output above and below this value. If its value is 999, all residues energy will be in output.

***
    -o final_contrib_energy.dat, --output final_contrib_energy.dat

Final output file containing binding energy contribution of each residue.

***
    -om energyMapIn.dat, --outmap energyMapIn.dat

It is input file for `energy2bfac`. It can be used with `energy2bfac` to map energy on structure for visualization.

***

[ ↑ top ↑](#top)

***

### <a name="energy2bfac"> </a> energy2bfac

This tool maps the binding energy contribution of each residue on the structure. The energy will be written in the B-factor field of the output PDB file/s. These PDB files can be used with any molecular visualizer and residues can be colored according to their energetic contribution. The molecular visualizer should support method to color residues by the B-factor values.

**WARNING:** `tpr/tpx` file should not contain PBC broken molecule. One may check by generating a PDB file with following command:`editconf -f topol.tpr -o check.pdb`. Check `check.pdb` file through visualization.

#### Usage

    Option     Filename  Type         Description
    ------------------------------------------------------------
      -s      topol.tpr  Input        Structure+mass(db): tpr tpb tpa gro g96 pdb
      -i decomp_energy.dat  Input        Generic data file
      -n      index.ndx  Input, Opt.  Index file
      -c    complex.pdb  Output, Opt. Protein data bank file
     -s1  subunit_1.pdb  Output, Opt. Protein data bank file
     -s2  subunit_2.pdb  Output, Opt. Protein data bank file
     
    Option       Type   Value   Description
    ------------------------------------------------------
    -[no]h       bool   yes     Print help info and quit
    -[no]version bool   no      Print version info and quit
    -nice        int    19      Set the nicelevel
    -b           time   0       First frame (ps) to read from trajectory
    -e           time   0       Last frame (ps) to read from trajectory
    -dt          time   0       Only use frame when t MOD dt = first time (ps)
    -tu          enum   ps      Time unit: fs, ps, ns, us, ms or s
    -[no]w       bool   no      View output .xvg, .xpm, .eps and .pdb files
    

***

#### Description

    -s topol.tpr

Input tpr/tpx file of molecule

***
    -n index.ndx

Input atomic index file. User will get choice to select atomic groups

***
    -i decomp_energy.dat

File containing energy contribution of each residue obtained from `MmPbSaDecomp.py`. One can use `-i energyMapIn.dat` directly for the input.

***
    -c complex.pdb

Output PDB file of molecular complex taken for the calculation, e.g. protein-ligand, protein-DNA or protein-protein.

***
    -s1 subunit_1.pdb

Output PDB file of first sub-unit in accordance with first atomic group chosen through index file. 

***
    -s2 subunit_2.pdb 

Output PDB file of second sub-unit in accordance with second atomic group chosen through index file. 

***
