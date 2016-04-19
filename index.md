---
layout: page
---

### Introduction


The development of g_mmpbsa package is initiated under [Open Source Drug Discovery Consortium (OSDD)][OSDD]{:target="_blank"}, which is a collaborative platform to design and discover new drugs for neglected tropical diseases such as Malaria, Tuberculosis, Leshmaniasis etc.

g_mmpbsa is developed using two widely used open source software i.e. [GROMACS][GROMACS]{:target="_blank"} and [APBS][APBS]{:target="_blank"} and it has similar user interface like other GROMACS tools.

The tool calculates components of binding energy using MM-PBSA method except the entropic term and energetic contribution of each residue to the binding using energy decomposition scheme.

The output from the tool is used further as input in python scripts which is provided in this package, to get the final binding energy and energetic contribution of each residue.

Kindly post problems and queries in [**g_mmpbsa forum**][forum]{:target="_blank"}, we will try our best to provide the solution.


***

#### Please always cite following two publications:

* Kumari _et al_ (2014) [g_mmpbsa - A GROMACS tool for high-throughput MM-PBSA calculations][g_mmpbsa paper]{:target="_blank"}. _J. Chem. Inf. Model._ 54:1951-1962.

* Baker _et al._ (2001) [Electrostatics of nanosystems: Application to microtubules and the ribosome][apbs paper]{:target="_blank"}. _Proc. Natl. Acad. Sci. USA_  98:10037-10041.

***


### Features v1.6


*   It is an open source tool and can be modified under the terms of the GNU public license

*   It is implemented through open source software [GROMACS][GROMACS]{:target="_blank"} and [APBS][APBS]{:target="_blank"}, making it accessible to large number of users.

*   Supports **GROMACS 4.5.x**, **4.6.x**, **GROMACS 5.0.x** and **GROMACS 5.1.x** versions.

*   Supports **APBS 1.2.x**, **1.3.x** and **1.4.x** versions

*   Inherits APBS capability of parallel computation using **OpenMP**. See details [here](How-to-Run.html#openmp).

*   Supports **external APBS** execuatble with **mpirun** for parallel computation on HPC. See details [here](How-to-Run.html#mpirun).

*   Options for van der Waal radii that are used for solvation free energy calculation using implicit solvent models.

*   Options for several non-polar solvation model such as SASA, SAV and Weeks–Chandler–Andersen (WCA).

*   Options for calculating contribution of each residue in the net binding energy.

*   Simultaneous computations of binding energy components and residue wise energy contribution, and thus it is computationally less expensive.

*   This tool can be modified and/or redistributed under terms of GNU public license.


***

#### Other Citations:

* Pronk _et al._ (2013) [GROMACS 4.5: a high-throughput and highly parallel open source molecular simulation toolkit][gromacs paper]{:target="_blank"}. _Bioinformatics_ 29:845-854.

* Eisenhaber _et al._ (1995) [The double cubic lattice method: Efficient approaches to numerical integration of surface area and volume and to dot surface contouring of molecular assemblies][sasa paper]{:target="_blank"}. _J. Comput. Chem._ 16:273-284.

* Wagoner _et al._ (2006) [Assessing implicit models for nonpolar mean solvation forces: The importance of dispersion and volume terms][wca paper]{:target="_blank"}. _Proc. Natl. Acad. Sci. USA_  103:8331-8336.

* * *

[OSDD]: http://www.osdd.net/
[GROMACS]: http://www.gromacs.org/
[APBS]: http://www.poissonboltzmann.org/
[forum]: https://groups.google.com/d/forum/g_mmpbsa
[g_mmpbsa paper]: http://pubs.acs.org/doi/abs/10.1021/ci500020m
[apbs paper]: http://www.pnas.org/content/98/18/10037.abstract
[gromacs paper]: http://bioinformatics.oxfordjournals.org/content/29/7/845.abstract
[sasa paper]: http://onlinelibrary.wiley.com/doi/10.1002/jcc.540160303/abstract
[wca paper]: http://www.pnas.org/content/103/22/8331.abstract
