### Introduction

The development of g_mmpbsa package is initiated under [Open Source Drug Discovery Consortium (OSDD)][OSDD], which is a collaborative platform to design and discover new drugs for neglected tropical diseases such as Malaria, Tuberculosis, Leshmaniasis etc.

g_mmpbsa is developed using two widely used open source software i.e. [GROMACS][GROMACS] and [APBS][APBS] and it has similar user interface like other GROMACS tools.

The tool calculates components of binding energy using MM-PBSA method except the entropic term and energetic contribution of each residue to the binding using energy decomposition scheme.

The output from the tool is used further as input in python scripts which is provided in this package, to get the final binding energy and energetic contribution of each residue.

#### For complete documentation, please visit [g_mmpbsa](http://rashmikumari.github.io/g_mmpbsa/).
#### Please post problems and queries in [g_mmpbsa forum](https://groups.google.com/d/forum/g_mmpbsa).

***

#### Please always cite following two publications:

* Kumari _et al_ (2014) [g_mmpbsa - A GROMACS tool for high-throughput MM-PBSA calculations][g_mmpbsa paper]. _J. Chem. Inf. Model._ 54:1951-1962.

* Baker _et al._ (2001) [Electrostatics of nanosystems: Application to microtubules and the ribosome][apbs paper]. _Proc. Natl. Acad. Sci. USA_  98:10037-10041.

***

### Features v2.0

*   It is an open source tool and can be modified under the terms of the GNU public license

*   It is implemented through open source software [GROMACS][GROMACS] and [APBS][APBS], making it accessible to large number of users.

*   Supports **GROMACS 2021**  and all previous versions.

*   Supports any **APBS** versions - now as an **external** tool

*   Inherits APBS capability of parallel computation using **OpenMP**.

*   Options for several non-polar solvation model such as SASA and SAV.

*   Options for calculating contribution of each residue in the net binding energy.

*   Simultaneous computations of binding energy components and residue wise energy contribution, and thus it is computationally less expensive.

*   This tool can be modified and/or redistributed under terms of GNU public license.

***

### Removed features from previous versions

*   Removed internal APBS linking, now always requires an external APBS program.

*   Removed van der Waal radii options. Now by default only Bondii radii is used. Any missing radius is taken from force-field parameters.

*   Weeks–Chandler–Andersen (WCA) non-polar solvation model is removed.


***

#### Other Citations:

* Pronk _et al._ (2013) [GROMACS 4.5: a high-throughput and highly parallel open source molecular simulation toolkit][gromacs paper]. _Bioinformatics_ 29:845-854.

* Eisenhaber _et al._ (1995) [The double cubic lattice method: Efficient approaches to numerical integration of surface area and volume and to dot surface contouring of molecular assemblies][sasa paper]. _J. Comput. Chem._ 16:273-284.

* Wagoner _et al._ (2006) [Assessing implicit models for nonpolar mean solvation forces: The importance of dispersion and volume terms][wca paper]. _Proc. Natl. Acad. Sci. USA_  103:8331-8336.

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
