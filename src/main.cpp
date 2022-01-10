/**
 * This file is part of g_mmpbsa.
 *
 * Authors: Rajendra Kumar, Rashmi Kumari and Andrew Lynn
 *
 * Copyright (C) 2013-2021 Rashmi Kumari and Andrew Lynn
 * Copyright (C) 2022- Rajendra Kumar and Rashmi Kumari
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

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>

#include <gromacs/trajectoryanalysis.h>
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include <gromacs/utility/exceptions.h>

using namespace gmx;

#define FLAG_DOTS       01
#define FLAG_VOLUME     02
#define FLAG_ATOM_AREA  04
extern int nsc_dclm ( rvec *coords, real *radius, int nat, int  densit,
                      int mode, real *value_of_area, double **at_area,
                      real *value_of_vol, double ** at_vol, real **lidots,
                      int *nu_dots, int index[] );

/**
* @brief Data structure to hold information about non-bonded interactions
* @author Rashmi Kumari
*/
typedef struct	{
    int nr_nb; /**< Number of non-bonded pairs excluding 1-2 and 1-3 pairs */
    int nr_14; /**<  Number of 1-4 pairs */
    int **pairNB;  /**< 2D array of of size `(nr_nb, 2)` containing all atom index of non-bonded pairs */
    int **pair14;  /**< 2D array of size `(nr_14, 2)` containing all 1-4 pairs atom index */
    int *pairtype; /**< array of pair-types for 1-4 Lenard-Jones parameters */
    bool *bItsA14; /**< To classify that first Atom of pairs is from subunit A in 1-4 pair list */
    bool *bItsA;   /**< To classify that first Atom of pairs is from subunit A in non-bonded pair list */
    bool *bItsB14; /**< To classify that second Atom of pairs is from subunit B in 1-4 pair list */
    bool *bItsB;   /**< To classify that second Atom of pairs is from subunit B in non-bonded pair list */

} t_non_bonded;

/**
* @brief Data structure to hold information about polar solvation energy calculations
* @author Rashmi Kumari
*/
typedef struct {
    real cfac;            /**< Factor by which molecular dimensions should expand to get a coarse grid dimensions */
    real gridspace;       /**< Fine grid spacing in Angstrom */
    real gmemceil;        /**<  Maximum memory (MB) available*/
    real fadd;            /**<  The amount (in Å) to add to molecular dimensions to get a fine grid dimensions */
    real ofrac;           /**<  Used in mg-para: Overlap in mesh during parallel calculations */
    bool bParallel;   /**<  Whether parallel calculation is required */
    bool bFocus;      /**<  Whether focus type is enabled */

    int mg_type;          /**<  multi-grid calculation type */
    ivec dime;            /**<  grid dimenstions */
    ivec pdime;           /**<  grid dimenstions for parallel calculations */
    rvec cglen;           /**<  coarse grid lengths */
    rvec cgcent;          /**<  center of coarse-grid box */
    rvec fglen;           /**<  fine grid lengths */
    rvec fgcent;          /**<  center of fine grid box */
    int pbsolver;         /**<  PB Solver type: linear or non-linear */
    int bcfl;             /**<  Type of boundary conditions */
    real pcharge;         /**<  Magnitude of positive charge */
    real ncharge;         /**<  Magnitude of negative charge */
    real prad;            /**<  Radius of positive charged ion */
    real nrad;            /**<  Radius of negative charged ion */
    real pconc;           /**<  Concentration of positive charge */
    real nconc;           /**<  Concentration of negative charge */
    real pdie;            /**<  Solute dielectric constant */
    real sdie;            /**<  Solvent dielectric constant */
    real vdie;            /**<  Vacuum or reference dielectric constant*/
    int srfm;             /**<  To construct the dielectric and ion-accessibility coefficients */
    int chgm;             /**<  To map the biomolecular point charges to the grid */
    real sdens;           /**<  Number of grid points per Å^2 for constructing the molecular surface or solvent accessible surface */
    real srad;            /**<  Radius (in Å) of solvent molecules */
    real swin;            /**<  Value for cubic spline window for spline-based surface definitions */
    real temp;            /**<  Temperature in K */

    // Apolar calculation
    real   press;       /**<  Solvent pressure proportionality term of SAV model */
    real   savconst;    /**<  Offset or constant of SAV model */
    real   savrad;      /**<  Solvent radius for SAV */
    real   gamma;       /**<  Surface tension proportionality term of SASA model */
    real   sasaconst;   /**<  Offset or constant of SASA mode */
    real   sasrad;      /**<  Solvent radius for SASA */

} t_pbsa_inputs;


/**
* @brief index for radius values in 2D array
* @ingroup PBSA_PREP
* @author Rashmi Kumari
*/
enum { eSASRAD, /**< For SASA model */
       eSAVRAD, /**< For SAV model */
     };

/**
* @brief index for PB solver
* @ingroup PBSA_INPUT
* @author Rashmi Kumari
*/
enum {eLPBE, /**< Linear solver */ eNPBE  /**< Non-linear solver */};
static const char *PBsolver[] = { "lpbe", "npbe", NULL };

/**
* @brief index of keywords to boundary condition keywords
* @ingroup PBSA_INPUT
* @author Rashmi Kumari
*/
enum {ezero, esdh, emdh};
static const char *bcfl_words[] = { "zero", "sdh", "mdh", NULL };

/**
* @brief index of keywords to method for mapping the biomolecular point charges to the grid
* @ingroup PBSA_INPUT
* @author Rashmi Kumari
*/
enum { espl2, espl4, emol, esmol};
static const char *srfm_words[] = { "spl2", "spl4",  "mol", "smol", NULL};

enum {espl0 = 2};
static const char *chgm_words[] = {"spl2", "spl4", "spl0", NULL};

enum { sacc };
static const char *APsrfm_words[] = { "sacc", NULL};

enum { mg_auto, mg_para };
static const char *mg_words[] = { "mg-auto", "mg-para", NULL};

typedef struct {
    const char *key;
    const char *author;
    const char *title;
    const char *journal;
    int volume, year;
    const char *pages;
} t_citerec;

void show_citation ( FILE *fp, const char *key )
{

    static const t_citerec citedb[] = {
        {
            "Abraham2015",
            "M. J. Abraham, T. Murtola, R. Schulz, S. Páll, J. C. Smith, B. Hess, E. Lindahl",
            "GROMACS: High performance molecular simulations through multi-level parallelism from "
            "laptops to supercomputers",
            "SoftwareX", 1, 2015, "19-25"
        },
        {
            "APBS2001",
            "N A Baker, D Sept, S Joseph, M J Holst, and J A McCammon",
            "Electrostatics of nanosystems: Application to microtubules and the ribosome",
            "Proc. Natl. Acad. Sci. USA",
            98, 2001, "10037-10041"
        },
        {
            "APBS2006",
            "J A Wagoner and N A Baker",
            "Assessing implicit models for nonpolar mean solvation forces: The importance of dispersion and volume terms",
            "Proc. Natl. Acad. Sci. USA",
            103, 2006, "8331-8336"
        },
        {
            "Eisenhaber95",
            "F Eisenhaber and P Lijnzaad and P Argos and C Sander and M Scharf",
            "The Double Cube Lattice Method: Efficient Approaches to Numerical Integration of Surface Area and Volume and to Dot Surface Contouring of Molecular Assemblies",
            "J. Comp. Chem.",
            16, 1995, "273-284"
        },
        {
            "gmmpbsa2014",
            "R Kumari, R Kumar, OSDD, and A Lynn",
            "g_mmpbsa - A GROMACS tool for high-throughput MM-PBSA calculations",
            "J. Chem. Inf. Model.",
            54, 2014, "1951–1962"
        },
    };

#define NSTR static_cast<int>(asize(citedb))

    int   index;
    char* author;
    char* title;
#define LINE_WIDTH 79

    if ( fp == nullptr ) {
        return;
    }

    for ( index = 0; index < NSTR && ( strcmp ( citedb[index].key, key ) != 0 ); index++ ) {}

    fprintf ( fp, "\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n" );
    if ( index < NSTR ) {
        /* Insert newlines */
        author = wrap_lines ( citedb[index].author, LINE_WIDTH, 0, FALSE );
        title  = wrap_lines ( citedb[index].title, LINE_WIDTH, 0, FALSE );
        fprintf ( fp, "%s\n%s\n%s %d (%d) pp. %s\n", author, title, citedb[index].journal,
                  citedb[index].volume, citedb[index].year, citedb[index].pages );
        sfree ( author );
        sfree ( title );
    } else {
        fprintf ( fp, "Entry %s not found in citation database\n", key );
    }
    fprintf ( fp, "-------- -------- --- Thank You --- -------- --------\n\n" );
    fflush ( fp );
}

/**
 * Copyright message to dispaly at the start of execution
 */
void CopyRightMsg()
{
    std::string copyright =
        "                                                                        \n"
        "                        :-)  g_mmpbsa (-:                               \n"
        "                                                                        \n"
        "                                                                        \n"
        "       Copyright (C) 2013 - 2021 Rashmi Kumari and Andrew Lynn          \n"
        "       Copyright (C) 2022- Rajendra Kumar and Rashmi Kumari             \n"
        "                                                                        \n"
        "g_mmpbsa is free software: you can redistribute it and/or modify        \n"
        "it under the terms of the GNU General Public License as published by    \n"
        "the Free Software Foundation, either version 3 of the License, or       \n"
        "(at your option) any later version.                                     \n"
        "                                                                        \n"
        "g_mmpbsa is distributed in the hope that it will be useful,             \n"
        "but WITHOUT ANY WARRANTY; without even the implied warranty of          \n"
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           \n"
        "GNU General Public License for more details.                            \n"
        "                                                                        \n"
        "You should have received a copy of the GNU General Public License       \n"
        "along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.       \n"
        "                                                                        \n"
        "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     \n"
        "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     \n"
        "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   \n"
        "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    \n"
        "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   \n"
        "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED\n"
        "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  \n"
        "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  \n"
        "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    \n"
        "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      \n"
        "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            \n"
        "                                                                        \n"
        "                                                                        \n"
        "                                                                        \n";
    fprintf ( stderr,"%s\n", copyright.c_str() );
}

/*! \brief
 * MM-PBSA analysis class
 */
class AnalysisMMPBSA : public TrajectoryAnalysisModule
{
public:
    AnalysisMMPBSA();

    virtual void initOptions ( IOptionsContainer  *options,  TrajectoryAnalysisSettings *settings );
    virtual void optionsFinished ( TrajectoryAnalysisSettings * settings );
    virtual void initAnalysis ( const TrajectoryAnalysisSettings &settings, const TopologyInformation &top );
    virtual void analyzeFrame ( int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata );
    virtual void finishAnalysis ( int nframes );
    virtual void writeOutput();

private:
    class ModuleData;

    // Input options
    real                             rvdw_;
    real                             pdie_;
    int                              ndots_;
    bool                             bDIFF_;
    bool                             bMM_;
    bool                             bDCOMP_;
    bool                             bPBSA_;
    bool                             bFocus_;
    bool                             bIncl14_;
    bool                             bVerbose_;
    std::string                      fnMDP_;
    std::string                      fnDist_;
    std::string                      fnVacMM_;
    std::string                      fnDecompMM_;
    std::string                      fnDecompPol_;
    std::string                      fnDecompAPol_;
    std::string                      fnPolar_;
    std::string                      fnAPolar_;

    // atom selections
    Selection                        selA_;
    Selection                        selB_;
    int                              numSelGroups_ = 1; // number of selection
    int                              *isize_;   // size of the selections
    const int                        **index_; // atom indices in selections

    // Varaible related to topology
    const gmx_mtop_t                 *mtop_;
    const gmx_localtop_t             *localtop_;
    AtomsDataPtr                      atoms_;
    
    // flags showing residue belonging to subunit A or B 
    std::vector<bool>                 bResA_ = { false }; // wheather a residue in subunit A
    std::vector<bool>                 bResB_ = { false }; // wheather a residue in subunit B
    
    // Vacuum MM energy variables
    t_non_bonded                      paramNonBond_;
    std::vector<real>                 EEnergyFrame_ = { 0 };
    std::vector<real>                 VdwEnergyFrame_ = { 0 };
    FILE                              *fDecompMM_;
    AnalysisData                      mmEnergyData_;
    
    // PBSA common variables
    t_pbsa_inputs                     pbsaInputKwords_;
    bool                              bPolar_ = false;
    bool                              bApolar_ = false;
    bool                              bSASA_ = false;
    bool                              bSAV_ = false;
    
    // Apolar variables
    real                              *radii_[2];          // radius for SASA/SAV calculation
    real                              totalArea_[3];    // total area
    real                              totalVolume_[3];  // volume
    double                            *atomArea_[3];    // area atom-wise
    double                            *atomVolume_[3];  // volume atom-wise
    real                              *apolarAtomsEnergy_[3]; // Apolar energy atom-wise
    FILE                              *fDecompAPol_;
    AnalysisData                      apolarEnergyData_;
    
    // Polar solvation energy variables
    std::string                       fnPQR_;
    std::string                       fnPolAPBS_;
    std::string                       fnApbsOut_;
    FILE                              *fDecompPol_;
    std::string                       apbsCommand_;
    real                              polarEnergyFrame_[3];    // polar energy of each group
    real                              *atomsPolarEnergy_[3];    // polar energy atom-wise
    AnalysisData                      polarEnergyData_;

    // Common functions
    void prepareOutputFiles ( const TrajectoryAnalysisSettings *settings );
    void writeOutputFrame ( int frnr, real time, TrajectoryAnalysisModuleData *pdata );

    // vacuum MM energy functions
    void buildNonBondedPairList();
    void readPBSAInputs();
    void vaccumMMFull ( rvec *x );
    void vaccumMMWithoutExclusions ( rvec *x );
    
    // PBSA common functions
    void assignRadius();
    std::vector<real> decomposeSolvationEnergy ( real **atomsEnergy );
    
    // polar solvation energy functions
    void psize ( rvec *x, int group );
    void createPolarInputForAPBS();
    void makePQR ( rvec *x, int group );
    void executeAPBS ( int group );
};




void
AnalysisMMPBSA::initOptions ( IOptionsContainer          *options,
                              TrajectoryAnalysisSettings *settings )
{
    static const char *const desc[] = {
        "g_mmpbsa calculates relative binding free energy using the MM-PBSA method for "
        "bio-molecular associations such as protein-protein, protein-ligand, protein-DNA etc. "
        "It calculates three components of the binding energy in separate files, so that",
        "user will have choice to calculate MM, PB and SA energy values according to their ",
        "objective. It also calculates contribution of each residue to the net binding energy",
        "and provides information about important contributing residues to the molecular",
        "association.\n\n",
        "For more detail, see please visit <http://rashmikumari.github.io/g_mmpbsa>"
    };

    settings->setHelpText ( desc );

    options->addOption ( FileNameOption ( "i" )
                         .legacyType ( efMDP ).inputFile()
                         .store ( &fnMDP_ ).defaultBasename ( "pbsa" )
                         .description ( "Input parameters for mmpbsa calculations." ) );

    options->addOption ( FileNameOption ( "mm" )
                         .filetype ( eftPlot ).outputFile()
                         .store ( &fnVacMM_ ).defaultBasename ( "energy_MM" )
                         .description ( "Vaccum MM energy as a function of time" ) );

    options->addOption ( FileNameOption ( "pol" )
                         .filetype ( eftPlot ).outputFile()
                         .store ( &fnPolar_ ).defaultBasename ( "polar" )
                         .description ( "Polar solvation energy as a function of time" ) );

    options->addOption ( FileNameOption ( "apol" )
                         .filetype ( eftPlot ).outputFile()
                         .store ( &fnAPolar_ ).defaultBasename ( "apolar" )
                         .description ( "Apolar solvation energy as a function of time" ) );

    options->addOption ( FileNameOption ( "mmcon" )
                         .filetype ( eftGenericData ).outputFile()
                         .store ( &fnDecompMM_ ).defaultBasename ( "contrib_MM" )
                         .description ( "Vacuum MM energy contribution to binding" ) );

    options->addOption ( FileNameOption ( "pcon" )
                         .filetype ( eftGenericData ).outputFile()
                         .store ( &fnDecompPol_ ).defaultBasename ( "contrib_pol" )
                         .description ( "Polar solvation energy contribution to binding" ) );


    options->addOption ( FileNameOption ( "apcon" )
                         .filetype ( eftGenericData ).outputFile()
                         .store ( &fnDecompAPol_ ).defaultBasename ( "contrib_apol" )
                         .description ( "Apolar solvation energy contribution to binding" ) );

    options->addOption ( FileNameOption ( "o" )
                         .filetype ( eftPlot ).outputFile()
                         .store ( &fnDist_ ).defaultBasename ( "avedist" )
                         .description ( "Average distances from reference group" ) );

    options->addOption ( SelectionOption ( "unit1" )
                         .store ( &selA_ ).required()
                         .description ( "Select protein or first group" ) );

    options->addOption ( SelectionOption ( "unit2" )
                         .store ( &selB_ )
                         .description ( "Select ligand or second group" ) );

    options->addOption ( RealOption ( "rvdw" ).store ( &rvdw_ ).defaultValue ( 0.1 )
                         .description ( "Default van der Waal radius (in nm) if radius not found for any atom-types)" ) );

    options->addOption ( RealOption ( "pdie" ).store ( &pdie_ ).defaultValue ( 1 )
                         .description ( "Dielectric constant of solute. Should be same as of polar solvation" ) );

    options->addOption ( IntegerOption ( "ndots" ).store ( &ndots_ ).defaultValue ( 24 )
                         .description ( "Number of dots per sphere in the calculation of SASA, more dots means more accuracy" ) );

    options->addOption ( BooleanOption ( "mme" ).store ( &bMM_ ).defaultValue ( true )
                         .description ( "To calculate vacuum molecular mechanics energy" ) );

    options->addOption ( BooleanOption ( "pbsa" ).store ( &bPBSA_ ).defaultValue ( false )
                         .description ( "To calculate polar and/or non-polar solvation energy" ) );

    options->addOption ( BooleanOption ( "diff" ).store ( &bDIFF_ ).defaultValue ( true )
                         .description ( "Calculate the energy difference between two group otherwise only calculates for one group" ) );

    options->addOption ( BooleanOption ( "decomp" ).store ( &bDCOMP_ ).defaultValue ( false )
                         .description ( "Number of dots per sphere in the calculation of SASA, more dots means more accuracy" ) );

    options->addOption ( BooleanOption ( "incl_14" ).store ( &bIncl14_ ).defaultValue ( false )
                         .description ( "Include 1-4 atom-pairs, exclude 1-2 and 1-3 atom pairs during MM calculation. Should be \"yes\" when groups are bonded with each other." ) );

    options->addOption ( BooleanOption ( "focus" ).store ( &bFocus_ ).defaultValue ( false )
                         .description ( "To enable focusing on the specfic region of molecule, group of atoms must be provided in index file" ) );

    options->addOption ( BooleanOption ( "silent" ).store ( &bVerbose_ ).defaultValue ( false )
                         .description ( "Display messages, output and errors from external APBS program" ) );

    settings->setFlag ( TrajectoryAnalysisSettings::efRequireTop );
    settings->setFlag ( settings->efNoUserPBC,true );
    settings->setFlag ( settings->efNoUserRmPBC, true );
}

void AnalysisMMPBSA::optionsFinished ( TrajectoryAnalysisSettings * )
{

    // Sanity check for input options
    
    if ( ( !bDIFF_ ) && ( bDCOMP_ ) ) {
        printf ( "\n\nWARNING: For single group calculations, decompositon cannot be used, switching it off!   \n\n" );
        bDCOMP_ = false;
    }

    if ( ( fnDecompMM_.empty() ) && ( bDCOMP_ ) && ( bMM_ ) ) {
        GMX_THROW ( InconsistentInputError ( "Decompositon requested, however. -mmcon option is missing. Aborting!!!\n" ) );
    }

    if ( ( !bPBSA_ ) && ( !bMM_ ) )   {
        GMX_THROW ( InconsistentInputError ( "No calculation opted. Use either \"-pbsa\" or \"-mm\" option.\n" ) );
    }

    if ( ( bMM_ ) && ( !bDIFF_ ) && ( !bIncl14_ ) ) {
        printf ( "\n\nWARNING: For single group calculations, 1-4 interactions are also included. \n\n" );
        bIncl14_ = TRUE;

    }

    if ( ( bPBSA_ ) && ( fnMDP_.empty() ) ) {
        GMX_THROW ( InconsistentInputError ( "Input parameter file for the PBSA calculation is missing, Use \"-i\" option\n" ) );
    }

    if ( bPBSA_ ) {
        readPBSAInputs();

        if ( ( bPolar_ ) && ( fnPolar_.empty() ) )  {
            GMX_THROW ( InconsistentInputError ( "Polar output file is missing. Use \"-pol\" option\n" ) );
        }

        if ( ( bPolar_ ) && ( fnDecompPol_.empty() ) && ( bDCOMP_ ) )  {
            GMX_THROW ( InconsistentInputError ( "Decompositon requested, however. -pcon option is missing. Aborting!!!\n" ) );
        }

        if ( ( bApolar_ ) && ( pbsaInputKwords_.gamma != 0 ) )
            bSASA_ = true;

        if ( ( bApolar_ ) && ( pbsaInputKwords_.press != 0 ) )
            bSAV_ = true;

        if ( ( bApolar_ ) && ( fnAPolar_.empty() ) )  {
            GMX_THROW ( InconsistentInputError ( "Apolar output file is missing. Use \"-apol\" option\n" ) );
        }

        if ( ( bApolar_ ) && ( fnDecompAPol_.empty() ) && ( bDCOMP_ ) )  {
            GMX_THROW ( InconsistentInputError ( "Decompositon requested, however. -apcon option is missing. Aborting!!!\n" ) );
        }
    }

}


void
AnalysisMMPBSA::initAnalysis ( const TrajectoryAnalysisSettings &settings,
                               const TopologyInformation         &topInfo )
{

    // more sanity checks for input options
    if ( !selA_.isValid() ) {
        GMX_THROW ( InconsistentInputError ( "No atoms are selected in first group!" ) );
    }

    if ( ( !selB_.isValid() ) && ( bDIFF_ ) ) {
        GMX_THROW ( InconsistentInputError ( "-diff is requested but no atoms are selected in second group!" ) );
    }

    // Read topology informations
    mtop_ = topInfo.mtop();
    localtop_ = topInfo.expandedTopology();
    atoms_ = topInfo.copyAtoms();

    // Building complex atom indices
    if ( bDIFF_ ) {
        snew ( isize_, 3 );
        snew ( index_, 3 );

        isize_[0] = selA_.atomCount();
        isize_[1] = selB_.atomCount();
        isize_[2] = isize_[0] + isize_[1];

        index_[0] = selA_.atomIndices().data();
        index_[1] = selB_.atomIndices().data();

        int *tmpIndex;
        snew ( tmpIndex, isize_[0] + isize_[1] );
        for ( int i = 0; i < isize_[0]; i++ )
            tmpIndex[i] = index_[0][i];
        for ( int i = 0; i < isize_[1]; i++ )
            tmpIndex[i + isize_[0]] = index_[1][i];

        index_[2] = tmpIndex;

        numSelGroups_ = 3;
    } else {
        snew ( isize_, 1 );
        snew ( index_, 1 );
    }

    // building bResA_ and bResB_ array here
    if ( bDCOMP_ ) {
        int nres = 0 ; //nres = total number of residue
        std::vector<int> resnmr ( atoms_->nres ); //resnmr = residue number
        std::vector<int> ResIstart ( atoms_->nres ); // Atom indeices at which each residue start
        int prev_res, curr_res;

        bResA_.resize ( atoms_->nres );
        bResB_.resize ( atoms_->nres );
        std::fill ( bResA_.begin(), bResA_.end(), false );
        std::fill ( bResB_.begin(), bResB_.end(), false );

        resnmr[0] = atoms_->resinfo[atoms_->atom[index_[2][0]].resind].nr;
        nres = 1;
        ResIstart[0] = index_[2][0];
        prev_res = atoms_->atom[index_[2][0]].resind;
        bResA_[atoms_->atom[index_[2][0]].resind] = TRUE;
        for ( int i=0; i<isize_[2]; i++ ) {
            curr_res = atoms_->atom[index_[2][i]].resind;
            if ( curr_res != prev_res ) {
                nres++;
                resnmr[nres-1] = atoms_->resinfo[curr_res].nr;
                ResIstart[nres-1] = index_[2][i];

                if ( i<isize_[0] )
                    bResA_[atoms_->atom[index_[2][i]].resind] = TRUE;
                else
                    bResB_[atoms_->atom[index_[2][i]].resind] = TRUE;

                prev_res = curr_res;
            }
        }
    }


    // Prepare output data and files
    prepareOutputFiles ( &settings );

    // build non-bonded pairs if requested
    if ( bIncl14_ ) buildNonBondedPairList();

    // prepare for PBSA calculations
    if ( bPBSA_ ) {

        assignRadius(); // assign radius to all atoms

        if ((bApolar_) && ( bDCOMP_ )) {
            for ( int i = 0; i < 3; i++ ) { // energy atom-wise memory allocation
                snew ( apolarAtomsEnergy_[i], isize_[i] );
            }
        }

        if ( bPolar_ )    {
            //Creating name for temporary data file
            char buf[256];
            strcpy ( buf, "pbXXXXX" );
            gmx_tmpnam ( buf );
            fnPQR_ = ( std::string ) buf + "A.pqr";
            fnPolAPBS_ = ( std::string ) buf + "A.in";
            fnApbsOut_ = ( std::string ) buf + "A.out";
            remove ( buf );

            // Getting path from $APBS environment
            char *apbs_env = NULL;
            apbs_env = std::getenv ( "APBS" );
            if ( apbs_env == NULL ) {
                printf ( "\n WARNING: APBS Environment variable is not defined... Looking in PATH...\n" );
                if ( system ( "which apbs > /dev/null 2>&1" ) )
                    GMX_THROW ( InconsistentInputError ( "APBS not found! Define APBS environment variable or install apbs!!" ) );
                else
                    printf ( "\n ... Found apbs in $PATH... Continuing...\n\n" );
                apbsCommand_ = "apbs";
            } else {
                apbsCommand_ = apbs_env;
                if ( !gmx_fexist ( apbsCommand_ ) )
                    GMX_THROW ( InconsistentInputError ( "APBS not found with the defined APBS environment variable.!!!" ) );
            }

            // Building apbs command
            apbsCommand_ += " " + fnPolAPBS_ + " --output-file=" + fnApbsOut_;
            if ( bVerbose_ )
                apbsCommand_ += " >/dev/null 2>&1";

            if ( bDCOMP_ )    {
                for ( int i =0; i < numSelGroups_; i++ )
                    snew ( atomsPolarEnergy_[i], isize_[i] );
            }

        }
    }
}


void
AnalysisMMPBSA::analyzeFrame ( int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata )
{
    // Electrostatic energy calculation
    if ( bMM_ ) {
        if ( bIncl14_ ) {
            vaccumMMFull ( fr.x );
        } else {
            vaccumMMWithoutExclusions ( fr.x );
        }
    }

    // Apolar calculation
    if ( bApolar_ ) {
        real tmparea; // it is required and discarded
        for ( int i = 0; i < numSelGroups_; i++ ) {

            // Calculate area
            nsc_dclm ( fr.x, radii_[eSASRAD], isize_[i], ndots_, FLAG_ATOM_AREA,
                       &totalArea_[i],  &atomArea_[i], NULL, NULL, NULL, NULL,
                       ( int * ) index_[i] );

            // Calculate volume
            nsc_dclm ( fr.x, radii_[eSAVRAD], isize_[i],  ndots_, FLAG_VOLUME,
                       &tmparea, NULL, &totalVolume_[i], &atomVolume_[i], NULL, NULL,
                       ( int * ) index_[i] );
        }

    }

    // polar calculation
    if ( bPolar_ ) {
        for ( int i = 0; i < numSelGroups_; i++ ) {
            psize ( fr.x, i );
            createPolarInputForAPBS();
            makePQR ( fr.x, i );
            executeAPBS ( i );
        }
    }

    // write output data for each frame
    writeOutputFrame ( frnr, fr.time, pdata );

}

void AnalysisMMPBSA::writeOutputFrame ( int frnr, real time, gmx::TrajectoryAnalysisModuleData* pdata )
{
    // Vacumm MM energy ouput
    if ( bMM_ ) {
        AnalysisDataHandle   mmh         = pdata->dataHandle ( mmEnergyData_ );

        // write to -mm energy_MM.xvg
        mmh.startFrame ( frnr, time );
        if ( bDIFF_ ) {
            for ( int i=0; i<3; i++ )  {
                mmh.setPoint ( ( i*2 ), VdwEnergyFrame_[i] );
                mmh.setPoint ( ( i*2 )+1, EEnergyFrame_[i] );
            }
            mmh.setPoint ( 6, ( VdwEnergyFrame_[2] + EEnergyFrame_[2] ) - ( EEnergyFrame_[0]+VdwEnergyFrame_[0]+EEnergyFrame_[1]+VdwEnergyFrame_[1] ) );
        } else {
            mmh.setPoint ( 0, VdwEnergyFrame_[0] );
            mmh.setPoint ( 1, EEnergyFrame_[0] );
            mmh.setPoint ( 2, EEnergyFrame_[0] + VdwEnergyFrame_[0] );
        }
        mmh.finishFrame();


        // write to -mmcon file
        if ( bDCOMP_ ) {
            fprintf ( fDecompMM_, "%15.3lf", time );

            for ( int i = 0; i < atoms_->nres; i++ )
                if ( ( bResA_[i] ) || ( bResB_[i] ) )
                    fprintf ( fDecompMM_, "%15.3lf", ( EEnergyFrame_[i+3] + VdwEnergyFrame_[i+3] ) /2 );
            fprintf ( fDecompMM_, "\n" );
            fflush ( fDecompMM_ );
        }
    }

    // Apolar solvation energy output
    if ( bApolar_ ) {
        AnalysisDataHandle  aphand         = pdata->dataHandle ( apolarEnergyData_ );

        // write to -apol apolar.xvg file after calculating energy from area/volume
        aphand.startFrame ( frnr, time );
        if ( bDIFF_ ) {
            if ( bSASA_ ) {
                aphand.setPoint ( 0, ( totalArea_[0] * pbsaInputKwords_.gamma * 100 ) + pbsaInputKwords_.sasaconst );
                aphand.setPoint ( 1, ( totalArea_[1] * pbsaInputKwords_.gamma * 100 ) + pbsaInputKwords_.sasaconst );
                aphand.setPoint ( 2, ( totalArea_[2] * pbsaInputKwords_.gamma * 100 ) + pbsaInputKwords_.sasaconst );
            } else {
                aphand.setPoint ( 0, 0.0 );
                aphand.setPoint ( 1, 0.0 );
                aphand.setPoint ( 2, 0.0 );
            }

            if ( bSAV_ ) {
                aphand.setPoint ( 3, ( totalVolume_[0] * pbsaInputKwords_.press * 1000 ) + pbsaInputKwords_.savconst );
                aphand.setPoint ( 4, ( totalVolume_[1] * pbsaInputKwords_.press * 1000 ) + pbsaInputKwords_.savconst );
                aphand.setPoint ( 5, ( totalVolume_[2] * pbsaInputKwords_.press * 1000 ) + pbsaInputKwords_.savconst );
            } else {
                aphand.setPoint ( 3, 0.0 );
                aphand.setPoint ( 4, 0.0 );
                aphand.setPoint ( 5, 0.0 );
            }

        } else {
            real energy = 0;
            if ( bSASA_ )  {
                energy = ( totalArea_[0] * pbsaInputKwords_.gamma * 100 ) + pbsaInputKwords_.sasaconst;
                aphand.setPoint ( 0, ( totalArea_[0] * pbsaInputKwords_.gamma * 100 ) + pbsaInputKwords_.sasaconst );
            } else {
                aphand.setPoint ( 0, 0.0 );
            }

            if ( bSAV_ ) {
                energy += ( totalVolume_[0] * pbsaInputKwords_.press * 1000 ) + pbsaInputKwords_.savconst;
                aphand.setPoint ( 1, ( totalVolume_[0] * pbsaInputKwords_.press * 1000 ) + pbsaInputKwords_.savconst );
            } else {
                aphand.setPoint ( 1, 0.0 );
            }

            aphand.setPoint ( 2, energy );
        }
        aphand.finishFrame();

        // write to -apcon file
        if ( bDCOMP_ ) {
            for ( int i = 0; i < 3; i++ ) { // calculate atom-wise energies
                for ( int j = 0; j < isize_[i]; j++ ) {
                    apolarAtomsEnergy_[i][j] = 0.0;

                    if ( bSASA_ )
                        apolarAtomsEnergy_[i][j] += ( atomArea_[i][j] * pbsaInputKwords_.gamma * 100 ) + pbsaInputKwords_.sasaconst;

                    if ( bSAV_ )
                        apolarAtomsEnergy_[i][j] += ( atomVolume_[i][j] * pbsaInputKwords_.press * 1000 ) + pbsaInputKwords_.savconst;
                }
            }

            // calculate residue-wise energies
            std::vector<real> resEnergy =  decomposeSolvationEnergy ( apolarAtomsEnergy_ );

            // write residue-wise energies to file
            fprintf ( fDecompAPol_, "%15.3lf", time );
            for ( int i=0; i < atoms_->nres; i++ ) {
                if ( ( bResA_[i] ) || ( bResB_[i] ) )
                    fprintf ( fDecompAPol_, "%15.3lf", resEnergy[i] );
            }
            fprintf ( fDecompAPol_, "\n" );
            fflush ( fDecompAPol_ );
        }

        // free memory of atomArea_ atomVolume_, these will be again allocated in next frame
        for ( int i = 0; i < numSelGroups_; i++ ) {
            sfree ( atomArea_[i] );
            sfree ( atomVolume_[i] );
        }
    }

    if ( bPolar_ ) {
        AnalysisDataHandle  phand         = pdata->dataHandle ( polarEnergyData_ );
        
        // write to -pol polar.xvg file
        phand.startFrame ( frnr, time );
        for ( int i = 0; i < numSelGroups_; i++ )
            phand.setPoint ( i, polarEnergyFrame_[i] );
        phand.finishFrame();

        // write to -pcon file
        if ( bDCOMP_ ) {
            // calculate residue-wise energies
            std::vector<real> resEnergy =  decomposeSolvationEnergy ( atomsPolarEnergy_ );

            // write residue-wise energies to file
            fprintf ( fDecompPol_, "%15.3lf", time );
            for ( int i=0; i < atoms_->nres; i++ ) {
                if ( ( bResA_[i] ) || ( bResB_[i] ) )
                    fprintf ( fDecompPol_, "%15.3lf", resEnergy[i] );
            }
            fprintf ( fDecompPol_, "\n" );
            fflush ( fDecompPol_ );
        }
    }
}


void AnalysisMMPBSA::writeOutput()
{
}


void
AnalysisMMPBSA::finishAnalysis ( int /*nframes*/ )
{
    show_citation ( stdout, "Abraham2015" );
    if ( bPolar_ )	{
        show_citation ( stdout, "APBS2001" );
    }

    if ( bApolar_ )	{
        if ( ( pbsaInputKwords_.gamma != 0 ) || ( pbsaInputKwords_.press != 0 ) )	{
            show_citation ( stdout, "Eisenhaber95" );
        }
    }
    show_citation ( stdout, "gmmpbsa2014" );
}

void AnalysisMMPBSA::buildNonBondedPairList()
{

    int i, j=0, k,l, ii, jj;
    const int *index;
    bool *globalIndex;
    int isize=0, splitIndex=0;
    gmx_bool bExclude=FALSE, bLJ14=FALSE;

    paramNonBond_.nr_nb = 1;
    paramNonBond_.nr_14 = 1;
    paramNonBond_.pairNB = ( int** ) malloc ( sizeof ( int* ) );
    paramNonBond_.pair14 = ( int** ) malloc ( sizeof ( int* ) );
    paramNonBond_.pairtype = ( int* ) malloc ( sizeof ( int ) );

    if ( bDIFF_ ) {
        snew ( paramNonBond_.bItsA,1 );
        snew ( paramNonBond_.bItsB,1 );
        snew ( paramNonBond_.bItsA14,1 );
        snew ( paramNonBond_.bItsB14,1 );

        index = index_[2];
        isize = isize_[2];
        splitIndex = isize_[0];
    } else {
        index = index_[0];
        isize = isize_[0];
    }

    // prepare global atom indices for easy access to parameters
    snew ( globalIndex, atoms_->nr );
    for ( i = 0; i < atoms_->nr; i ++ )
        globalIndex[i] = false;

    for ( i = 0; i < isize; i++ )
        globalIndex[index[i]] = true;

    // list 1-4 atom pairs
    k = 0;
    while ( k < localtop_->idef.il[F_LJ14].size() ) {
        ii = localtop_->idef.il[F_LJ14].iatoms[k+1];
        jj = localtop_->idef.il[F_LJ14].iatoms[k+2];

        if ( ( globalIndex[ii] ) && ( globalIndex[jj] ) ) {
            paramNonBond_.pairtype = ( int* ) realloc ( paramNonBond_.pairtype, ( paramNonBond_.nr_14*sizeof ( int ) ) );
            paramNonBond_.pairtype[paramNonBond_.nr_14-1] = k;

            paramNonBond_.pair14 = ( int** ) realloc ( paramNonBond_.pair14, ( paramNonBond_.nr_14*sizeof ( int* ) ) );
            paramNonBond_.pair14[paramNonBond_.nr_14-1] = ( int* ) malloc ( 2*sizeof ( int ) );
            paramNonBond_.pair14[paramNonBond_.nr_14-1][0] = ii;
            paramNonBond_.pair14[paramNonBond_.nr_14-1][1] = jj;
            bLJ14=TRUE;

            if ( bDIFF_ )	{
                srenew ( paramNonBond_.bItsA14, paramNonBond_.nr_14 );
                srenew ( paramNonBond_.bItsB14, paramNonBond_.nr_14 );
                paramNonBond_.bItsA14[paramNonBond_.nr_14-1] = FALSE;
                paramNonBond_.bItsB14[paramNonBond_.nr_14-1] = FALSE;

                if ( ( ii < splitIndex ) && ( jj  < splitIndex ) ) {
                    paramNonBond_.bItsA14[paramNonBond_.nr_14-1] = TRUE;
                    paramNonBond_.bItsB14[paramNonBond_.nr_14-1] = FALSE;
                }

                if ( ( ii >= splitIndex ) && ( jj >= splitIndex ) ) {
                    paramNonBond_.bItsB14[paramNonBond_.nr_14-1] = TRUE;
                    paramNonBond_.bItsA14[paramNonBond_.nr_14-1] = FALSE;
                }
            }
            paramNonBond_.nr_14++;
        }
        k= k+3;
    }


    // list all atom-pairs and excludes 1-2 and 1-3 atom-pairs
    float progress=0;
    for ( i = 0; i < isize; i++ ) {
        for ( j=0; j < i; j++ )	{
            bExclude=FALSE;

            for ( l = 0; l < localtop_->excls[index[i]].size(); l++ ) {
                if ( index[j] == localtop_->excls[index[i]][l] ) {
                    bExclude=TRUE;
                    break;
                }
            }

            if ( !bExclude )	{
                paramNonBond_.pairNB = ( int** ) realloc ( paramNonBond_.pairNB, ( paramNonBond_.nr_nb*sizeof ( int* ) ) );
                paramNonBond_.pairNB[paramNonBond_.nr_nb-1] = ( int* ) malloc ( 2*sizeof ( int ) );
                paramNonBond_.pairNB[paramNonBond_.nr_nb-1][0] = index[i];
                paramNonBond_.pairNB[paramNonBond_.nr_nb-1][1] = index[j];

                if ( bDIFF_ )	{
                    srenew ( paramNonBond_.bItsA,paramNonBond_.nr_nb );
                    srenew ( paramNonBond_.bItsB,paramNonBond_.nr_nb );
                    paramNonBond_.bItsA[paramNonBond_.nr_nb-1] = FALSE;
                    paramNonBond_.bItsB[paramNonBond_.nr_nb-1] = FALSE;
                    if ( ( i <  splitIndex ) && ( j  <  splitIndex ) ) {
                        paramNonBond_.bItsA[paramNonBond_.nr_nb-1] = TRUE;
                        paramNonBond_.bItsB[paramNonBond_.nr_nb-1] = FALSE;
                    }
                    if ( ( i  >= splitIndex ) && ( j  >= splitIndex ) ) {
                        paramNonBond_.bItsB[paramNonBond_.nr_nb-1] = TRUE;
                        paramNonBond_.bItsA[paramNonBond_.nr_nb-1] = FALSE;
                    }
                }
                paramNonBond_.nr_nb++;
            }
        }
        progress = ( ( float ) i/ ( float ) isize ) * 100;
        fprintf ( stderr,"\r %5.0f %% completed...",progress );
        fflush ( stdout );
    }
    printf ( "\n Finished pair generation....\nTotal %d 1-4 pairs and %d non-bonded pairs generated.\n\n",paramNonBond_.nr_14-1,paramNonBond_.nr_nb-1 );

    paramNonBond_.nr_14 = paramNonBond_.nr_14-1;
    paramNonBond_.nr_nb = paramNonBond_.nr_nb-1;
}

void AnalysisMMPBSA::vaccumMMFull ( rvec *x )
{
    rvec dx;
    real colmb_factor = 138.935485;
    double qi, qj, c6, c12, c6j, c12j,rij, c6ij, c12ij;
    int itypeA, itypeB, ntype = mtop_->atomtypes.nr;
    int i, j,k,l,n;
    int atomA, atomB, resA, resB;
    real TempEE, TempVdw;
    int nres = atoms_->nres;

    std::vector<real> EERes ( nres, 0.0 ), VdwRes ( nres, 0.0 );

    if ( bDIFF_ ) {
        if ( EEnergyFrame_.size() != nres+3 )
            EEnergyFrame_.resize ( nres+3 );
        if ( VdwEnergyFrame_.size() != nres+3 )
            VdwEnergyFrame_.resize ( nres+3 );
    }

    std::fill ( EEnergyFrame_.begin(), EEnergyFrame_.end(), 0.0 );
    std::fill ( VdwEnergyFrame_.begin(), VdwEnergyFrame_.end(), 0.0 );

    // Energy of all previously listed atom-pairs except 1-2, 1-3 and 1-4 pairs
    for ( i=0; i<paramNonBond_.nr_nb; i++ )	{
        atomA = paramNonBond_.pairNB[i][0];
        atomB = paramNonBond_.pairNB[i][1];
        resA = atoms_->atom[atomA].resind;
        resB = atoms_->atom[atomB].resind;

        rij = sqrt ( pow ( ( x[atomA][0]-x[atomB][0] ),2 )+pow ( ( x[atomA][1]-x[atomB][1] ),2 ) + pow ( ( x[atomA][2]-x[atomB][2] ),2 ) );

        itypeA = atoms_->atom[atomA].type;
        itypeB = atoms_->atom[atomB].type;

        if ( itypeA<=itypeB )	{
            c6 = localtop_->idef.iparams[itypeA*ntype+itypeB].lj.c6;
            c12 = localtop_->idef.iparams[itypeA*ntype+itypeB].lj.c12;
        } else	{
            c6 = localtop_->idef.iparams[itypeB*ntype+itypeA].lj.c6;
            c12 = localtop_->idef.iparams[itypeB*ntype+itypeA].lj.c12;
        }
        TempEE = ( ( colmb_factor/pdie_ ) * ( atoms_->atom[atomA].q*atoms_->atom[atomB].q ) ) /rij;
        TempVdw = ( c12/pow ( rij,12 ) ) - ( c6/pow ( rij,6 ) );

        if ( bDIFF_ )	{
            if ( paramNonBond_.bItsA[i] ) {
                EEnergyFrame_[0] += TempEE;
                VdwEnergyFrame_[0] += TempVdw;
            }
            if ( paramNonBond_.bItsB[i] ) {
                EEnergyFrame_[1] += TempEE;
                VdwEnergyFrame_[1] += TempVdw;
            }
            EEnergyFrame_[2] += TempEE;
            VdwEnergyFrame_[2] += TempVdw;
        } else {
            EEnergyFrame_[0] += TempEE;
            VdwEnergyFrame_[0] += TempVdw;
        }

        if ( bDCOMP_ )
            if ( ( !paramNonBond_.bItsA[i] ) && ( !paramNonBond_.bItsB[i] ) ) {
                EERes[resA] += TempEE;
                EERes[resB] += TempEE;
                VdwRes[resA] += TempVdw;
                VdwRes[resB] += TempVdw;
            }
    }
    
    // Energy of previously listed 1-4 atom pairs
    for ( i=0; i<paramNonBond_.nr_14; i++ ) {
        atomA = paramNonBond_.pair14[i][0];
        atomB = paramNonBond_.pair14[i][1];
        resA = atoms_->atom[atomA].resind;
        resB = atoms_->atom[atomB].resind;

        rij = sqrt ( pow ( ( x[atomA][0]-x[atomB][0] ),2 )+pow ( ( x[atomA][1]-x[atomB][1] ),2 ) + pow ( ( x[atomA][2]-x[atomB][2] ),2 ) );

        c6 = localtop_->idef.iparams[localtop_->idef.il[F_LJ14].iatoms[paramNonBond_.pairtype[i]]].lj14.c6A;
        c12 = localtop_->idef.iparams[localtop_->idef.il[F_LJ14].iatoms[paramNonBond_.pairtype[i]]].lj14.c12A;

        TempEE = mtop_->ffparams.fudgeQQ * ( ( colmb_factor/pdie_ ) * ( atoms_->atom[atomA].q*atoms_->atom[atomB].q ) ) /rij;
        TempVdw = ( c12/pow ( rij,12 ) ) - ( c6/pow ( rij,6 ) );

        if ( bDIFF_ ) {
            if ( paramNonBond_.bItsA14[i] ) {
                EEnergyFrame_[0] += TempEE;
                VdwEnergyFrame_[0] += TempVdw;
            }
            if ( paramNonBond_.bItsB14[i] ) {
                EEnergyFrame_[1] += TempEE;
                VdwEnergyFrame_[1] += TempVdw;
            }
            EEnergyFrame_[2] += TempEE;
            VdwEnergyFrame_[2] += TempVdw;
        } else {
            EEnergyFrame_[0] += TempEE;
            VdwEnergyFrame_[0] += TempVdw;
        }
        if ( bDCOMP_ )
            if ( ( !paramNonBond_.bItsA14[i] ) && ( !paramNonBond_.bItsB14[i] ) ) {
                EERes[resA] += TempEE;
                EERes[resB] += TempEE;
                VdwRes[resA] += TempVdw;
                VdwRes[resB] += TempVdw;
            }
    }

    if ( bDCOMP_ )
        for ( i=0; i<nres; i++ ) {
            EEnergyFrame_[i+3] = EERes[i];
            VdwEnergyFrame_[i+3] = VdwRes[i];
        }

}

void AnalysisMMPBSA::vaccumMMWithoutExclusions ( rvec *x )
{
    int i, j=0;
    real colmb_factor = 138.935485;
    double c6, c12, rij;
    int itypeA, itypeB, ntype = mtop_->atomtypes.nr;
    int resA, resB;
    real TempEE, TempVdw;
    int nres = atoms_->nres;

    std::vector<real> EERes ( nres, 0.0 ), VdwRes ( nres, 0.0 );

    if ( EEnergyFrame_.size() != nres+3 )
        EEnergyFrame_.resize ( nres+3 );
    if ( VdwEnergyFrame_.size() != nres+3 )
        VdwEnergyFrame_.resize ( nres+3 );

    std::fill ( EEnergyFrame_.begin(), EEnergyFrame_.end(), 0.0 );
    std::fill ( VdwEnergyFrame_.begin(), VdwEnergyFrame_.end(), 0.0 );

    for ( i=0; i<isize_[0]; i++ ) {
        for ( j=0; j<isize_[1]; j++ )	{

            resA = atoms_->atom[index_[0][i]].resind;
            resB = atoms_->atom[index_[1][j]].resind;

            rij = sqrt ( pow ( ( x[index_[0][i]][0]-x[index_[1][j]][0] ),2 )+pow ( ( x[index_[0][i]][1]-x[index_[1][j]][1] ),2 ) + pow ( ( x[index_[0][i]][2]-x[index_[1][j]][2] ),2 ) );

            itypeA = atoms_->atom[index_[0][i]].type;
            itypeB = atoms_->atom[index_[1][j]].type;

            if ( itypeA<=itypeB )	{
                c6 = localtop_->idef.iparams[itypeA*ntype+itypeB].lj.c6;
                c12 = localtop_->idef.iparams[itypeA*ntype+itypeB].lj.c12;
            } else	{
                c6 = localtop_->idef.iparams[itypeB*ntype+itypeA].lj.c6;
                c12 = localtop_->idef.iparams[itypeB*ntype+itypeA].lj.c12;
            }
            TempEE = ( ( colmb_factor/pdie_ ) * ( atoms_->atom[index_[0][i]].q * atoms_->atom[index_[1][j]].q ) ) /rij;
            TempVdw = ( c12/pow ( rij,12 ) ) - ( c6/pow ( rij,6 ) );

            EEnergyFrame_[2] += TempEE;
            VdwEnergyFrame_[2] += TempVdw;

            if ( bDCOMP_ ) {
                EERes[resA] += TempEE;
                EERes[resB] += TempEE;
                VdwRes[resA] += TempVdw;
                VdwRes[resB] += TempVdw;
            }
        }
    }

    if ( bDCOMP_ )
        for ( i=0; i<nres; i++ ) {
            EEnergyFrame_[i+3] = EERes[i];
            VdwEnergyFrame_[i+3] = VdwRes[i];
        }
}


void AnalysisMMPBSA::prepareOutputFiles ( const TrajectoryAnalysisSettings *settings )
{
    // Vacuum MM energy
    if ( bMM_ )  {
        int mmColumnCount = ( bDIFF_ ) ? 7 : 3;
        std::vector<std::string> mmLegends;

        // plot legends
        if ( bDIFF_ ) {
            mmLegends.push_back ( std::string ( selA_.name() ) +  " VdW Energy" );
            mmLegends.push_back ( std::string ( selA_.name() ) +  " Elec. Energy" );
            mmLegends.push_back ( std::string ( selB_.name() ) +  " VdW Energy" );
            mmLegends.push_back ( std::string ( selB_.name() ) +  " Elec. Energy" );
            if ( bIncl14_ ) {
                mmLegends.push_back ( std::string ( selA_.name() ) + "+" + std::string ( selB_.name() ) + " VdW Energy" );
                mmLegends.push_back ( std::string ( selA_.name() ) + "+" + std::string ( selB_.name() ) + " Elec. Energy" );
            } else {
                mmLegends.push_back ( std::string ( selA_.name() ) + "-" + std::string ( selB_.name() ) + " VdW Energy" );
                mmLegends.push_back ( std::string ( selA_.name() ) + "-" + std::string ( selB_.name() ) + " Elec. Energy" );
            }
            mmLegends.push_back ( std::string ( selA_.name() ) + "-" + std::string ( selB_.name() ) + " Total Energy" );
        } else {
            mmLegends.push_back ( std::string ( selA_.name() ) +  " VdW Energy" );
            mmLegends.push_back ( std::string ( selA_.name() ) +  " Elec. Energy" );
            mmLegends.push_back ( std::string ( selA_.name() ) +  " Total Energy" );
        }

        // plot file setup
        mmEnergyData_.setColumnCount ( 0, mmColumnCount );
        {
            AnalysisDataPlotModulePointer plotm ( new AnalysisDataPlotModule ( settings->plotSettings() ) );
            plotm->setFileName ( fnVacMM_ );
            plotm->setTitle ( "Vaccum MM Energy)" );
            plotm->setXAxisIsTime();
            plotm->setYLabel ( "Energy (kJ/mol)" );
            for ( size_t i = 0; i < mmColumnCount; ++i ) {
                plotm->appendLegend ( mmLegends[i] );
            }
            mmEnergyData_.addModule ( plotm );
        }

        // decompositon file setup
        if ( bDCOMP_ ) {
            fDecompMM_ = gmx_ffopen ( fnDecompMM_, "w" );
            fprintf ( fDecompMM_, "# Time\t" );
            for ( int i=0; i < atoms_->nres; i++ ) {
                if ( ( bResA_[i] ) || ( bResB_[i] ) ) {
                    fprintf ( fDecompMM_,"%s-%d\t", * ( atoms_->resinfo[i].name ), atoms_->resinfo[i].nr );
                }
            }
            fprintf ( fDecompMM_,"\n" );
        }
    }

    // non-polar solvation energy
    if ( bApolar_ ) {
        int apolrColumnCount = ( bDIFF_ ) ? 6 : 3;
        std::vector<std::string> apolarLegends;

        // plot legends
        if ( bDIFF_ ) {
            apolarLegends.push_back ( std::string ( selA_.name() ) +  "-Surf-ten energy" );
            apolarLegends.push_back ( std::string ( selB_.name() ) +  "-Surf-ten energy" );
            apolarLegends.push_back ( std::string ( selA_.name() ) + "+" + std::string ( selB_.name() )+  "-Surf-ten energy" );

            apolarLegends.push_back ( std::string ( selA_.name() ) +  "-Press-Vol energy" );
            apolarLegends.push_back ( std::string ( selB_.name() ) +  "-Press-Vol energy" );
            apolarLegends.push_back ( std::string ( selA_.name() ) + "+" + std::string ( selB_.name() )+  "-Press-Vol energy" );
        } else {
            apolarLegends.push_back ( std::string ( selA_.name() ) +  "-Surf-ten energy" );
            apolarLegends.push_back ( std::string ( selA_.name() ) +  "-Press-Vol energy" );
            apolarLegends.push_back ( std::string ( selA_.name() ) +  "Total energy" );
        }

        // plot file setup
        apolarEnergyData_.setColumnCount ( 0, apolrColumnCount );
        {
            AnalysisDataPlotModulePointer plotm ( new AnalysisDataPlotModule ( settings->plotSettings() ) );
            plotm->setFileName ( fnAPolar_ );
            plotm->setTitle ( "Apolar solvation energy)" );
            plotm->setXAxisIsTime();
            plotm->setYLabel ( "Energy (kJ/mol)" );
            for ( size_t i = 0; i < apolrColumnCount; ++i ) {
                plotm->appendLegend ( apolarLegends[i] );
            }
            apolarEnergyData_.addModule ( plotm );
        }

        // decompositon file setup
        if ( bDCOMP_ ) {
            fDecompAPol_ = gmx_ffopen ( fnDecompAPol_, "w" );
            fprintf ( fDecompAPol_, "# Time\t" );
            for ( int i=0; i < atoms_->nres; i++ ) {
                if ( ( bResA_[i] ) || ( bResB_[i] ) ) {
                    fprintf ( fDecompAPol_,"%s-%d\t", * ( atoms_->resinfo[i].name ), atoms_->resinfo[i].nr );
                }
            }
            fprintf ( fDecompAPol_,"\n" );
        }
    }

    // Polar solvation energy
    if ( bPolar_ ) {
        int polrColumnCount = ( bDIFF_ ) ? 3 : 1;
        std::vector<std::string> polarLegends;

        // plot legends
        if ( bDIFF_ ) {
            polarLegends.push_back ( std::string ( selA_.name() ) +  " PB energy" );
            polarLegends.push_back ( std::string ( selB_.name() ) +  " PB energy" );
            polarLegends.push_back ( std::string ( selA_.name() ) + "=" + std::string ( selB_.name() )+  " PB energy" );
        } else {
            polarLegends.push_back ( "PB energy" );
        }

        // plot file setup
        polarEnergyData_.setColumnCount ( 0, polrColumnCount );
        {
            AnalysisDataPlotModulePointer plotm ( new AnalysisDataPlotModule ( settings->plotSettings() ) );
            plotm->setFileName ( fnPolar_ );
            plotm->setTitle ( "Polar solvation energy)" );
            plotm->setXAxisIsTime();
            plotm->setYLabel ( "Energy (kJ/mol)" );
            for ( size_t i = 0; i < polrColumnCount; ++i ) {
                plotm->appendLegend ( polarLegends[i] );
            }
            polarEnergyData_.addModule ( plotm );
        }

        // decompositon file setup
        if ( bDCOMP_ ) {
            fDecompPol_ = gmx_ffopen ( fnDecompPol_, "w" );
            fprintf ( fDecompPol_, "# Time\t" );
            for ( int i=0; i < atoms_->nres; i++ ) {
                if ( ( bResA_[i] ) || ( bResB_[i] ) ) {
                    fprintf ( fDecompPol_,"%s-%d\t", * ( atoms_->resinfo[i].name ), atoms_->resinfo[i].nr );
                }
            }
            fprintf ( fDecompPol_,"\n" );
        }
    }
}

void AnalysisMMPBSA::readPBSAInputs()
{
    warninp_t wi;
    gmx_bool bAllowWarnings=FALSE;
    int maxwarning = 99;


    // Start reading input file
    wi = init_warning ( bAllowWarnings, maxwarning );
    gmx::TextInputFile     stream ( fnMDP_ );
    std::vector<t_inpfile> inp = read_inpfile ( &stream, fnMDP_.c_str(), wi );

    //To check for polar solvation calculation
    std::string polar = get_estr ( &inp, "polar", nullptr );
    if ( polar == "yes" ) bPolar_ = true;

    std::string apolar = get_estr ( &inp, "apolar", nullptr );
    if ( apolar == "yes" ) bApolar_ = true;


    if ( bPolar_ ) {
        //Psize keywords
        pbsaInputKwords_.cfac      = get_ereal ( &inp, "cfac", 2, wi );
        pbsaInputKwords_.gridspace = get_ereal ( &inp, "gridspace", 0.5, wi );
        pbsaInputKwords_.gmemceil  = get_ereal ( &inp, "gmemceil", 400, wi );
        pbsaInputKwords_.fadd      = get_ereal ( &inp, "fadd", 20, wi );
        pbsaInputKwords_.ofrac     = get_ereal ( &inp, "ofrac", 0.1, wi );

        //Polar Keywords
        pbsaInputKwords_.mg_type   = get_eeenum ( &inp, "mg-type", mg_words, wi );
        pbsaInputKwords_.pcharge   = get_ereal ( &inp, "pcharge", 0,    wi );
        pbsaInputKwords_.ncharge   = get_ereal ( &inp, "ncharge", 0,    wi );
        pbsaInputKwords_.prad      = get_ereal ( &inp, "prad",    0,    wi );
        pbsaInputKwords_.nrad      = get_ereal ( &inp, "nrad",    0,    wi );
        pbsaInputKwords_.pconc     = get_ereal ( &inp, "pconc",   0,    wi );
        pbsaInputKwords_.nconc     = get_ereal ( &inp, "nconc",   0,    wi );
        pbsaInputKwords_.pdie      = get_ereal ( &inp, "pdie",    4,    wi );
        pbsaInputKwords_.sdie      = get_ereal ( &inp, "sdie",    78.4, wi );
        pbsaInputKwords_.vdie      = get_ereal ( &inp, "vdie",    1,    wi );
        pbsaInputKwords_.srad      = get_ereal ( &inp, "srad",    1.4,  wi ); // same for SASA and SAV
        pbsaInputKwords_.swin      = get_ereal ( &inp, "swin",    0.30, wi );
        pbsaInputKwords_.sdens     = get_ereal ( &inp, "sdens",   10,   wi );
        pbsaInputKwords_.temp      = get_ereal ( &inp, "temp",    300,  wi );
        pbsaInputKwords_.srfm      = get_eeenum ( &inp, "srfm", srfm_words, wi );
        pbsaInputKwords_.chgm      = get_eeenum ( &inp, "chgm", chgm_words, wi );
        pbsaInputKwords_.bcfl      = get_eeenum ( &inp, "bcfl", bcfl_words, wi );
        pbsaInputKwords_.pbsolver  = get_eeenum ( &inp, "PBsolver", PBsolver, wi );



    }

    if ( bApolar_ ) {
        pbsaInputKwords_.gamma     = get_ereal ( &inp, "gamma",    0.030096,  wi );
        pbsaInputKwords_.sasaconst = get_ereal ( &inp, "sasconst",    0,  wi );
        pbsaInputKwords_.sasrad = get_ereal ( &inp, "sasrad",    1.4,  wi );

        pbsaInputKwords_.press     = get_ereal ( &inp, "press",       0,  wi );
        pbsaInputKwords_.savconst  = get_ereal ( &inp, "savconst",    0,  wi );
        pbsaInputKwords_.savrad = get_ereal ( &inp, "savrad",    1.29,  wi );
    }

}

void AnalysisMMPBSA::assignRadius()
{
    std::map<std::string, real> radiusDef = {
        {"o",  1.520},
        {"s",  1.830},
        {"n",  1.550},
        {"c",  1.700},
        {"h",  1.200},
        {"p",  1.800},
        {"f",  1.470},
        {"i",  2.060},
        {"cl", 1.770}, // Chlorine atomtype
        {"br", 1.920},
        {"ca", 1.770},
        {"cb", 1.770},
        {"cc", 1.770},
        {"cn", 1.770},
        {"cr", 1.770},
        {"cv", 1.770},
        {"cw", 1.770},
        {"c*", 1.770},
        {"cd", 1.770},
        {"ha", 1.000},
        {"h4", 1.000},
        {"h5", 1.000},
        {"c0", 2.310}, // Calcium atomtype
        {"mw", 0.050}, // virtual-sites pr sigma-holes
    };

    snew ( radii_[eSASRAD], atoms_->nr );
    snew ( radii_[eSAVRAD], atoms_->nr );
    if ( atoms_->pdbinfo == nullptr ) {
        snew ( atoms_->pdbinfo, atoms_->nr );
    } else {
        srenew ( atoms_->pdbinfo, atoms_->nr );
    }


    int itype, ntype = mtop_->atomtypes.nr;
    real c6, c12, sig6, rad = -1;
    std::string atomtype, atomname, atomname2, atomtype2;
    for ( int i = 0; i < atoms_->nr; i++ ) {    // get charge and radius for all atoms
        rad = -1;

        // calculate radius from force-field parameters
        // if radius is not found later on, force-field radius will be used
        itype = atoms_->atom[i].type;
        c12 = localtop_->idef.iparams[itype * ntype + itype].lj.c12;
        c6 =  localtop_->idef.iparams[itype * ntype + itype].lj.c6;
        if ( ( c6 != 0 ) && ( c12 != 0 ) ) {
            sig6 = c12 / c6;
            rad = 0.5 * pow ( sig6, 1.0 / 6.0 );
        } else {
            rad = rvdw_;
        }

        rad *= 10; //Conversion of nano meter to angstroms

        // Try to find radius based on atomname and atomtype
        atomtype = * ( atoms_->atomtype[i] );
        std::transform ( atomtype.begin(), atomtype.end(), atomtype.begin(), [] ( unsigned char c ) {
            return std::tolower ( c );
        } );

        atomname = * ( atoms_->atomname[i] );
        std::transform ( atomname.begin(), atomname.end(), atomname.begin(), [] ( unsigned char c ) {
            return std::tolower ( c );
        } );

        // First assign on the basis of element, i.e. first charecters of atomname
        atomname2 = atomname.at ( 0 );
        if ( radiusDef.count ( atomname2 ) > 0 ) { // atomname - mostly single charecters atomname
            rad = radiusDef[atomname2];
        } else { // atomname - mostly for two charecters atomname. Cl and Ca will be assigned by atomtype
            if ( radiusDef.count ( atomname ) > 0 )
                rad = radiusDef[atomname];
        }

        // re-assign based on atom-type
        // IMPORTANT: it overrides atomname based assignment
        if ( radiusDef.count ( atomtype ) > 0 )  { // reassign based on atomtype
            rad = radiusDef[atomtype];
        } else {
            if ( atomtype.length() > 2 ) { // reassign based on first two charecters of atomtype match
                atomtype2 = atomtype[0]+atomtype[1];
                if ( radiusDef.count ( atomtype2 ) > 0 )
                    rad = radiusDef[atomtype];
            }
        }

        // Assigned charge and radius to ocuppancy and bfactor field
        atoms_->pdbinfo[i].occup = atoms_->atom[i].q;
        atoms_->pdbinfo[i].bfac = rad;

        radii_[eSASRAD][i] = ( rad + pbsaInputKwords_.sasrad ) /10;
        radii_[eSAVRAD][i] = ( rad + pbsaInputKwords_.savrad ) /10;
    }

}

std::vector<real> AnalysisMMPBSA::decomposeSolvationEnergy ( real **atomsEnergy )
{
    std::vector<real> ResOnly ( atoms_->nres, 0.0 ), ResComplex ( atoms_->nres, 0.0 ), ResEnergy ( atoms_->nres, 0.0 );

    //Residues in Complex
    for ( int i=0; i < isize_[2]; i++ )
        ResComplex[atoms_->atom[index_[2][i]].resind] += atomsEnergy[2][i];

    //Residues in A only
    for ( int i=0; i < isize_[0]; i++ )
        ResOnly[atoms_->atom[index_[0][i]].resind] += atomsEnergy[0][i];

    //Residues in B only
    for ( int i=0; i < isize_[1]; i++ )
        ResOnly[atoms_->atom[index_[1][i]].resind] += atomsEnergy[1][i];

    for ( int i=0; i < atoms_->nres; i++ )
        ResEnergy[i] = ( ResComplex[i] - ResOnly[i] );


    return ResEnergy;
}

void AnalysisMMPBSA::psize ( rvec *x, int group )
{
    int i,j;
    rvec minlen, maxlen;
    rvec olen, clen, flen, cen;
    std::vector<int> n ( 3 ), np ( 3 ), tn ( 3 ), nsmall ( 3 );
    real nsmem, gmem, zofac;
    real r, np_float;

    minlen[XX] = 9999;
    minlen[YY] = 9999;
    minlen[ZZ] = 9999;
    maxlen[XX] = -9999;
    maxlen[YY] = -9999;
    maxlen[ZZ] = -9999;

    for ( i=0; i < isize_[group]; i++ )	{
        r = atoms_->pdbinfo[index_[group][i]].bfac;
        for ( j=0; j<DIM; j++ )	{
            if ( ( x[index_[group][i]][j]*10 )-r < minlen[j] )
                minlen[j] = ( x[index_[group][i]][j]*10 ) - r;

            if ( ( x[index_[group][i]][j]*10 )+r > maxlen[j] )
                maxlen[j] = ( x[index_[group][i]][j]*10 ) + r;

        }
    }

    for ( i=0; i<DIM; i++ )		{
        olen[i] = maxlen[i] - minlen[i];
        clen[i] = pbsaInputKwords_.cfac * olen[i];
        flen[i] = pbsaInputKwords_.fadd + olen[i];
        if ( flen[i]>clen[i] )
            flen[i] = clen[i];
        cen[i] = ( maxlen[i] + minlen[i] ) /2;

        tn[i] = ( int ) flen[i]/pbsaInputKwords_.gridspace + 0.5;
        n[i] = 32* ( ( int ) ( ( tn[i] - 1 ) / 32.0 + 0.5 ) ) + 1;
        nsmall[i] = 32* ( ( int ) ( ( tn[i] - 1 ) / 32.0 + 0.5 ) ) + 1;

        if ( nsmall[i] < 33 )
            nsmall[i] = 33;
    }

    //To Check the available memory
    gmem = 200.0 * n[XX] * n[YY] * n[ZZ] / 1024 / 1024;
    while ( 1 )	{
        nsmem = 200.0 * nsmall[XX] * nsmall[YY] * nsmall[ZZ] / 1024 / 1024;
        if ( nsmem<pbsaInputKwords_.gmemceil )
            break;
        else	{
            std::vector<int>::iterator maxElm = std::max_element ( nsmall.begin(), nsmall.end() );
            i = std::distance ( nsmall.begin(), maxElm );
            nsmall[i] = 32 * ( ( nsmall[i] - 1 ) /32 - 1 ) + 1;
            if ( nsmall[i] <= 0 )	{
                GMX_THROW ( InconsistentInputError ( "You picked a memory ceiling that is too small\n" ) );
            }


        }
    }

    // Calculating pdime => np
    if ( gmem >= pbsaInputKwords_.gmemceil )	{
        zofac = 1 + 2 * pbsaInputKwords_.ofrac;
        for ( i=0; i<DIM; i++ )	{
            np_float = n[i]/ ( float ) nsmall[i];
            if ( np_float > 1 )
                np[i] = ( int ) ( zofac*n[1]/nsmall[i] + 1.0 );
        }
    }


    if ( gmem >= pbsaInputKwords_.gmemceil )
        pbsaInputKwords_.mg_type = mg_para;
    else
        pbsaInputKwords_.mg_type = mg_auto;


    copy_rvec ( clen, pbsaInputKwords_.cglen );
    copy_rvec ( cen, pbsaInputKwords_.cgcent );
    copy_rvec ( flen, pbsaInputKwords_.fglen );
    copy_rvec ( cen, pbsaInputKwords_.fgcent );

    for ( i=0; i<DIM; i++ )	{
        pbsaInputKwords_.dime[i] = nsmall[i];
    }

    if ( pbsaInputKwords_.mg_type == mg_para ) {
        for ( i=0; i<DIM; i++ )	{
            pbsaInputKwords_.pdime[i] = np[i];
        }
    }
}

void AnalysisMMPBSA::createPolarInputForAPBS()
{
    FILE *fIn;
    fIn = gmx_ffopen ( fnPolAPBS_, "w" );

    fprintf ( fIn, "read\n    mol pqr %s\nend\n", fnPQR_.c_str() );

    if ( pbsaInputKwords_.mg_type == mg_para )	{
        fprintf ( fIn, "\nelec name mol1\n    mg-para\n" );
        fprintf ( fIn, "    dime  %d %d %d\n", pbsaInputKwords_.dime[XX], pbsaInputKwords_.dime[YY], pbsaInputKwords_.dime[ZZ] );
        fprintf ( fIn, "    pdime  %d %d %d\n", pbsaInputKwords_.pdime[XX], pbsaInputKwords_.pdime[YY], pbsaInputKwords_.pdime[ZZ] );
        fprintf ( fIn, "    ofrac %g\n", pbsaInputKwords_.ofrac );
    } else	{
        fprintf ( fIn, "\nelec name mol1\n    mg-auto\n" );
        fprintf ( fIn, "    dime  %d %d %d\n", pbsaInputKwords_.dime[XX], pbsaInputKwords_.dime[YY], pbsaInputKwords_.dime[ZZ] );
    }
    fprintf ( fIn, "    cglen %6.3lf %6.3lf %6.3lf\n", pbsaInputKwords_.cglen[XX], pbsaInputKwords_.cglen[YY], pbsaInputKwords_.cglen[ZZ] );
    fprintf ( fIn, "    fglen %6.3lf %6.3lf %6.3lf\n", pbsaInputKwords_.fglen[XX], pbsaInputKwords_.fglen[YY], pbsaInputKwords_.fglen[ZZ] );
    fprintf ( fIn, "    cgcent %6.3lf %6.3lf %6.3lf\n", pbsaInputKwords_.cgcent[XX], pbsaInputKwords_.cgcent[YY], pbsaInputKwords_.cgcent[ZZ] );
    fprintf ( fIn, "    fgcent %6.3lf %6.3lf %6.3lf\n", pbsaInputKwords_.fgcent[XX], pbsaInputKwords_.fgcent[YY], pbsaInputKwords_.fgcent[ZZ] );
    fprintf ( fIn, "    mol 1\n" );
    fprintf ( fIn, "    %s\n", PBsolver[pbsaInputKwords_.pbsolver] );
    fprintf ( fIn, "    bcfl %s\n", bcfl_words[pbsaInputKwords_.bcfl] );
    fprintf ( fIn, "    ion %.1g %.3g %.3g\n", pbsaInputKwords_.pcharge, pbsaInputKwords_.pconc, pbsaInputKwords_.prad );
    fprintf ( fIn, "    ion %.1g %.3g %.3g\n", pbsaInputKwords_.ncharge, pbsaInputKwords_.nconc, pbsaInputKwords_.nrad );
    fprintf ( fIn, "    pdie %g\n", pbsaInputKwords_.pdie );
    fprintf ( fIn, "    sdie %g\n", pbsaInputKwords_.sdie );
    fprintf ( fIn, "    srfm %s\n", srfm_words[pbsaInputKwords_.srfm] );
    fprintf ( fIn, "    chgm %s\n", chgm_words[pbsaInputKwords_.chgm] );
    fprintf ( fIn, "    sdens %g\n", pbsaInputKwords_.sdens );
    fprintf ( fIn, "    srad %g\n", pbsaInputKwords_.srad );
    fprintf ( fIn, "    swin %g\n", pbsaInputKwords_.swin );
    fprintf ( fIn, "    temp %g\n", pbsaInputKwords_.temp );
    if ( bDCOMP_ )
        fprintf ( fIn, "    calcenergy comps\n" );
    else
        fprintf ( fIn, "    calcenergy total\n" );
    fprintf ( fIn, "end\n" );

    if ( pbsaInputKwords_.mg_type == mg_para )	{
        fprintf ( fIn, "\nelec name mol2\n    mg-para\n" );
        fprintf ( fIn, "    dime  %d %d %d\n", pbsaInputKwords_.dime[XX], pbsaInputKwords_.dime[YY], pbsaInputKwords_.dime[ZZ] );
        fprintf ( fIn, "    pdime  %d %d %d\n", pbsaInputKwords_.pdime[XX], pbsaInputKwords_.pdime[YY], pbsaInputKwords_.pdime[ZZ] );
        fprintf ( fIn, "    ofrac %g\n", pbsaInputKwords_.ofrac );
    } else	{
        fprintf ( fIn, "\nelec name mol2\n    mg-auto\n" );
        fprintf ( fIn, "    dime  %d %d %d\n", pbsaInputKwords_.dime[XX], pbsaInputKwords_.dime[YY], pbsaInputKwords_.dime[ZZ] );
    }

    fprintf ( fIn, "    cglen %6.3lf %6.3lf %6.3lf\n", pbsaInputKwords_.cglen[XX], pbsaInputKwords_.cglen[YY], pbsaInputKwords_.cglen[ZZ] );
    fprintf ( fIn, "    fglen %6.3lf %6.3lf %6.3lf\n", pbsaInputKwords_.fglen[XX], pbsaInputKwords_.fglen[YY], pbsaInputKwords_.fglen[ZZ] );
    fprintf ( fIn, "    cgcent %6.3lf %6.3lf %6.3lf\n", pbsaInputKwords_.cgcent[XX], pbsaInputKwords_.cgcent[YY], pbsaInputKwords_.cgcent[ZZ] );
    fprintf ( fIn, "    fgcent %6.3lf %6.3lf %6.3lf\n", pbsaInputKwords_.fgcent[XX], pbsaInputKwords_.fgcent[YY], pbsaInputKwords_.fgcent[ZZ] );
    fprintf ( fIn, "    mol 1\n" );
    fprintf ( fIn, "    %s\n", PBsolver[pbsaInputKwords_.pbsolver] );
    fprintf ( fIn, "    bcfl %s\n", bcfl_words[pbsaInputKwords_.bcfl] );
    fprintf ( fIn, "    ion %.1g %.3g %.3g\n", pbsaInputKwords_.pcharge, pbsaInputKwords_.pconc, pbsaInputKwords_.prad );
    fprintf ( fIn, "    ion %.1g %.3g %.3g\n", pbsaInputKwords_.ncharge, pbsaInputKwords_.nconc, pbsaInputKwords_.nrad );
    fprintf ( fIn, "    pdie %g\n", pbsaInputKwords_.pdie );
    fprintf ( fIn, "    sdie %g\n", pbsaInputKwords_.vdie );
    fprintf ( fIn, "    srfm %s\n", srfm_words[pbsaInputKwords_.srfm] );
    fprintf ( fIn, "    chgm %s\n", chgm_words[pbsaInputKwords_.chgm] );
    fprintf ( fIn, "    sdens %g\n", pbsaInputKwords_.sdens );
    fprintf ( fIn, "    srad %g\n", pbsaInputKwords_.srad );
    fprintf ( fIn, "    swin %g\n", pbsaInputKwords_.swin );
    fprintf ( fIn, "    temp %g\n", pbsaInputKwords_.temp );
    if ( bDCOMP_ )
        fprintf ( fIn, "    calcenergy comps\n" );
    else
        fprintf ( fIn, "    calcenergy total\n" );
    fprintf ( fIn, "end\n" );
    fprintf ( fIn, "print elecEnergy mol1 - mol2 end\n" );
    fprintf ( fIn, "quit\n" );
    gmx_ffclose ( fIn );
}

void AnalysisMMPBSA::makePQR ( rvec* x, int group )
{
    FILE *fPQR;
    int i;
    char *resname, *atomname;
    int resnmr;
    real q, r;
    static const char *pdb = "%-6s%5u  %-5.4s%4.4s  %4d    %8.3f%8.3f%8.3f";
    //set_pdb_wide_format(TRUE);
    fPQR = gmx_ffopen ( fnPQR_, "w" );
    //write_pdbfile_indexed(fPQR,NULL,&(top->atoms),x,ePBC,box,' ',-1,isize,index,NULL,TRUE);
    for ( i = 0; i < isize_[group]; i++ ) {
        atomname = * ( atoms_->atomname[index_[group][i]] );
        resname = * ( atoms_->resinfo[atoms_->atom[index_[group][i]].resind].name );
        resnmr = atoms_->resinfo[atoms_->atom[index_[group][i]].resind].nr;
        q = atoms_->pdbinfo[index_[group][i]].occup;
        r = atoms_->pdbinfo[index_[group][i]].bfac;
        //%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f";
        fprintf ( fPQR, pdb, "ATOM", i + 1, atomname, resname, resnmr, 10 * x[index_[group][i]][XX], 10 * x[index_[group][i]][YY], 10 * x[index_[group][i]][ZZ] );
        fprintf ( fPQR, "%8.3f%8.3f\n", q, r );
    }
    gmx_ffclose ( fPQR );
}

void AnalysisMMPBSA::executeAPBS ( int group )
{
    std::string apbsCommand;
    const char *apbs_env = NULL;
    FILE *fApbsOut;
    char **data=NULL;
    bool bID_A=FALSE, bID_B =FALSE, bWCA_end=FALSE;
    int nlines, i;
    int mol1_start = 0, mol2_start = 0;
    int mol1_lastID = 0, mol2_lastID = 0;

    double totEnergy1, totEnergy2;
    double *atEnergy1=NULL, *atEnergy2=NULL;
    int at_count = 0;

    /* Executing APBS command */
    if ( 0 != system ( apbsCommand_.c_str() ) )
        GMX_THROW ( InternalError ( "Failed to execute command: " + apbsCommand_ ) );

    gmx::TextInputFile inputStream ( fnApbsOut_ );
    std::string line;
    std::vector<std::string> tempVector;
    int id, atomIndex;
    while ( inputStream.readLine ( &line ) ) {
        line = gmx::stripString ( line );
        
        if ( startsWith ( line, "id 2" ) )
            id = 2;
        
        if ( startsWith ( line, "id 4" ) )
            id = 4;
        
        if ( startsWith ( line, "end" ) )
            id = 0;
        
        if ( ( id == 2 ) && startsWith ( line, "totEnergy" ) ) {
            tempVector = splitString ( line );
            polarEnergyFrame_[group] = std::stod ( tempVector[1] );
        }
        
        if ( ( id == 2 ) && startsWith ( line, "atom" ) && ( bDCOMP_ ) ) {
            tempVector = splitString ( line );
            atomIndex = std::stoi ( tempVector[1] );
            atomsPolarEnergy_[group][atomIndex] = std::stod ( tempVector[2] );
        }
        
        if ( ( id == 4 ) && startsWith ( line, "totEnergy" ) ) {
            tempVector = splitString ( line );
            polarEnergyFrame_[group] -= std::stod ( tempVector[1] );
        }
        
        if ( ( id == 4 ) && startsWith ( line, "atom" ) && ( bDCOMP_ ) ) {
            tempVector = splitString ( line );
            atomIndex = std::stoi ( tempVector[1] );
            atomsPolarEnergy_[group][atomIndex] -= std::stod ( tempVector[2] );
        }
    }
    inputStream.close();

    // remove intermediate files
    if ( gmx_fexist ( fnApbsOut_ ) )
        remove ( fnApbsOut_.c_str() );

    if ( gmx_fexist ( fnPQR_ ) )
        remove ( fnPQR_.c_str() );

    if ( gmx_fexist ( fnPolAPBS_ ) )
        remove ( fnPolAPBS_.c_str() );

    if ( gmx_fexist ( "io.mc" ) )
        remove ( "io.mc" );
}


AnalysisMMPBSA::AnalysisMMPBSA()
{
    registerAnalysisDataset ( &mmEnergyData_, "mmEnergy" );
    registerAnalysisDataset ( &apolarEnergyData_, "apolarEnergy" );
    registerAnalysisDataset ( &polarEnergyData_, "polarEnergy" );
}


/*! \brief
 * The main function for the analysis template.
 */
int
main ( int argc, char *argv[] )
{
    CopyRightMsg();
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AnalysisMMPBSA> ( argc, argv );
}
