import sys
import os
def main():
    options = {'run': 'Run MM/PBSA calculation using trajectory and tpr file',
               'apbs' : 'apbs program used in PB calculation',
               'energy2bfac': 'Generate PDB file where decomposed binding energy is written in B-factor column',
               'average': 'Calculate final average binding energy including all energy terms. Supports multiple complexes at once.',
               'correlation'   : 'Same as average, but also calculates correlation between predicted and experimental binding energies',
               'decompose': 'Calculate final average decomposed energies of residues with plots',
              }

    program = sys.argv[0]
    os.environ['APBS'] = f'{program} apbs' # set-up APBS env which will be later used in mmpbsa
    
    if len(sys.argv)<=1:
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] not in options:
        print(' ERROR: "{0}" is not an accepted option.\n' .format(sys.argv[1]))
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] == 'run':
        from .g_mmpbsa import mmpbsa
        mmpbsa([program + ' run'] + sys.argv[2:])
        
    if sys.argv[1] == 'apbs':
        from .g_mmpbsa import apbs
        apbs(['apbs'] + sys.argv[2:])
    
    if sys.argv[1] == 'energy2bfac':
        from .g_mmpbsa import energy2bfac
        energy2bfac([program + ' energy2bfac'] + sys.argv[2:])
        
    if sys.argv[1] == 'average':
        from . import MmPbSaStat
        sys.argv = sys.argv[1:]
        MmPbSaStat.main()
        
    if sys.argv[1] == 'correlation':
        from . import MmPbSaStat_correlation
        sys.argv = sys.argv[1:]
        MmPbSaStat_correlation.main()
        
    if sys.argv[1] == 'decompose':
        from . import MmPbSaDecomp
        sys.argv = sys.argv[1:]
        MmPbSaDecomp.main()
        
def show_help(options):
    print(' ==============================')
    print(' Usage:')
    print(' g_mmpbsa <Option>\n')
    print(' ---------------------')
    print(' Use following options:')
    print(' ---------------------\n')

    for tool in options:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print(' ==============================')


if __name__=="__main__":
    main()
