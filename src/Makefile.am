#############################################################
# Input to automake
# To generate Makefile.in
# Contains information about compile and linking the program
############################################################

CCLD = $(F77) $(FFLAGS)
bin_PROGRAMS = g_mmpbsa energy2bfac
g_mmpbsa_SOURCES= apbs_main.c  energy_mm.c  g_mmpbsa.h  InputPBSA.c  main.c  nsc.c  PbsaPrep.c  psize.c  radius.c  routines.c  routines.h
energy2bfac_SOURCES= energy2bfac.c
INCLUDES = @gmx_include@ @apbs_include@
AM_LDFLAGS = @lib_path_list@ 
g_mmpbsa_LDADD = @lib_list@
energy2bfac_CC= $(CC)
energy2bfac_LDADD = @lib_list@
