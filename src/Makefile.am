#############################################################
# Input to automake
# To generate Makefile.in
# Contains information about compile and linking the program
############################################################

CCLD = $(F77) $(FFLAGS) $(LDFLAGS) ${LIBS}
bin_PROGRAMS = g_mmpbsa energy2bfac
g_mmpbsa_SOURCES= main.c apbs_main.c  energy_mm.c  g_mmpbsa.h  InputPBSA.c nsc.c  PbsaPrep.c  psize.c  radius.c  routines.c  routines.h
energy2bfac_SOURCES= energy2bfac.c
energy2bfac_CC= $(CC)
