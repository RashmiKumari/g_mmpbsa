#
# This file is part of g_mmpbsa.
#
# Authors: Rashmi Kumari and Andrew Lynn
# Contribution: Rajendra Kumar
#
# Copyright (C) 2013, 2014, 2015 Rashmi Kumari and Andrew Lynn
#
# g_mmpbsa is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# g_mmpbsa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#


option(EXT_APBS "compile without APBS" OFF)
option(APBS_INSTALL "Path to APBS installation directory" "/opt/apbs-1.3-source")
option(APBS_SRC "Path to APBS source directory" "/opt/apbs")
option(APBS14 "Enable compiling and linking with APBS-1.4.x" OFF)
option(BUILD_STATIC "Build standalone static executables" OFF)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -std=gnu99 -fomit-frame-pointer -finline-functions -Wno-unused -funroll-all-loops")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

if (BUILD_STATIC)
	SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
	set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)
	set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
	set(CMAKE_SHARED_LIBRARY_C_FLAGS)
	set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
	set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
	set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
endif()

##############################################################
################### CHECKING GROMACS PATH ####################
##############################################################
if (NOT DEFINED GMX_PATH AND (NOT DEFINED GMX_LIB OR NOT DEFINED GMX_INCLUDE) )
message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs \nOR\n -DGMX_LIB=/path/to/gromacs/lib -DGMX_INCLUDE=/path/to/gromacs/include")
endif()

if(DEFINED GMX_PATH)
	set(GMX_INCLUDE "${GMX_PATH}/include")
	set(GMX_LIB "${GMX_PATH}/lib")
endif()

include_directories("${GMX_INCLUDE}/gromacs")
find_library(LIBGMX gmx HINTS ${GMX_LIB})

if (NOT LIBGMX)
	message(FATAL_ERROR "\n libgmx.a or libgmx.so not found in ${GMX_LIB} directory \n")
endif()

list(APPEND GMX_LIBS "-lpthread")
list(APPEND GMX_LIBS "-L${GMX_LIB}")
list(APPEND GMX_LIBS "-lgmx")
list(APPEND GMX_LIBS "-ldl")
list(APPEND GMX_LIBS "-lm")

##############################################################

##############################################################
############ FFTW3F LIBRARY for GROMACS-4.6.x ################
##############################################################
find_library(LIBFFTW fftw3f HINTS ${FFTW_LIB})
if(NOT LIBFFTW)
	message(FATAL_ERROR "\nFFTW library file libfftw3f.so or libfftw3f.a not found... \nUSE:\n -DFFTW_LIB=/path/to/fftw3/lib/")
endif()
list(APPEND GMX_LIBS ${LIBFFTW})
#############################################################


##############################################################
######################## APBS-1.3 PATH  ######################
##############################################################
if (NOT EXT_APBS AND NOT APBS14)
	ENABLE_LANGUAGE(Fortran)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

	if(BUILD_STATIC)
		SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -static")
	endif()

	# Add definition to compile internal APBS code
	ADD_DEFINITIONS(-DINT_APBS)	

	if(NOT APBS_INSTALL)
		message(FATAL_ERROR "\nPath to APBS installation directory is not provided. \nUSE:\n -DAPBS_INSTALL=/path/to/APBS/installation/directory")
	endif()
	
	# Libraries from installation directory
	find_library(APBS_LIB1 apbsmainroutines HINTS "${APBS_INSTALL}/lib")
	find_library(APBS_LIB2 apbs HINTS "${APBS_INSTALL}/lib")
	
	if (NOT APBS_LIB1)
		message(FATAL_ERROR "\n libapbsmainroutines.a not found in ${APBS_INSTALL}/lib directory \n")
	endif()
	
	if (NOT APBS_LIB2)
		message(FATAL_ERROR "\n libapbs.a not found in ${APBS_INSTALL}/lib directory \n")
	endif()

	list(APPEND APBS_LIBS "-L${APBS_INSTALL}/lib")
	list(APPEND APBS_LIBS "-lapbsmainroutines")
	list(APPEND APBS_LIBS "-lapbs")

	# Libraries from source directory
	if (NOT APBS_SRC)
		message(FATAL_ERROR "\nPath to APBS source directory is not defined....\nUSE:\n -DAPBS_SRC=/path/to/apbs-1.3-source \n")
	endif()

	find_library(APBS_BLAS apbsblas HINTS "${APBS_SRC}/contrib/lib")
	find_library(MALOC maloc HINTS "${APBS_SRC}/contrib/lib")
	find_library(APBSGEN apbsgen HINTS "${APBS_SRC}/src/generic/.libs")
	find_library(Z z HINTS "${APBS_SRC}/contrib/zlib/.libs")
	
	if (NOT APBS_BLAS)
		message(FATAL_ERROR "\n libapbsblas.a not found in \"${APBS_SRC}/contrib/lib\" directory \n")
	endif()
	if (NOT MALOC)
		message(FATAL_ERROR "\n libmaloc.a not found in \"${APBS_SRC}/contrib/lib\" directory \n")
	endif()
	if (NOT APBSGEN)
		message(FATAL_ERROR "\n libapbsgen.a not found in \"${APBS_SRC}/src/generic/.libs\" directory \n")
	endif()
	if (NOT Z)
		message(FATAL_ERROR "\n libz.a not found in \"${APBS_SRC}/contrib/zlib/.libs\" directory \n")
	endif()

	
	list(APPEND APBS_LIBS "-L${APBS_SRC}/contrib/lib")
	list(APPEND APBS_LIBS "-L${APBS_SRC}/contrib/zlib/.libs")
	list(APPEND APBS_LIBS "-L${APBS_SRC}/src/generic/.libs")
	list(APPEND APBS_LIBS "-lapbsblas")
	list(APPEND APBS_LIBS "-lmaloc")
	list(APPEND APBS_LIBS "-lapbsgen")
	list(APPEND APBS_LIBS "-lz")

	# Headers
	include_directories("${APBS_SRC}/src/aaa_inc")
	include_directories("${APBS_SRC}/contrib/include")
endif()
################### APBS-1.3 SECTION END #####################



##############################################################
######################## APBS-1.4 PATH  ######################
##############################################################
if (NOT EXT_APBS AND APBS14)
	
	# Add definition to compile internal APBS code
	ADD_DEFINITIONS(-DINT_APBS)	

	if(NOT APBS_INSTALL)
		message(FATAL_ERROR "\nPath to APBS installation directory is not provided. \nUSE:\n -DAPBS_INSTALL=/path/to/APBS/installation/directory")
	endif()
	
	# Libraries from installation directory
	find_library(MALOC maloc HINTS "${APBS_INSTALL}/lib")
	find_library(GEOFLOW apbs_geoflow HINTS "${APBS_INSTALL}/lib")
	find_library(MG apbs_mg HINTS "${APBS_INSTALL}/lib")
	find_library(PMGC apbs_pmgc HINTS "${APBS_INSTALL}/lib")
	find_library(GENERIC apbs_generic HINTS "${APBS_INSTALL}/lib")
	
	if (NOT MALOC)
		message(FATAL_ERROR "\n libmaloc.a or libmaloc.so not found in ${APBS_INSTALL}/lib directory \n")
	endif()
	
	if (NOT GEOFLOW)
		message(FATAL_ERROR "\n libapbs_geoflow.a  or libapbs_geoflow.so not found in ${APBS_INSTALL}/lib directory \n")
	endif()
	
	if (NOT MG)
		message(FATAL_ERROR "\n libapbs_mg.a  or libapbs_mg.so not found in ${APBS_INSTALL}/lib directory \n")
	endif()
	
	if (NOT PMGC)
		message(FATAL_ERROR "\n libapbs_pmgc.a  or libapbs_pmgc.so not found in ${APBS_INSTALL}/lib directory \n")
	endif()
	
	if (NOT GENERIC)
		message(FATAL_ERROR "\n libapbs_generic.a  or libapbs_generic.so not found in ${APBS_INSTALL}/lib directory \n")
	endif()

	list(APPEND APBS_LIBS "-L${APBS_INSTALL}/lib")
	list(APPEND APBS_LIBS "-lmaloc")
	list(APPEND APBS_LIBS "-lapbs_geoflow")
	list(APPEND APBS_LIBS "-lapbs_mg")
	list(APPEND APBS_LIBS "-lapbs_pmgc")
	list(APPEND APBS_LIBS "-lapbs_generic")

	# Headers
	include_directories("${APBS_INSTALL}/include")
	include_directories("${APBS_INSTALL}/include/apbs")
endif()

################### APBS-1.4 SECTION END #####################





##############################################################
########### Compiling and linking Without APBS ###############
##############################################################
if(EXT_APBS)
	add_executable(g_mmpbsa main.c g_mmpbsa.h energy_mm.c InputPBSA.c nsc.c PbsaPrep.c psize.c  radius.c extrn_apbs.c ExtractData.c ExtractData.h)
	target_link_libraries(g_mmpbsa ${LIBPTHREAD} ${GMX_LIBS} ${LIBDL} ${LIBM})
endif()
##############################################################


##############################################################
########### Compiling and linking With APBS-1.3 ##############
##############################################################
if(NOT EXT_APBS AND NOT APBS14)
	add_executable(g_mmpbsa main.c g_mmpbsa.h energy_mm.c InputPBSA.c nsc.c PbsaPrep.c psize.c radius.c apbs13_main.c apbs13_routines.c apbs13_routines.h)
	SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE Fortran)
	target_link_libraries(g_mmpbsa ${LIBPTHREAD} ${GMX_LIBS} ${LIBDL} ${APBS_LIBS} ${LIBM})
endif()
##############################################################


##############################################################
########### Compiling and linking With APBS-1.4 ##############
##############################################################
if(NOT EXT_APBS AND APBS14)
	add_executable(g_mmpbsa main.c g_mmpbsa.h energy_mm.c InputPBSA.c nsc.c PbsaPrep.c psize.c radius.c apbs14.h apbs14_main.c apbs14_routines.c apbs14_routines.h)
	SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE CXX)
	target_link_libraries(g_mmpbsa ${LIBPTHREAD} ${GMX_LIBS} ${LIBDL} ${APBS_LIBS} ${LIBM})
endif()
##############################################################



################# energy2bfac ###############################
add_executable(energy2bfac energy2bfac.c)
target_link_libraries(energy2bfac ${LIBPTHREAD} ${GMX_LIBS} ${LIBDL} ${LIBM})


install(TARGETS g_mmpbsa DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
install(TARGETS energy2bfac DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
