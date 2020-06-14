#
# This file is part of g_mmpbsa.
#
# Authors: Rashmi Kumari and Andrew Lynn
# Contribution: Rajendra Kumar
# Contribution: Ievgen Sliusar
#
# Copyright (C) 2013-2016 Rashmi Kumari and Andrew Lynn
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
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc -static-libgfortran")
endif(BUILD_STATIC)

set(CMAKE_PREFIX_PATH "" CACHE STRING "Extra locations to search for external libraries and tools (give directory without lib, bin, or include)")


##############################################################
################### CHECKING GROMACS PATH ####################
##############################################################
find_package(GROMACS COMPONENTS "libgromacs")
if (NOT GROMACS_FOUND)
	message(STATUS "Continue checking for other Gromacs version library file")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
	find_package(GROMACS COMPONENTS "libgmx")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
	message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs\n  or  \n export CMAKE_PREFIX_PATH=\"/path/to/gromacs\"\n" )
endif(NOT GROMACS_FOUND)

message(STATUS "Found Gromacs include directory: " ${GROMACS_INCLUDE_DIRS})

if(${GROMACS_VERSION} VERSION_LESS "4.5.0")
	message(FATAL_ERROR, "\n Not implemented for Gromacs version: ${GROMACS_VERSION} \n")
endif()

if(${GROMACS_VERSION} VERSION_GREATER "5.0" AND ${GROMACS_VERSION} VERSION_LESS "5.1")
	set(HAVE_GROMACS50 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS50")
endif(${GROMACS_VERSION} VERSION_GREATER "5.0" AND ${GROMACS_VERSION} VERSION_LESS "5.1")

if(${GROMACS_VERSION} VERSION_GREATER "5.1" AND ${GROMACS_VERSION} VERSION_LESS "2016")
	set(HAVE_GROMACS51 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS51")
endif(${GROMACS_VERSION} VERSION_GREATER "5.1" AND ${GROMACS_VERSION} VERSION_LESS "2016")

if(NOT ${GROMACS_VERSION} VERSION_LESS "2016")
	set(HAVE_GROMACS2016 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS2016")
endif(NOT ${GROMACS_VERSION} VERSION_LESS "2016")

if(${GROMACS_VERSION} VERSION_GREATER "5.0")
	ENABLE_LANGUAGE(CXX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GROMACS_CXX_STANDARD}")
endif(${GROMACS_VERSION} VERSION_GREATER "5.0")

include_directories( ${GROMACS_INCLUDE_DIRS} )
add_definitions( ${GROMACS_DEFINITIONS} )
##############################################################


if (NOT EXT_APBS AND HAVE_GROMACS2016)
	message(FATAL_ERROR "\nOnly external APBS invocation is supported for Gromacs-${GROMACS_VERSION}\n")
endif (NOT EXT_APBS AND HAVE_GROMACS2016)

##############################################################
######################## APBS-1.3 PATH  ######################
##############################################################
if (NOT EXT_APBS AND NOT APBS14)

	if(${GROMACS_VERSION} VERSION_EQUAL "5.0" OR ${GROMACS_VERSION} VERSION_GREATER "5.0")
		message(WARNING "\nTNG format trajectory support in Gromacs-${GROMACS_VERSION} conflicts the APBS-1.3 library. \n 1) Re-install Gromacs by using \"-DGMX_USE_TNG=off\" option to disable TNG trajectory support. \n\tor\n 2) Use APBS-1.4.x.\n\t IGNORE this WARNING, IF any of the above option is already considered, ")
	endif()

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
###### GROMACS-45/46: Compiling and linking Without APBS #####
##############################################################
if(EXT_APBS AND NOT HAVE_GROMACS50 AND NOT HAVE_GROMACS51 AND NOT HAVE_GROMACS2016)
	add_executable(g_mmpbsa GMX46/main.c GMX46/g_mmpbsa.h GMX46/energy_mm.c GMX46/InputPBSA.c GMX46/nsc.c GMX46/PbsaPrep.c GMX46/psize.c  GMX46/radius.c GMX46/extrn_apbs.c ExtractData.c ExtractData.h)
	target_link_libraries(g_mmpbsa ${GROMACS_LIBRARIES})
endif()
##############################################################


##############################################################
##### GROMACS-45/46: Compiling and linking With APBS-1.3 #####
##############################################################
if(NOT EXT_APBS AND NOT APBS14 AND NOT HAVE_GROMACS50 AND NOT HAVE_GROMACS51 AND NOT HAVE_GROMACS2016)
	add_executable(g_mmpbsa GMX46/main.c GMX46/g_mmpbsa.h GMX46/energy_mm.c GMX46/InputPBSA.c GMX46/nsc.c GMX46/PbsaPrep.c GMX46/psize.c GMX46/radius.c apbs13/apbs_main.c apbs13/apbs_routines.c apbs13/apbs_routines.h)
	SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE Fortran)
	target_link_libraries(g_mmpbsa ${GROMACS_LIBRARIES} ${APBS_LIBS} )
endif()
##############################################################


##############################################################
##### GROMACS-45/46: Compiling and linking With APBS-1.4 #####
##############################################################
if(NOT EXT_APBS AND APBS14 AND NOT HAVE_GROMACS50 AND NOT HAVE_GROMACS51 AND NOT HAVE_GROMACS2016)
	add_executable(g_mmpbsa GMX46/main.c GMX46/g_mmpbsa.h GMX46/energy_mm.c GMX46/InputPBSA.c GMX46/nsc.c GMX46/PbsaPrep.c GMX46/psize.c GMX46/radius.c apbs14/apbs.h apbs14/apbs_main.c apbs14/apbs_routines.c apbs14/apbs_routines.h)
	SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE CXX)
	target_link_libraries(g_mmpbsa ${GROMACS_LIBRARIES} ${APBS_LIBS} "-fopenmp")
endif()
##############################################################


##############################################################
####### GROMACS-50: Compiling and linking Without APBS #######
##############################################################
if(EXT_APBS AND (HAVE_GROMACS50 OR HAVE_GROMACS51) )
	add_executable(g_mmpbsa GMX50/main.c GMX50/g_mmpbsa.h GMX50/energy_mm.c GMX50/InputPBSA.c GMX50/nsc.c GMX50/PbsaPrep.c GMX50/psize.c  GMX50/radius.c GMX50/extrn_apbs.c ExtractData.c ExtractData.h)
	SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE CXX)
	target_link_libraries(g_mmpbsa ${GROMACS_LIBRARIES})
endif()
##############################################################


##############################################################
###### GROMACS-50 Compiling and linking With APBS-1.3 ########
##############################################################
if(NOT EXT_APBS AND NOT APBS14 AND (HAVE_GROMACS50 OR HAVE_GROMACS51) )
	add_executable(g_mmpbsa GMX50/main.c GMX50/g_mmpbsa.h GMX50/energy_mm.c GMX50/InputPBSA.c GMX50/nsc.c GMX50/PbsaPrep.c GMX50/psize.c GMX50/radius.c apbs13/apbs_main.c apbs13/apbs_routines.c apbs13/apbs_routines.h)
	SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE CXX)
	target_link_libraries(g_mmpbsa ${GROMACS_LIBRARIES} ${APBS_LIBS} "-lgfortran -lz")
endif()
##############################################################


##############################################################
####### GROMACS-50: Compiling and linking With APBS-1.4 ######
##############################################################
if(NOT EXT_APBS AND APBS14 AND (HAVE_GROMACS50 OR HAVE_GROMACS51) )
	add_executable(g_mmpbsa GMX50/main.c GMX50/g_mmpbsa.h GMX50/energy_mm.c GMX50/InputPBSA.c GMX50/nsc.c GMX50/PbsaPrep.c GMX50/psize.c GMX50/radius.c apbs14/apbs.h apbs14/apbs_main.c apbs14/apbs_routines.c apbs14/apbs_routines.h)
	SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE CXX)
	target_link_libraries(g_mmpbsa ${GROMACS_LIBRARIES} ${APBS_LIBS} "-lgfortran")
endif()
##############################################################


##############################################################
###### GROMACS-20xx: Compiling and linking Without APBS ######
##############################################################
if(EXT_APBS AND HAVE_GROMACS2016)
	add_executable(g_mmpbsa GMX2K/main.cxx GMX2K/g_mmpbsa.h GMX2K/energy_mm.cxx GMX2K/InputPBSA.cxx GMX2K/nsc.cxx GMX2K/PbsaPrep.cxx GMX2K/psize.cxx GMX2K/radius.cxx GMX2K/extrn_apbs.cxx GMX2K/ported/readinp.cxx GMX2K/ported/warninp.cxx ExtractData.c ExtractData.h)
	SET_SOURCE_FILES_PROPERTIES(ExtractData.c PROPERTIES LANGUAGE CXX)
	#SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE CXX)
	target_link_libraries(g_mmpbsa ${GROMACS_LIBRARIES})
endif()
##############################################################


##############################################################
###### GROMACS-2020 and newer: do header fixups         ######
##############################################################
if(HAVE_GROMACS2016)
	include(FixupHeader)
	set(FIXUP_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/GMX2K/fixups)
	set(CMAKE_REQUIRED_INCLUDES "${GROMACS_INCLUDE_DIRS}")
	FixupHeaderCXX("gromacs/commandline/cmdlineinit.h")
	FixupHeaderCXX("gromacs/commandline/filenm.h")
	FixupHeaderCXX("gromacs/commandline/pargs.h")
	FixupHeaderCXX("gromacs/fileio/confio.h")
	FixupHeaderCXX("gromacs/fileio/oenv.h")
	FixupHeaderCXX("gromacs/fileio/filetypes.h")
	FixupHeaderCXX("gromacs/fileio/readinp.h")
	FixupHeaderCXX("gromacs/fileio/warninp.h")
	FixupHeaderCXX("gromacs/fileio/xvgr.h")
	FixupHeaderCXX("gromacs/math/vec.h")
	FixupHeaderCXX("gromacs/pbcutil/pbc.h")
	FixupHeaderCXX("gromacs/pbcutil/rmpbc.h")
	FixupHeaderCXX("gromacs/topology/index.h")
	FixupHeaderCXX("gromacs/utility/arrayref.h")
	FixupHeaderCXX("gromacs/utility/arraysize.h")
	FixupHeaderCXX("gromacs/utility/cstringutil.h")
	FixupHeaderCXX("gromacs/utility/exceptions.h")
	FixupHeaderCXX("gromacs/utility/fatalerror.h")
	FixupHeaderCXX("gromacs/utility/futil.h")
endif()
##############################################################


################## energy2bfac ###############################
if(HAVE_GROMACS2016)
	add_executable(energy2bfac GMX2K/energy2bfac.cxx)
	target_link_libraries(energy2bfac ${GROMACS_LIBRARIES})
elseif(NOT HAVE_GROMACS50 AND NOT HAVE_GROMACS51)
	add_executable(energy2bfac GMX46/energy2bfac.c)
	target_link_libraries(energy2bfac ${GROMACS_LIBRARIES})
else()
	add_executable(energy2bfac GMX50/energy2bfac.c)
	SET_TARGET_PROPERTIES(energy2bfac PROPERTIES LINKER_LANGUAGE CXX)
	target_link_libraries(energy2bfac ${GROMACS_LIBRARIES})
endif()


install(TARGETS g_mmpbsa DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
install(TARGETS energy2bfac DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
