#
# This file is part of g_mmpbsa.
#
# Authors: Rashmi Kumari and Andrew Lynn
# Contribution: Rajendra Kumar
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

ENABLE_LANGUAGE(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
    find_package(GROMACS COMPONENTS "libgromacs_d")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
    find_package(GROMACS COMPONENTS "libgromacs_mpi")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
	message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs\n  or  \n export CMAKE_PREFIX_PATH=\"/path/to/gromacs\"\n" )
endif(NOT GROMACS_FOUND)

message(STATUS "Found Gromacs include directory: " ${GROMACS_INCLUDE_DIRS})

if(${GROMACS_VERSION} VERSION_LESS "2023")
	message(FATAL_ERROR, "\n Not implemented for Gromacs version: ${GROMACS_VERSION} \n")
endif()
##############################################################


##############################################################
################### Including GROMACS SRC ####################
##############################################################
if(GMX_SRC)
    if(EXISTS "${GMX_SRC}/src/gromacs")
        set(GROMACS_SRC "${GMX_SRC}/src/")
        message(STATUS "Found GROMACS Source directory for header files: ${GMX_SRC}/src/")
    else()
        message(FATAL_ERROR "\n ${GMX_SRC}/src/gromacs not found in GROMACS SOURCE PATH. \nUSE:\n -DGMX_SRC=/path/to/gromacs/source\n" )
    endif()
else()
    message(FATAL_ERROR "\nGROMACS SOURCE PATH is not defined.....\nUSE:\n -DGMX_SRC=/path/to/gromacs/source\n" )
endif()

FILE(GLOB INCLUDE_DIRS RELATIVE ${GROMACS_SRC} ${GROMACS_SRC}/gromacs/*/include)
FOREACH(INC_DIR ${INCLUDE_DIRS})
  include_directories(${GROMACS_SRC}${INC_DIR})
ENDFOREACH()
##############################################################


include_directories( ${GROMACS_SRC} )
include_directories( ${GROMACS_INCLUDE_DIRS} )
add_definitions( ${GROMACS_DEFINITIONS} )

add_executable(g_mmpbsa main.cpp nsc.cpp)
SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 17)
SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES COMPILE_FLAGS "-fopenmp")
if( NOT ("${GROMACS_LIBRARIES}" MATCHES "openmp"))
	SET_TARGET_PROPERTIES(g_mmpbsa PROPERTIES LINK_LIBRARIES "-fopenmp")
endif()
target_link_libraries(g_mmpbsa ${GROMACS_LIBRARIES})
	
add_executable(energy2bfac energy2bfac.cpp)
SET_TARGET_PROPERTIES(energy2bfac PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 17)
target_link_libraries(energy2bfac ${GROMACS_LIBRARIES})

install(TARGETS g_mmpbsa DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
install(TARGETS energy2bfac DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
