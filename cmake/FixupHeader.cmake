# COPYRIGHT (c) 2016 Obsidian Research Corporation.
# Licensed under BSD (MIT variant) or GPLv2. See COPYING.

# Execute a header fixup based on NOT_NEEDED for HEADER

# The buildlib includes alternate header file shims for several scenarios, if
# the build system detects a feature is present then it should call FixupHeader
# with the test as true. If false then the shim header will be installed.

# Typically the shim header will replace a missing header with stubs, or it
# will augment an existing header with include_next.
function(DoFixup not_needed header)
  cmake_parse_arguments(ARGS "NO_SHIM" "" "" ${ARGN})

  set(SRC "${FIXUP_INCLUDE}/${header}")
  if (NOT EXISTS ${SRC})
    # NO_SHIM lets cmake succeed if the header exists in the system but no
    # shim is provided, but this will always fail if the shim is needed but
    # does not exist.
    if (NOT ARGS_NO_SHIM OR NOT "${not_needed}")
      message(FATAL_ERROR "Fixup header ${header} is not present")
    endif()
  endif()

  set(DEST "${BUILD_INCLUDE}/${header}")
  if (NOT "${not_needed}")
    if(CMAKE_VERSION VERSION_LESS "2.8.12")
      get_filename_component(DIR ${DEST} PATH)
    else()
      get_filename_component(DIR ${DEST} DIRECTORY)
    endif()
    file(MAKE_DIRECTORY "${DIR}")

    message(STATUS "Installing fixup header: ${header}")
    file(COPY ${SRC} DESTINATION ${DIR})
  else()
    file(REMOVE ${DEST})
  endif()
endfunction()

function(HeaderToVar header ovar)
  string(TOUPPER "${header}" hvar)
  string(REPLACE "/" "_" hvar "${hvar}")
  string(REPLACE "." "_" hvar "${hvar}")
  set("${ovar}" "HAVE_${hvar}" PARENT_SCOPE)
endfunction()

include(CheckIncludeFile)
function(FixupHeader header)
  HeaderToVar("${header}" hvar)
  CHECK_INCLUDE_FILE("${header}" "${hvar}")
  DoFixup("${${hvar}}" "${header}")
  if ("${${hvar}}")
    add_definitions(-D${hvar})
  endif()
endfunction()

include(CheckIncludeFileCXX)
function(FixupHeaderCXX header)
  HeaderToVar("${header}" hvar)
  CHECK_INCLUDE_FILE_CXX("${header}" "${hvar}")
  DoFixup("${${hvar}}" "${header}")
  if ("${${hvar}}")
    add_definitions(-D${hvar})
  endif()
endfunction()
