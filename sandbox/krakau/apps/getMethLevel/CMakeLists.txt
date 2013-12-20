# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/example/apps/getMethLevel/CMakeLists.txt
#
# CMakeLists.txt file for getMethLevel.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_example_apps_getMethLevel)
message (STATUS "Configuring sandbox/example/apps/getMethLevel")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

find_package (Boost REQUIRED)

if (NOT Boost_FOUND)
    message (STATUS "  Boost or not found: Not building getmethlevel.")
    return ()
endif (NOT Boost_FOUND )

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (${Boost_INCLUDE_DIRS})


# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (getMethLevel getMethLevel.cpp
                             getMethLevel.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (getMethLevel -lboost_regex ${SEQAN_LIBRARIES} ${Boost_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (getMethLevel)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install getMethLevel in ${PREFIX}/bin directory
install (TARGETS getMethLevel
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/getMethLevel for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (getMethLevel)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:getMethLevel")
  set (CPACK_PACKAGE_NAME "getMethLevel")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "getMethLevel")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (getMethLevel "getMethLevel")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:getMethLevel")
