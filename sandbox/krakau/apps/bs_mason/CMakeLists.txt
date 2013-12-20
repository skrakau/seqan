# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/krakau/apps/bs_mason/CMakeLists.txt
#
# CMakeLists.txt file for bs_mason.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_krakau_apps_bs_mason)
message (STATUS "Configuring sandbox/krakau/apps/bs_mason")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (bs_mason bs_mason.cpp 
                         bs_mason.h
                         simulate_454.h
                         simulate_454_base_calling.h
                         simulate_illumina.h
                         simulate_illumina_bs.h
                         simulate_sanger.h
                         store_config.h
                         util.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (bs_mason ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (bs_mason)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install bs_mason in ${PREFIX}/bin directory
install (TARGETS bs_mason
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/bs_mason for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (bs_mason)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:bs_mason")
  set (CPACK_PACKAGE_NAME "bs_mason")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "bs_mason")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Sabrina")
  set (CPACK_PACKAGE_VENDOR "Sabrina")

  seqan_configure_cpack_app (bs_mason "bs_mason")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:bs_mason")
