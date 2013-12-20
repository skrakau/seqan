# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/krakau/apps/bam_sam_conv/CMakeLists.txt
#
# CMakeLists.txt file for bam_sam_conv.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_krakau_apps_bam_sam_conv)
message (STATUS "Configuring sandbox/krakau/apps/bam_sam_conv")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (bam_sam_conv bam_sam_conv.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (bam_sam_conv ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (bam_sam_conv)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install bam_sam_conv in ${PREFIX}/bin directory
install (TARGETS bam_sam_conv
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/bam_sam_conv for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (bam_sam_conv)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:bam_sam_conv")
  set (CPACK_PACKAGE_NAME "bam_sam_conv")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "bam_sam_conv")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (bam_sam_conv "bam_sam_conv")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:bam_sam_conv")

set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} bam_sam_conv CACHE INTERNAL "")

