# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/krakau/apps/merge_VCFs/CMakeLists.txt
#
# CMakeLists.txt file for merge_VCFs.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_krakau_apps_merge_VCFs)
message (STATUS "Configuring sandbox/krakau/apps/merge_VCFs")

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
add_executable (merge_VCFs merge_VCFs.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (merge_VCFs ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (merge_VCFs)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install merge_VCFs in ${PREFIX}/bin directory
install (TARGETS merge_VCFs
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/merge_VCFs for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (merge_VCFs)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:merge_VCFs")
  set (CPACK_PACKAGE_NAME "merge_VCFs")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "merge_VCFs")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (merge_VCFs "merge_VCFs")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:merge_VCFs")



set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} merge_VCFs CACHE INTERNAL "")
