# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/krakau/apps/benchmark_calls/CMakeLists.txt
#
# CMakeLists.txt file for benchmark_calls.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_krakau_apps_benchmark_calls)
message (STATUS "Configuring sandbox/krakau/apps/benchmark_calls")

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
add_executable (benchmark_calls benchmark_calls.cpp 
                                benchmark_calls.h
                                benchmark_calls_vs_snpstore.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (benchmark_calls ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (benchmark_calls)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install benchmark_calls in ${PREFIX}/bin directory
install (TARGETS benchmark_calls
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/benchmark_calls for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (benchmark_calls)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:benchmark_calls")
  set (CPACK_PACKAGE_NAME "benchmark_calls")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "benchmark_calls")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Sabrina")
  set (CPACK_PACKAGE_VENDOR "Sabrina")

  seqan_configure_cpack_app (benchmark_calls "benchmark_calls")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:benchmark_calls")


set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} benchmark_calls CACHE INTERNAL "")
