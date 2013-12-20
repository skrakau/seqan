# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/example/apps/ex_app/CMakeLists.txt
#
# CMakeLists.txt file for ex_app.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_krakau_apps_snp_meth_store)
message (STATUS "Configuring sandbox/krakau/apps/snp_meth_store")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# Search BOOST, snp_store is not built if Boost is not found.
find_package (Boost REQUIRED)

if (NOT Boost_FOUND OR NOT ZLIB_FOUND)
    message (STATUS "  Boost or zlib not found: Not building snp_store.")
    return ()
endif (NOT Boost_FOUND OR NOT ZLIB_FOUND)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (${Boost_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (snp_meth_store snp_meth_store.cpp 
                               snp_meth_store.h 
                               meths.h 
                               bs_one_calling.h
                               bs_alphabets.h
                               util.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (snp_meth_store ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})


# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

MESSAGE("Boost information:")
MESSAGE(STATUS "** Search Boost root: ${BOOST_ROOT}")
MESSAGE("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}") 
MESSAGE("  Boost_LIBRARIES: ${Boost_LIBRARIES}") 
MESSAGE("  Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}") 

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (snp_meth_store)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install ex_app in ${PREFIX}/bin directory
install (TARGETS snp_meth_store
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/ex_app for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (snp_meth_store)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:snp_meth_store")
  set (CPACK_PACKAGE_NAME "snp_meth_store")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "snp_meth_store")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Sabrina")
  set (CPACK_PACKAGE_VENDOR "Sabrina")

  seqan_configure_cpack_app (snp_meth_store "smp_meth_store")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:snp_meth_store")



set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} snp_meth_store CACHE INTERNAL "")
