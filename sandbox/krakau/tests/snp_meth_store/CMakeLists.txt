# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/krakau/tests/snp_meth_store/CMakeLists.txt
#
# CMakeLists.txt file for the snp_meth_store module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_krakau_tests_snp_meth_store)
message (STATUS "Configuring sandbox/krakau/tests/snp_meth_store")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# Search BOOST, snp_store is not built if Boost is not found.
find_package (Boost)

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

# Update the list of file names below if you add source files to your test.
add_executable (test_snp_meth_store
                test_snp_meth_store.cpp
                test_snp_meth_store.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_snp_meth_store ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_snp_meth_store COMMAND $<TARGET_FILE:test_snp_meth_store>)
