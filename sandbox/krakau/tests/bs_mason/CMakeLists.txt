# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/krakau/tests/bs_mason/CMakeLists.txt
#
# CMakeLists.txt file for the bs_mason module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_krakau_tests_bs_mason)
message (STATUS "Configuring sandbox/krakau/tests/bs_mason")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# Search BOOST, bs_mason is not built if Boost is not found.
find_package (Boost)

if (NOT Boost_FOUND OR NOT ZLIB_FOUND)
    message (STATUS "  Boost or zlib not found: Not building bs_mason.")
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
add_executable (test_bs_mason
                test_bs_mason.cpp
                test_bs_mason.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_bs_mason ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_bs_mason COMMAND $<TARGET_FILE:test_bs_mason>)
