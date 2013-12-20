// ==========================================================================
//                                  post_pro
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef SANDBOX_KRAKAU_TESTS_POST_PRO_TEST_POST_PRO_H_
#define SANDBOX_KRAKAU_TESTS_POST_PRO_TEST_POST_PRO_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

// TODO change
struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    CharString readFileName;
    CharString readFileName2;
    CharString samFileName;
    CharString refFileName;
    CharString outputFileName;
    unsigned intervalOffset;
    double minMapq;
    unsigned max3Errors;        // max. allowed errors in 3 letter alphabet (real errors)
    unsigned max4Errors;        // max. allowed errors in 4 letter alphabet (corresponding to mapper settings)
    double maxScore;
    unsigned maxBasePenalty;    // limit the penalty for a single base 

    int minScore;
    bool outputSingleMates;
 
    AppOptions() :
        intervalOffset(3),
        minMapq(0),
        max3Errors(4),
        max4Errors(2),
        maxScore(1000000),  // TODO: what would be reasonable?
        maxBasePenalty(-3),  // scaled to single penalties
        minScore(0),
        outputSingleMates(true),    // Output also read whose mate didn't map & if no match mate pair found, output mates single 
        verbosity(1)
    {}
};



// A test for strings.
SEQAN_DEFINE_TEST(test_post_pro_strings_example1)
{
    using namespace seqan;

    // Define some constant test data for comparison...
    CharString const STRING1 = "test 1";
    CharString const STRING2 = "test 2";

    // Append to a string and make equality assertion on the result.
    CharString myStr = "test ";
    append(myStr, "1");
    SEQAN_ASSERT_EQ(STRING1, myStr);

    // Demonstration of other assertions.
    SEQAN_ASSERT_GT(STRING2, myStr);
    SEQAN_ASSERT_GEQ(STRING2, myStr);
    SEQAN_ASSERT_LT(myStr, STRING2);
    SEQAN_ASSERT_LEQ(STRING2, STRING2);
}

template<typename TReadGaps, typename TContigGaps, typename TFragmentStore, typename TId, typename TOptions>
void testReAlign4SetUp1(TReadGaps &readGaps, TContigGaps &contigGaps, TFragmentStore &store, TId &id, TOptions &options)
{
    typedef typename TFragmentStore::TAlignedReadStore                          TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                             TAlignedRead;

    id = length(store.alignedReadStore);

    // Contig gaps (empty)
    String<GapAnchor<unsigned> > contigGapAnchors;
    clear(contigGapAnchors);
    Dna5String contigInf = "ACGTACGTACGTAAAACCCCAAAACCCCAAAACCCCAAAACCCC";
    assignSource(contigGaps, contigInf);

    // Read gaps (empty)
    String<GapAnchor<unsigned> >  readGapAnchors;
    String<Dna5Q> readSeq = "ACGTgCGTACGTAAAACCCCAAAACCCCAAAACCCCAAAACCCC";
    String<int> quals;
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    appendValue(quals, 10);     // mismatch position
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    SEQAN_ASSERT_EQ(length(readSeq), length(quals));

    assignQualities(readSeq, quals);
    assignSource(readGaps, readSeq);

    // Store
    appendRead(store, readSeq);
    TId readId = length(store.readSeqStore) - 1;
    TAlignedRead alignedRead = TAlignedRead(id, readId, 0, 0, length(readSeq), readGapAnchors);
    appendValue(store.alignedReadStore, alignedRead);
    
    resize(store.alignQualityStore, length(store.alignedReadStore));
    // matePairStore ???
    
    // Options

}

template<typename TReadGaps, typename TContigGaps, typename TFragmentStore, typename TId, typename TOptions>
void testReAlign4SetUp2(TReadGaps &readGaps, TContigGaps &contigGaps, TFragmentStore &store, TId &id, TOptions &options)
{
    typedef typename TFragmentStore::TAlignedReadStore                          TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                             TAlignedRead;
 
    id = length(store.alignedReadStore);

    // Contig gaps (empty)
    String<GapAnchor<unsigned> > contigGapAnchors;
    Dna5String contigInf = "ACGTACGTACGTAAAACCCCAAAACCCCAAAACCCCAAAACCCC";
    assignSource(contigGaps, contigInf);

    // Read gaps (empty)
    String<GapAnchor<unsigned> >  readGapAnchors;
    String<Dna5Q> readSeq = "ACGTgCGTACGTAAAACCCCAAAACCCCAAAACCCCAAAACCCC";
    String<int> quals;
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    appendValue(quals, 60);     // mismatch position
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    SEQAN_ASSERT_EQ(length(readSeq), length(quals));

    assignQualities(readSeq, quals);
    assignSource(readGaps, readSeq);

    // Store
    appendRead(store, readSeq);
    TId readId = length(store.readSeqStore) - 1;
    TAlignedRead alignedRead = TAlignedRead(id, readId, 0, 0, length(readSeq), readGapAnchors);
    appendValue(store.alignedReadStore, alignedRead);
    
    resize(store.alignQualityStore, length(store.alignedReadStore));
    // matePairStore ???
    
    // Options
}

template<typename TReadGaps, typename TContigGaps, typename TFragmentStore, typename TId, typename TOptions>
void testReAlign4SetUp3(TReadGaps &readGaps, TContigGaps &contigGaps, TFragmentStore &store, TId &id, TOptions &options)
{
    typedef typename TFragmentStore::TAlignedReadStore                          TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                             TAlignedRead;
 
    id = length(store.alignedReadStore);

    // Contig gaps (empty)
    String<GapAnchor<unsigned> > contigGapAnchors;
    Dna5String contigInf = "ACGTACGTACGTAAAACCCCAAAACCCCAAAACCCCAAAACCCC";
    assignSource(contigGaps, contigInf);

    // Read gaps (empty)
    String<GapAnchor<unsigned> >  readGapAnchors;
    String<Dna5Q> readSeq = "ACGTgCGTACGTATAACCCCAAAACCCCAAAACCCCAAAACCCC";
    String<int> quals;
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    appendValue(quals, 60);     // mismatch position
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40); //
    appendValue(quals, 40);
    appendValue(quals, 40);
    appendValue(quals, 40);

    SEQAN_ASSERT_EQ(length(readSeq), length(quals));

    assignQualities(readSeq, quals);
    assignSource(readGaps, readSeq);

    // Store
    appendRead(store, readSeq);
    TId readId = length(store.readSeqStore) - 1;
    TAlignedRead alignedRead = TAlignedRead(id, readId, 0, 0, length(readSeq), readGapAnchors);
    appendValue(store.alignedReadStore, alignedRead);
    
    resize(store.alignQualityStore, length(store.alignedReadStore));
    // matePairStore ???
    
    // Options
}

// Test realigning and score 
SEQAN_DEFINE_TEST(test_reAlign4)
{
    using namespace seqan;

    // Setup datastructures
    FragmentStore<MyFragmentStoreConfig> store;
    AppOptions options;
    
    Gaps<String<Dna5Q>, AnchorGaps<String<GapAnchor<unsigned> > > > readGaps;
    Gaps<String<Dna5Q>, AnchorGaps<String<GapAnchor<unsigned> > > > contigGaps;
    unsigned id;
  
    /////////////////////////////////////
    // Fill with test case
    testReAlign4SetUp1(readGaps, contigGaps, store, id, options);
    std::cerr << "4: contigGaps: " << contigGaps << std::endl;
    // Call desired function
    reAlign4(readGaps, contigGaps, store, id, options);

    std::cout << "store.alignQualityStore[0].score: " << store.alignQualityStore[0].score << std::endl;
    std::cout << std::endl;
    // Test results 
    //SEQAN_ASSERT_EQ(store.alignQualityStore[0].errors, 1);
    // Score:
    // 59 * log( (10*(1-10^(-40/10)) )
    //double score = 59 * std::log10( 10*(1.0-pow(10, -40/10)) ) + std::log10( 10*(pow(10, -10/10)/3.0) );
    //score -= 60;
    //score *= -10;
    //SEQAN_ASSERT_EQ(store.alignQualityStore[0].score, score);

    /////////////////////////////////////
    // Fill with test case
    clear(store.alignedReadStore);
    clear(store.alignQualityStore);
    clearReads(store);

    testReAlign4SetUp2(readGaps, contigGaps, store, id, options);

    // Call desired function
    reAlign4(readGaps, contigGaps, store, id, options);

    std::cout << "store.alignQualityStore[0].score: " << store.alignQualityStore[0].score << std::endl;
    // Test results 
    //SEQAN_ASSERT_EQ(store.alignQualityStore[0].errors, 1);
    // Score:
    // 59 * log( (10*(1-10^(-40/10)) )
    //score = 59 * std::log10( 10*(1.0-pow(10, -40/10)) ) + std::log10( 10*(pow(10, -60/10)/3.0) );
    //score -= 60;
    //score *= -10;
    //SEQAN_ASSERT_EQ(store.alignQualityStore[0].score, score);
    
    /////////////////////////////////////
    // Fill with test case
    clear(store.alignedReadStore);
    clear(store.alignQualityStore);
    clearReads(store);

    testReAlign4SetUp3(readGaps, contigGaps, store, id, options);

    // Call desired function
    reAlign4(readGaps, contigGaps, store, id, options);

    std::cout << "store.alignQualityStore[0].score: " << store.alignQualityStore[0].score << std::endl;
}

template<typename TFragmentStore, typename TId, typename TOptions>
void testVerifyReadSetUp1(TFragmentStore &store, TId &id, TOptions &options)
{
    // We need:
    // readSeqStore
    // alignedReadStore: readId
    // alignQualityStore: score, errors
}


// Test verifying reads and mapq
SEQAN_DEFINE_TEST(test_verifyRead)
{
    using namespace seqan;

    // Setup datastructures
    FragmentStore<MyFragmentStoreConfig> store;
    AppOptions options;
    unsigned id;
  
    /////////////////////////////////////
    // One hit only
    testVerifyReadSetUp1(store, id, options);

    //VerifiedRead verifiedRead = verifyRead(store, options);


    /////////////////////////////////////
    // Unique hit
    


    /////////////////////////////////////
    // Nonunique hit

}


template<typename TFragmentStore, typename TId, typename TOptions>
void testSetUpVM1(TFragmentStore &store, TId &id, TOptions &options)
{
  
}


// Test verifying reads and mapq
SEQAN_DEFINE_TEST(test_verifyMates)
{
    using namespace seqan;

    // Setup datastructures
    FragmentStore<MyFragmentStoreConfig> store;
    AppOptions options;
     
    /////////////////////////////////////
    // One hit only
    unsigned id;
    testSetUpVM1(store, id, options);


    //verifyMates();



    /////////////////////////////////////
    // Unique hit
    


    /////////////////////////////////////
    // Nonunique hit


    
}




#endif  // SANDBOX_KRAKAU_TESTS_POST_PRO_TEST_POST_PRO_H_
