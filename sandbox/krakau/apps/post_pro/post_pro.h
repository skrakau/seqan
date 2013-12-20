#ifndef SANDBOX_KRAKAU_APPS_POST_PRO_POST_PRO_H_
#define SANDBOX_KRAKAU_APPS_POST_PRO_POST_PRO_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fstream>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/bam_io.h>
#include <seqan/score.h> 

//using namespace std;
using namespace seqan;

class Times
{
public:
    double time_all;
    double time_globalAlignment;
    double time_writeBsAlignment;
 
    static Times & instance()
    {
        static Times times;
        return times;
    }

private:
    Times() :
        time_all(0),
        time_globalAlignment(0),
        time_writeBsAlignment(0)
    {}
};


// Realign with 4-letter alphabet and set alignedReadStore entries
template<typename TReadGaps, typename TContigGaps, typename TFragmentStore, typename TId, typename TBsScoreCTLeft, typename TBsScoreCTRight, typename TBsScoreGALeft, typename TBsScoreGARight, typename TOptions>
inline int
reAlign4(TReadGaps &readGaps, 
         TContigGaps &contigGaps, 
         TFragmentStore &store, 
         TId &id, 
         TBsScoreCTLeft &scoringSchemeCTLeft, 
         TBsScoreCTRight &scoringSchemeCTRight, 
         TBsScoreGALeft &scoringSchemeGALeft, 
         TBsScoreGARight &scoringSchemeGARight, 
         TOptions &options)
{ 
    typedef typename Iterator<TReadGaps>::Type      TReadGapsIterator;
    typedef typename Iterator<TContigGaps>::Type    TContigGapsIterator;

    typedef double TValue;
    TValue scoreA;

#ifdef POST_PRO_PROFILE 
    double timeStamp = sysTime();
#endif 
    // Do alignment
    //std::cout << " Do alignment: " << std::endl;
    int band = 6;   // dep. on interOffset?
    if (getMateNo(store, store.alignedReadStore[id].readId) != -1 )     // If paired
    {
        if (getMateNo(store, store.alignedReadStore[id].readId) == 0 && store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
        {
            scoreA = globalAlignment(contigGaps, readGaps, scoringSchemeCTLeft, AlignConfig<true, false, false, true>(), -band, band);
        }    
        else if (getMateNo(store, store.alignedReadStore[id].readId) == 0 && store.alignedReadStore[id].beginPos > store.alignedReadStore[id].endPos)
        {
            scoreA = globalAlignment(contigGaps, readGaps, scoringSchemeGALeft, AlignConfig<true, false, false, true>(), -band, band);
        }  
        else if (getMateNo(store, store.alignedReadStore[id].readId) == 1 && store.alignedReadStore[id].beginPos > store.alignedReadStore[id].endPos)
        {
            scoreA = globalAlignment(contigGaps, readGaps, scoringSchemeCTRight, AlignConfig<true, false, false, true>(), -band, band);
        }     
        else //if (getMateNo(store, store.alignedReadStore[id].readId) == 1 && store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
        {
            scoreA = globalAlignment(contigGaps, readGaps, scoringSchemeGARight, AlignConfig<true, false, false, true>(), -band, band);
        } 
    }
    else
    {
        if (store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
        {
            scoreA = globalAlignment(contigGaps, readGaps, scoringSchemeCTLeft, AlignConfig<true, false, false, true>(), -band, band);
        }
        else
        {
            scoreA = globalAlignment(contigGaps, readGaps, scoringSchemeGALeft, AlignConfig<true, false, false, true>(), -band, band);
        }
    }
#ifdef POST_PRO_PROFILE
    Times::instance().time_globalAlignment += (sysTime() - timeStamp);
#endif
    /*
        std::cout << "align: (reAlign4 1)" << std::endl;
        std::cout << contigGaps << std::endl;
        std::cout << readGaps << std::endl;
    */

    if (store.readNameStore[store.alignedReadStore[id].readId] == "myRead.000000001"
        || store.readNameStore[store.alignedReadStore[id].readId] == "myRead.000424162" 
            || store.readNameStore[store.alignedReadStore[id].readId] == "myRead.000424163" 
            || store.readNameStore[store.alignedReadStore[id].readId] == "myRead.000424164")
    {
        std::cout << "align: (reAlign4 1)" << std::endl;
        std::cout << contigGaps << std::endl;
        std::cout << readGaps << std::endl;
    }
    // Get number of errors (mismatches, indels)
    // Remove gaps at beginning and end of read caused by inaccurate begin and end positions of contig infix 
    TReadGapsIterator itR = begin(readGaps);
    TContigGapsIterator itC = begin(contigGaps);
    while (isGap(itR))
    {
        ++itR;
        ++itC;
    }
    setBeginPosition(readGaps, beginPosition(readGaps));

    if (itC != begin(contigGaps)) // If read has gaps at beginning -> itC != begin(contigGaps) -> set new clipped position; otherwise do nothing 
        setClippedBeginPosition(contigGaps, position(itC));
    TReadGapsIterator itREnd = end(readGaps);
    unsigned countEndGaps = 0;
    while(isGap(--itREnd))    // start behind last position?
        ++countEndGaps;
    
    setEndPosition(readGaps, endPosition(readGaps)); // ? Set the end position of the clipped gapped sequence, given a source position
    if (countEndGaps > 0) setEndPosition(contigGaps, endPosition(contigGaps)-countEndGaps); // Only shorten if there are gaps at end of read, not if gaps at end of contig 


    itREnd = end(readGaps);
    itR = begin(readGaps);      // does this make sense?
    itC = begin(contigGaps);    
    unsigned gaps = 0;
    unsigned mismatches = 0;
    unsigned matches = 0;

   
    if (store.readNameStore[store.alignedReadStore[id].readId] == "myRead.000000001"
        || store.readNameStore[store.alignedReadStore[id].readId] == "myRead.000424162" 
            || store.readNameStore[store.alignedReadStore[id].readId] == "myRead.000424163" 
            || store.readNameStore[store.alignedReadStore[id].readId] == "myRead.000424164")
    {
        std::cout << "align: (reAlign4 2)" << std::endl;
        std::cout << contigGaps << std::endl;
        std::cout << readGaps << std::endl;
    }

    // Set beginPos and endPos new
    if (store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
    {
        if (store.alignedReadStore[id].beginPos < options.intervalOffset) 
            store.alignedReadStore[id].beginPos = 0 + beginPosition(contigGaps);     
        else 
            store.alignedReadStore[id].beginPos = store.alignedReadStore[id].beginPos - options.intervalOffset + beginPosition(contigGaps);     

        store.alignedReadStore[id].endPos = store.alignedReadStore[id].beginPos + endPosition(contigGaps) - beginPosition(contigGaps);  //store.alignedReadStore[id].endPos + options.intervalOffset - countEndGaps;
    }
    else
    {
        if (store.alignedReadStore[id].beginPos < options.intervalOffset) 
            store.alignedReadStore[id].endPos = 0 + beginPosition(contigGaps);
        else
            store.alignedReadStore[id].endPos = store.alignedReadStore[id].endPos - options.intervalOffset + beginPosition(contigGaps);
        store.alignedReadStore[id].beginPos = store.alignedReadStore[id].endPos + endPosition(contigGaps) - beginPosition(contigGaps); //store.alignedReadStore[id].beginPos + options.intervalOffset - countEndGaps;     
    }
    //std::cout<< "Contig positions: " << store.alignedReadStore[id].beginPos << "  end:" << store.alignedReadStore[id].endPos << std::endl; 
    // Count errors
    if (getMateNo(store, store.alignedReadStore[id].readId) != -1 )     // If paired
    {
        if ((getMateNo(store, store.alignedReadStore[id].readId) == 0 && store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos) ||    // Mapped against CT ref
            (getMateNo(store, store.alignedReadStore[id].readId) == 1 && store.alignedReadStore[id].beginPos > store.alignedReadStore[id].endPos) )     
            for (; !atEnd(itR) && !atEnd(itC); ++itR, ++itC)    
            {
                if (isGap(itR) || isGap(itC)) ++gaps;
                else if (*itR == *itC || (*itR == 'T' && *itC == 'C')) ++matches;
                else ++mismatches;
            }
        else    // Mapped against GA ref
            for (; !atEnd(itR) && !atEnd(itC); ++itR, ++itC)    
            {
                if (isGap(itR) || isGap(itC)) ++gaps;
                else if (*itR == *itC || (*itR == 'A' && *itC == 'G') ) ++matches;
                else ++mismatches;
            }
    }
    else
    {
        if (store.alignedReadStore[id].beginPos < store.alignedReadStore[id].endPos)
            for (; !atEnd(itR) && !atEnd(itC); ++itR, ++itC)    
            {
                if (isGap(itR) || isGap(itC)) ++gaps;
                else if (*itR == *itC || (*itR == 'T' && *itC == 'C') ) ++matches;
                else ++mismatches;
            }
        else
            for (; !atEnd(itR) && !atEnd(itC); ++itR, ++itC)    
            {
                if (isGap(itR) || isGap(itC)) ++gaps;
                else if (*itR == *itC || (*itR == 'A' && *itC == 'G')) ++matches;
                else ++mismatches;
            }
    }
    // Rescale score and update store entries
    resize(store.alignQualityStore, length(store.alignedReadStore), Generous());
    store.alignQualityStore[id].errors = gaps + mismatches;
    //store.alignQualityStore[id].score = (scoreA - (double)(mismatches+matches+gaps)*1)*(-10); // scale back to original scaling and convert to phred scale
    store.alignQualityStore[id].score = scoreA*(-10);
    if (store.alignedReadStore[id].readId < 10)
    {
        std::cout << "Score before scaling: " << scoreA << " scale factor: " << -10  
            << "  Score after scaling: " << store.alignQualityStore[id].score
            <<  " #gaps: " << gaps 
            << "  #mismatches: " << mismatches << std::endl;
    }
    store.alignedReadStore[id].gaps = _dataAnchors(readGaps);

    return 0;
}


template <typename TStream, typename TNameStore, typename TNameStoreCache, typename TFragmentStore, typename TContigGaps, typename TId, typename TScore, typename TOptions>
inline int
writeBsAlignment(TStream & stream,
           BamIOContext<TNameStore, TNameStoreCache> const & context,
           TFragmentStore &store, 
           TContigGaps &contigGaps,
           TId &bestId,
           TScore &mapq,
           BamAlignmentRecord &record,
           TOptions &/*options*/)
{
#ifdef POST_PRO_PROFILE 
    double timeStamp = sysTime();
#endif 
    typedef typename TFragmentStore::TContigStore                               TContigStore;
    typedef typename TFragmentStore::TContigPos                                 TContigPos;
    typedef typename Value<TContigStore>::Type                                  TContig;
    typedef typename TContig::TId                                               TContigId;
    typedef typename TContig::TContigSeq                                        TContigSeq;
    //typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> >                    TContigGaps;

    typedef typename TFragmentStore::TAlignedReadStore                          TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                             TAlignedRead;
    typedef typename TAlignedRead::TGapAnchors                                  TReadGapAnchors;
    typedef typename TFragmentStore::TReadSeqStore                              TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq                                   TReadSeq;
    typedef Gaps<TReadSeq, AnchorGaps<TReadGapAnchors> >                        TReadGaps;

    // Create (mate independent) record entries
    record.qName = store.readNameStore[store.alignedReadStore[bestId].readId];
    if (store.alignedReadStore[bestId].beginPos < store.alignedReadStore[bestId].endPos) 
        record.beginPos = store.alignedReadStore[bestId].beginPos;
    else 
    {
        record.beginPos = store.alignedReadStore[bestId].endPos;
        record.flag |= 0x0010;          
    }
    record.mapQ = mapq;

    CharString md;
    String<CigarElement<> > cigar;
    
    TReadSeq readSeq = store.readSeqStore[store.alignedReadStore[bestId].readId]; 

    //TContigPos beginInf;    // TODO unnecessary
    //TContigPos endInf;
    if (store.alignedReadStore[bestId].beginPos < store.alignedReadStore[bestId].endPos)
    {
        //beginInf = store.alignedReadStore[bestId].beginPos;
        //endInf = store.alignedReadStore[bestId].endPos;
        record.seq = readSeq;
    }
    else
    {
        //beginInf = store.alignedReadStore[bestId].endPos;
        //endInf = store.alignedReadStore[bestId].beginPos;
        Dna5String tmpReadSeq = readSeq;
        record.seq = Dna5StringReverseComplement(tmpReadSeq); 
    }
    TReadGaps readGaps(record.seq, store.alignedReadStore[bestId].gaps);

    //TContigGaps contigGaps(infix(store.contigStore[store.alignedReadStore[bestId].contigId].seq, beginInf, endInf), contigGapAnchors);  //

    int printPos = 10;
    if (record.beginPos < printPos)
    {
        std::cout << "align:   errors: " << static_cast<unsigned int>(store.alignQualityStore[bestId].errors) << std::endl;
        //std::cout << "Contig: " <<  contigGaps.data_cutBegin << "  " << contigGaps.data_cutEnd << "  " << contigGaps.data_viewCutBegin << "  " << contigGaps.data_viewCutEnd << std::endl;
        std::cout << "contig: " << contigGaps << std::endl;
        std::cout << "read:   " << readGaps << std::endl;
        std::cout << "quals:  ";
    }

    //std::cout<< "Contig positions: " << store.alignedReadStore[bestId].beginPos << "  end:" << store.alignedReadStore[bestId].endPos << std::endl;
  
    getCigarString(cigar, contigGaps, readGaps);
    getMDString(md, contigGaps, readGaps);
    record.cigar = cigar;

    clear(record.qual);
    resize(record.qual, length(readSeq));
    unsigned avgQual = 0;
    if (store.alignedReadStore[bestId].beginPos < store.alignedReadStore[bestId].endPos)
    {
        unsigned i = 0;
        for (unsigned j = 0; j < length(readSeq); ++j, ++i)
        {
            record.qual[i] = (char)(getQualityValue(readSeq[j]) + 33);
            if (record.beginPos < printPos)
            {
                std::cout << getQualityValue(readSeq[j]) << ",";
            }
            avgQual += getQualityValue(readSeq[j]);
        }
        avgQual = avgQual/(i+1);
    }
    else
    {
        unsigned i = 0;
        for (int j = length(readSeq) -1; j >= 0; --j, ++i)
        {
            record.qual[i] = (char)(getQualityValue(readSeq[j]) + 33);
            if (record.beginPos < printPos)
            {
                std::cout << getQualityValue(readSeq[j]) << ",";
            }
           avgQual += getQualityValue(readSeq[j]);
        }
        avgQual = avgQual/(i+1);
    }
    if (record.beginPos < printPos)
    {
        std::cout << std::endl;
        std::cout << ".............." << '\n';
        std::cout << "mapq: " << static_cast<int>(record.mapQ) << " avgQual: " << avgQual <<  std::endl;
        std::cout << ".............." << std::endl; 
        std::cout << '\n';
    }


    // TODO Use existing function to write record ?
    // Problem with tags...
    int res = 0;
    
#define SEQAN_PUT_TAB                           \
    do {                                        \
        res = streamPut(stream, '\t');      \
        if (res != 0)                       \
            return res;                     \
    } \
    while (false)

    res = streamPut(stream, record.qName);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    res = streamPut(stream, record.flag);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (store.alignedReadStore[bestId].contigId == TAlignedRead::INVALID_ID)
        res = streamPut(stream, '*');
    else
        res = streamPut(stream, nameStore(context)[store.alignedReadStore[bestId].contigId]);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (store.alignedReadStore[bestId].contigId == TAlignedRead::INVALID_ID)
        res = streamPut(stream, '*');
    else
        res = streamPut(stream, record.beginPos + 1);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    res = streamPut(stream, static_cast<int>(record.mapQ));
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (length(record.cigar) == 0u)
    {
        res = streamPut(stream, '*');
        if (res != 0)
            return res;
    }
    else
    {
        for (unsigned i = 0; i < length(record.cigar); ++i)
        {
            res = streamPut(stream, record.cigar[i].count);
            if (res != 0)
                return res;

            res = streamPut(stream, record.cigar[i].operation);
            if (res != 0)
                return res;
        }
    }

    SEQAN_PUT_TAB;

    if (record.rNextId == BamAlignmentRecord::INVALID_REFID)
        res = streamPut(stream, '*');
    else if ((TContigId)store.alignedReadStore[bestId].contigId == (TContigId)record.rNextId)
        res = streamPut(stream, '=');
    else
        res = streamPut(stream, nameStore(context)[record.rNextId]);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (record.pNext == BamAlignmentRecord::INVALID_POS)
        res = streamPut(stream, '0');
    else
        res = streamPut(stream, record.pNext + 1);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (record.tLen == BamAlignmentRecord::INVALID_LEN)
        res = streamPut(stream, '*');
    else 
        res = streamPut(stream, record.tLen);

    if (res != 0)
        return res;

    SEQAN_PUT_TAB;

    if (empty(record.seq))
        res = streamPut(stream, '*');  // Case of empty seq string / "*".
    else
        res = streamPut(stream, record.seq);
    if (res != 0)
        return res;

    SEQAN_PUT_TAB;


    if (empty(record.qual))  // Case of empty quality string / "*".
        res = streamPut(stream, '*');
    else
        res = streamPut(stream, record.qual);
    if (res != 0)
        return res;

/*
    if (length(record.tags) > 0u)
    {
        SEQAN_PUT_TAB;
        CharString buffer;
        assignTagsBamToSam(buffer, record.tags);
        streamPut(stream, buffer);
    }
*/

    streamPut(stream, "\tNM:i:");
    streamPut(stream, static_cast<unsigned int>(store.alignQualityStore[bestId].errors));


    if (!empty(md))
    {
        streamPut(stream, "\tMD:Z:");
        streamPut(stream, md);
    }

    // for testing
    streamPut(stream, "\tDD:i:");
    streamPut(stream, static_cast<unsigned int>(avgQual));


/*
    if (alignedId < length(store.alignedReadTagStore) && !empty(store.alignedReadTagStore[alignedId]))
    {
        _streamPut(target, '\t');
        _streamWrite(target, store.alignedReadTagStore[alignedId]);
    }
*/ // TODO what about other tags? read into alignedReadTagStore ?? alignment independent information?

#ifdef POST_PRO_PROFILE
    Times::instance().time_writeBsAlignment += (sysTime() - timeStamp);
#endif

    return streamPut(stream, '\n');

#undef SEQAN_PUT_TAB
}

// TODO Does it make sense to count errors? (with score it would have to be exactly the same alignment then) threshold?

template<typename TFragmentStore, typename TId>
inline unsigned
countEqualHits(TFragmentStore &store, TId &id)
{
    // Count number of hits which have the same number of errors as second best hit 
    unsigned count = 0;
    unsigned errors = store.alignQualityStore[id].errors;
    for (unsigned i = 0; i < length(store.alignQualityStore); ++i)
            if (store.alignQualityStore[i].errors == errors)
                ++count;
    return count;
}

// Compute pseud worst score for read
// Assuming errors are at best quality positions
// No gaps assumed
// Take max. allowed errors + 1
template<typename TScore, typename TReadSeq, typename TOptions>
inline void 
computePseudoWorstScore(TScore &score, TReadSeq &readSeq, TOptions &options)
{
    // Get best qualitites
    
    String<int> bestQuals;
    String<int> lowQuals;
    unsigned max3Errors = floor(options.max3Error * length(readSeq)/ 100.0);
    resize(bestQuals, max3Errors + 1, 0);
    for (unsigned i = 0; i < length(readSeq); ++i)
    {
        // Quite inefficient, but in practice there will be only around 2 errors
        // So sorting will be fast
        if (getQualityValue(readSeq[i]) > back(bestQuals))
        {
            int j = max3Errors - 1;
            for(; j >= 0; --j)
            {
                if (getQualityValue(readSeq[i]) < bestQuals[j]) continue;
            }
            insert(bestQuals, j+1, getQualityValue(readSeq[i]));
            if (back(bestQuals) > 0 ) appendValue(lowQuals, back(bestQuals));             // move() ?
            eraseBack(bestQuals);
         }
        else
            appendValue(lowQuals, getQualityValue(readSeq[i]));
    }
    // Compute pseudo worst score
    score = 0;
    for (unsigned i = 0; i < length(bestQuals); ++i)
    {
        //long double e = pow(10, -(long double)bestQuals[i]/10.0);
        //score += std::log10(10.0*(e/3.0));  // Assume mismatch
        score += -1.0; //-4.0;  // Assume gap opening (worst score at the moment) TODO
        //std::cout << " mismatch: " << std::log10(10.0*(e/3.0)) << std::endl;
    }
    for (unsigned i = 0; i < length(lowQuals); ++i)
    {
            long double e = pow(10, -(long double)lowQuals[i]/10.0);
        if (e == 1) e = 0.9;    // e.g. in case of Ns: take as low probability into account
        score += std::log10(options.scoreMatch*(1.0 - e) + (options.scoreMismatch)*e);    // Assume match
    }
    SEQAN_ASSERT_EQ(length(readSeq), length(bestQuals)+length(lowQuals));
    //score -= length(readSeq);
    score *= -10;
}
// Use avg. quality, faster to compute and shouldn't make such a big difference, if we just use -1 for mismatches and quality only for matches
template<typename TScore, typename TReadSeq, typename TOptions>
inline void 
computePseudoWorstScore2(TScore &score, TReadSeq &readSeq, TOptions &options)
{
    double avgQual = 0;
    for (unsigned i = 0; i < length(readSeq); ++i)
    {
        avgQual += getQualityValue(readSeq[i]);
    }
    avgQual = avgQual/(double)length(readSeq);
    // Compute pseudo worst score
    unsigned max3Errors = floor(options.max3Error * length(readSeq)/ 100.0);
    unsigned pseudoErrors = max3Errors + 1;
    score = 0;
    score += pseudoErrors * (-1.0);   // Assumed mismatches
    double e = pow(10, -(long double)avgQual/10.0);
    score += (length(readSeq) - pseudoErrors) * std::log10(options.scoreMatch*(1.0 - e) + (options.scoreMismatch)*e);   // Assumed matches
    score -= length(readSeq);
    //score *= -10;
} 


struct VerifiedRead
{
    typedef FragmentStore<MyFragmentStoreConfig>::TMappingQuality          TScore;
    typedef FragmentStore<MyFragmentStoreConfig>::TAlignedReadStore        TAlignedReadStore;
    typedef Value<TAlignedReadStore>::Type                                 TAlignedReadStoreElement;
    typedef TAlignedReadStoreElement::TId                                  TId;

    TId alignedReadId;
    BamAlignmentRecord record;
    TScore mapq;

	VerifiedRead():
		alignedReadId(0),
		mapq(0) {}

	VerifiedRead(TId _alignedReadId, BamAlignmentRecord _record, TScore _mapq):
		alignedReadId(_alignedReadId),
		record(_record),
		mapq(_mapq) {}

};

template<typename TFragmentStore, typename TOptions>
inline VerifiedRead
verifyRead(TFragmentStore &store, TOptions &options)
{
    typedef typename TFragmentStore::TAlignQualityStore     TAlignQualityStore;
    typedef typename TFragmentStore::TMatePairStore         TMatePairStore;
    typedef typename Value<TMatePairStore>::Type            TMatePairStoreElement;
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedReadStoreElement;
    typedef typename TFragmentStore::TReadSeq               TReadSeq;
    typedef typename TFragmentStore::TMappingQuality        TScore;

    // Find best and second best hit (match mate hit) -> lowest score
    TScore bestScore = store.alignQualityStore[0].score;
    unsigned bestId = 0;
    TScore secBestScore = store.alignQualityStore[0].score;
    unsigned secBestId = 0;
    for (unsigned i = 1; i < length(store.alignedReadStore); ++i)      
    {
        if (bestScore > store.alignQualityStore[i].score)
        {
            secBestScore = bestScore;
            secBestId = bestId;
            bestScore = store.alignQualityStore[i].score;
            bestId = i;
        }
        else if (secBestScore > store.alignQualityStore[i].score)
        {
            secBestScore = store.alignQualityStore[i].score;
            secBestId = i;
        }
    }

    unsigned countHits = length(store.alignedReadStore);
    TScore mapq;
    computeMapq(mapq, bestId, countHits, secBestId, store, options);

    BamAlignmentRecord record;
    record.flag = 0;
    record.rNextId = BamAlignmentRecord::INVALID_REFID;
    record.pNext = BamAlignmentRecord::INVALID_POS;
    record.tLen = 0;

    VerifiedRead verifiedRead(bestId, record, mapq);
    return verifiedRead;
}

template<typename TFragmentStore, typename TId>
inline unsigned
countEqualMateHits(TFragmentStore &store, TId &mateId)
{
    // Count number of hits of one mate which have the same number of errors 
    unsigned count = 0;
    unsigned errors = store.alignQualityStore[mateId].errors;
    if (getMateNo(store, store.alignedReadStore[mateId].readId) == 0)
    {
        for (unsigned i = 0; i < length(store.alignQualityStore); ++i)
            if (getMateNo(store, store.alignedReadStore[i].readId) == 0 && store.alignQualityStore[i].errors == errors)
                ++count;
    }
    else if (getMateNo(store, store.alignedReadStore[mateId].readId) == 1)
    {
        for (unsigned i = 0; i < length(store.alignQualityStore); ++i)
            if (getMateNo(store, store.alignedReadStore[i].readId) == 1 && store.alignQualityStore[i].errors == errors)
                ++count;
    }
    return count;
}

template<typename TFragmentStore, typename TId>
inline unsigned
countEqualMateHits(TFragmentStore &store, TId &mateIdL, TId &mateIdR)
{
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedReadStoreElement;
    // Count number of match mate hits which have the same number of errors
    unsigned count = 0;
   
    unsigned errors = store.alignQualityStore[mateIdL].errors + store.alignQualityStore[mateIdR].errors;

    for (unsigned i = 0; i < length(store.alignQualityStore); ++i)
        if (getMateNo(store, store.alignedReadStore[i].readId) == 0 && store.alignedReadStore[i].pairMatchId != TAlignedReadStoreElement::INVALID_ID) // Only if left mate and has match mate
            if (store.alignQualityStore[i].errors + store.alignQualityStore[store.alignedReadStore[i].pairMatchId].errors == errors)
                ++count;

    return count;
}


template<typename TScore, typename TId, typename TFragmentStore, typename TOptions>
inline void
computeMapq(TScore &mapq, TId &bestId, unsigned countHits, TId &secBestId, TFragmentStore &store, TOptions &options)
{
    typedef typename TFragmentStore::TReadSeq               TReadSeq;

    unsigned countBest;
    unsigned countSecBest;
    // If only one pair hit: use single hit counts
    // If multiple pair hits: use pair hit counts 
    if (!empty(store.matePairStore))    
    {
        countBest = countEqualMateHits(store, bestId);
        countSecBest = countEqualMateHits(store, secBestId);
    }
    else
    {
        countBest = countEqualHits(store, bestId);
        countSecBest = countEqualHits(store, secBestId);
    }
    unsigned printPos = 10;
    unsigned pos = (store.alignedReadStore[bestId].beginPos < store.alignedReadStore[bestId].endPos)? store.alignedReadStore[bestId].beginPos:store.alignedReadStore[bestId].endPos;

    if (countHits == 1)
    {
        if ( pos< printPos)
        {
            std::cout << " Only one hit" << std::endl;
        }
        TScore worstScore;
        TReadSeq readSeq = store.readSeqStore[store.alignedReadStore[bestId].readId]; 
        computePseudoWorstScore2(worstScore, readSeq, options);
        mapq = worstScore - store.alignQualityStore[bestId].score;
        //std::cout << "**************computeMapq: worstScore :" << worstScore << " score: " << store.alignQualityStore[bestId].score << std::endl;
    }
    else if (countBest != 1) // Same edit distance
    {
        if (pos < printPos)
        {
            std::cout << " Non-unique" << std::endl;
        }
        mapq = 0;
    }
    else
    {
        mapq = store.alignQualityStore[secBestId].score  - store.alignQualityStore[bestId].score - 10*std::log10(countSecBest);
        if (pos < printPos)
        {
            std::cout << " Unique: " << std::endl;
            std::cout << "***: sec best score :" << store.alignQualityStore[secBestId].score  << " best score: " << store.alignQualityStore[bestId].score << "  count2: " << countSecBest << std::endl;
        }
        //  can be negative, since it's possible that there exist single hits, which are better than left or right mate match hit
    }

    // Q: look at mate pairs as one whole?
    // secBest -> secBest in secBestMate Hit ?
    // No, could be that one read is mapped correct, and other one not
}

//template<>
struct VerifiedMates
{
    typedef FragmentStore<MyFragmentStoreConfig>::TMappingQuality          TScore;
    typedef FragmentStore<MyFragmentStoreConfig>::TAlignedReadStore        TAlignedReadStore;
    typedef Value<TAlignedReadStore>::Type                                 TAlignedReadStoreElement;
    typedef TAlignedReadStoreElement::TId                                  TId;

    TId alignedReadIdL;
    TId alignedReadIdR;
    BamAlignmentRecord recordL;
    BamAlignmentRecord recordR;
    TScore mapqL;
    TScore mapqR;

	VerifiedMates():
		alignedReadIdL(0),
		alignedReadIdR(0),
		mapqL(0),
		mapqR(0) {}

	VerifiedMates(TId _alignedReadIdL, 
                  TId _alignedReadIdR, 
	              BamAlignmentRecord _recordL,
	              BamAlignmentRecord _recordR,
	              TScore _mapqL,
	              TScore _mapqR):
		alignedReadIdL(_alignedReadIdL),
		alignedReadIdR(_alignedReadIdR),
		recordL(_recordL),
		recordR(_recordR),
		mapqL(_mapqL), 
		mapqR(_mapqR) {}
};


template<typename TFragmentStore, typename TOptions>
inline VerifiedMates
verifyMates(TFragmentStore &store, TOptions &options)
{
    //std::cout << " VerifyMates...." << std::endl;
    typedef typename TFragmentStore::TAlignQualityStore     TAlignQualityStore;
    typedef typename TFragmentStore::TMatePairStore         TMatePairStore;
    typedef typename Value<TMatePairStore>::Type            TMatePairStoreElement;
    typedef typename TFragmentStore::TAlignedReadStore      TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type         TAlignedReadStoreElement;
    typedef typename TAlignedReadStoreElement::TId          TId;
    typedef typename TFragmentStore::TContigPos             TContigPos;
    typedef typename TFragmentStore::TMappingQuality        TScore;
    typedef typename TFragmentStore::TReadSeqStore          TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type             TReadSeq;
    
    // Find best and second best hit (match mate hit) -> lowest score
    TScore bestScoreL = store.alignQualityStore[0].score;
    TId bestIdL = 0;
    TScore secBestScoreL = store.alignQualityStore[0].score;
    TId secBestIdL = 0;
    TScore bestScoreR = store.alignQualityStore[0].score;
    TId bestIdR = 0;
    TScore secBestScoreR = store.alignQualityStore[0].score;
    TId secBestIdR = 0;

    unsigned countPairHits = 0;
    unsigned countHitsL = 0;
    unsigned countHitsR = 0;
    TScore bestMateScore = 1000000;   // Set to maximum
    TId bestMateId = 0;
    TScore secBestMateScore = 1000000;
    TId secBestMateId = 0;
    if (store.matePairStore[0].readId[0] != TMatePairStoreElement::INVALID_ID && store.matePairStore[0].readId[1] != TMatePairStoreElement::INVALID_ID) // If read is paired
    {
        for (TId i = 0; i < length(store.alignedReadStore); ++i)    
        {
            // Get single left and right scores
            if (getMateNo(store, store.alignedReadStore[i].readId) == 0 )
            {
                ++countHitsL;
                if (bestScoreL >= store.alignQualityStore[i].score)
                {
                    secBestScoreL = bestScoreL;
                    secBestIdL = bestIdL;
                    bestScoreL = store.alignQualityStore[i].score;
                    bestIdL = i;
                }
                else if (secBestScoreL >= store.alignQualityStore[i].score)
                {
                    secBestScoreL = store.alignQualityStore[i].score;
                    secBestIdL = i;
                }
            }
            else if (getMateNo(store, store.alignedReadStore[i].readId) == 1 )
            {
                ++countHitsR;
                if (bestScoreR >= store.alignQualityStore[i].score)
                {
                    secBestScoreR = bestScoreR;
                    secBestIdR = bestIdR;
                    bestScoreR = store.alignQualityStore[i].score;
                    bestIdR = i;
                }
                else if (secBestScoreR >= store.alignQualityStore[i].score)
                {
                    secBestScoreR = store.alignQualityStore[i].score;
                    secBestIdR = i;
                }
            }
            // Get mate match scores
            //if (store.alignedReadStore[i].pairMatchId == TAlignedReadStoreElement::INVALID_ID)
                //std::cout << "      i: " << i << "   pairMatchId INVALID  " <<  getMateNo(store, store.alignedReadStore[i].readId) << std::endl;

            if (getMateNo(store, store.alignedReadStore[i].readId) == 0 && store.alignedReadStore[i].pairMatchId != TAlignedReadStoreElement::INVALID_ID ) // Only if left mate and has match mate
            {
                ++countPairHits;

                double currMateScore = store.alignQualityStore[i].score + store.alignQualityStore[store.alignedReadStore[i].pairMatchId].score;
                //std::cout << " i: " << i << " currMatescore: " << currMateScore << "  bestMateScore " << bestMateScore << std::endl;
                if (currMateScore <= bestMateScore)
                {
                    secBestMateScore = bestMateScore;
                    secBestMateId = bestMateId;             // Id of left aligned read
                    bestMateScore = currMateScore;
                    bestMateId = i;
                    //std::cout << " i: " << i << std::endl;
                }
                else if (currMateScore <= secBestMateScore)
                {
                    secBestMateScore = currMateScore;
                    secBestMateId = i;
                }
            }
        }
    }
    //std::cout << "      countPairHits: " << countPairHits << "  countHitsL  " << countHitsL << "  countHitsR  " << countHitsR << std::endl;

    // Compute mapping qualities dependent on case
    unsigned count1 = 0;
    if (countPairHits >= 1) 
        count1 = countEqualMateHits(store, bestMateId, store.alignedReadStore[bestMateId].pairMatchId);
    //std::cout << "countPairHits: " << countPairHits << "  count1: " << count1 << std::endl;
    TScore mapqL;
    TScore mapqR;
    BamAlignmentRecord recordL;
    BamAlignmentRecord recordR;
    if (countPairHits == 1)  // Case1: Only one pair hit (single hits can be nonunique)
    {
        TId bestMateIdR = store.alignedReadStore[bestMateId].pairMatchId;
        computeMapq(mapqL, bestMateId, countHitsL, secBestIdL, store, options);
        computeMapq(mapqR, bestMateIdR, countHitsR, secBestIdR, store, options);

        if (store.readNameStore[store.alignedReadStore[bestMateId].readId] == "myRead.000000001" || 
                store.readNameStore[store.alignedReadStore[bestMateId].readId] == "myRead.000424162" || 
                store.readNameStore[store.alignedReadStore[bestMateId].readId] == "myRead.000424163" || 
                store.readNameStore[store.alignedReadStore[bestMateId].readId] == "myRead.000424164")
        {
            std::cout << "One hit: bestMateId: " << bestMateId << "  bestMateIdR: " << bestMateIdR << std::endl;
            std::cout << "mapqL: " << mapqL << "  mapqR: " << mapqR << std::endl;
        }

        mapqL = mapqL + mapqR;
        mapqR = mapqL;

        // Set record entries with matepair information
        recordL.flag = 0;
        recordR.flag = 0;
        recordL.flag |= 0x001;
        recordR.flag |= 0x001;
        recordL.flag |= 0x002;
        recordR.flag |= 0x002;
        recordL.flag |= 0x0040;  // Is left read
        recordR.flag |= 0x0080;  // Is right read
        recordL.rNextId = store.alignedReadStore[bestMateIdR].contigId;
        recordR.rNextId = store.alignedReadStore[bestMateId].contigId;

        if (store.alignedReadStore[bestMateId].beginPos < store.alignedReadStore[bestMateId].endPos) 
        {
            recordR.pNext = store.alignedReadStore[bestMateId].beginPos;
        }
        else 
        {
            recordR.flag |= 0x0020;
            recordR.pNext = store.alignedReadStore[bestMateId].endPos;
        }
        if (store.alignedReadStore[bestMateIdR].beginPos < store.alignedReadStore[bestMateIdR].endPos) 
        {
            recordL.pNext = store.alignedReadStore[bestMateIdR].beginPos;
        }
        else 
        {
            recordL.flag |= 0x0020; // Mate is reverse complemented
            recordL.pNext = store.alignedReadStore[bestMateIdR].endPos;
        }
        TContigPos beginPair = store.alignedReadStore[bestMateId].beginPos;
        TContigPos endPair = store.alignedReadStore[bestMateIdR].beginPos;
        if (beginPair < endPair)
        {
            recordL.tLen = endPair - beginPair; 
            recordR.tLen = -(endPair - beginPair);  
        }
        else
        {
            recordL.tLen = -(beginPair - endPair); 
            recordR.tLen = beginPair - endPair;        
        }

        //std::cout << "One hit:  mapqL: " << mapqL << "  errors: " << static_cast<unsigned>(store.alignQualityStore[bestMateId].errors) << "  score: " << store.alignQualityStore[bestMateId].score << std::endl;
        //std::cout << "  right:  mapqR: " << mapqR << "  errors: " << static_cast<unsigned>(store.alignQualityStore[bestMateIdR].errors) << "  score: " << store.alignQualityStore[bestMateIdR].score << std::endl;

        VerifiedMates verifiedMates(bestMateId, bestMateIdR, recordL, recordR, mapqL, mapqR);
        return verifiedMates;

    }
    else if (countPairHits > 1)   // Case 2: Multiple mate matches (no summing up of mapqs)
    {
        if (true) // this is realy the best...?
        {
            // Compare best mate match score against second best mate match (not single) 
            // -> Ensure that mate match get higher mapping quality than it would get in not paired mode
            // -> Increases prob. that score of each single mate is higher that that of second best mate match
            // (is mate match is classified as best, doesn't mean there are no single reads which are better than single mates)
            TId bestMateIdR = store.alignedReadStore[bestMateId].pairMatchId;
            TId secBestMateIdR = store.alignedReadStore[secBestMateId].pairMatchId;

            computeMapq(mapqL, bestMateId, countHitsL, secBestMateId, store, options);
            computeMapq(mapqR, bestMateIdR, countHitsR, secBestMateIdR, store, options);

            if (store.readNameStore[store.alignedReadStore[bestMateId].readId] == "myRead.000000001" || 
                store.readNameStore[store.alignedReadStore[bestMateId].readId] == "myRead.000424162" || 
                store.readNameStore[store.alignedReadStore[bestMateId].readId] == "myRead.000424163" || 
                store.readNameStore[store.alignedReadStore[bestMateId].readId] == "myRead.000424164")
            {
                std::cout << "Multi hits: bestMateId: " << bestMateId << "  bestMateIdR: " << bestMateIdR << std::endl;
                std::cout << "mapqL: " << mapqL << "  mapqR: " << mapqR << std::endl;
                std::cout << "Mutiple hits:  mapqL: " << mapqL << "  errors: " << static_cast<unsigned>(store.alignQualityStore[bestMateId].errors) << "  score: " << store.alignQualityStore[bestMateId].score << std::endl;
                std::cout << "  right:  mapqR: " << mapqR << "  errors: " << static_cast<unsigned>(store.alignQualityStore[bestMateIdR].errors) << "  score: " << store.alignQualityStore[bestMateIdR].score << std::endl;
            }

            // If pair hit is unique
            if (count1 == 1) 
            {
                mapqL = mapqL + mapqR;
                mapqR = mapqL;
            }
        

            // Set record entries with matepair information
            recordL.flag = 0;
            recordR.flag = 0;
            recordL.flag |= 0x001;
            recordR.flag |= 0x001;
            recordL.flag |= 0x002;
            recordR.flag |= 0x002;
            recordL.flag |= 0x0040;  // Is left read
            recordR.flag |= 0x0080;  // Is right read 
            recordL.rNextId = store.alignedReadStore[bestMateIdR].contigId;
            recordR.rNextId = store.alignedReadStore[bestMateId].contigId;

            if (store.alignedReadStore[bestMateId].beginPos < store.alignedReadStore[bestMateId].endPos) 
            {
                recordR.pNext = store.alignedReadStore[bestMateId].beginPos;
            }
            else 
            {
                recordR.flag |= 0x0020;
                recordR.pNext = store.alignedReadStore[bestMateId].endPos;
            }
            if (store.alignedReadStore[bestMateIdR].beginPos < store.alignedReadStore[bestMateIdR].endPos) 
            {
                recordL.pNext = store.alignedReadStore[bestMateIdR].beginPos;
            }
            else 
            {
                recordL.flag |= 0x0020; // Mate is reverse complemented
                recordL.pNext = store.alignedReadStore[bestMateIdR].endPos;
            }
            TContigPos beginPair = store.alignedReadStore[bestMateId].beginPos;
            TContigPos endPair = store.alignedReadStore[bestMateIdR].beginPos;
            if (beginPair < endPair)
            {
                recordL.tLen = endPair - beginPair; 
                recordR.tLen = -(endPair - beginPair);  
            }
            else
            {
                recordL.tLen = -(beginPair - endPair); 
                recordR.tLen = beginPair - endPair;        
            }

            VerifiedMates verifiedMates(bestMateId, bestMateIdR, recordL, recordR, mapqL, mapqR);
            return verifiedMates;

        }
        else    // There is a single hit which is better than the mate match
        {

        }
    }
    else
    {
        SEQAN_ASSERT_EQ(true, false);
        std::cout << " ....else" << std::endl;
        VerifiedMates verifiedMates;
        return verifiedMates;
    }
    /*else if (options.outputSingleMates && countPairHits == 0)   // Case 3: Only single matches // input of single reads by mapper possible? 
    {
        if (countHitsL > 0)
        {
            // Compute mapqs
            unsigned count1L = countEqualHits(store, bestIdL);
            unsigned count2L = countEqualHits(store, secBestIdL);
            if (countHitsL == 1)
                mapqL = 255 - bestScoreL; //-10*log10(1-pow(10, -bestScoreL/10.0)); 
            else if (count1L != 1)
                mapqL = 0;
            else
                mapqL = secBestScoreL - bestScoreL - 10*std::log10(count2L); 

            // Set record entries with matepair information
            recordL.flag = 0;
            recordL.flag |= 0x008;   // Mate is unmapped (or not paired with this one)        
            recordL.flag |= 0x0040;  // Is left read 
            recordL.rNextId = '*';
            recordL.pNext = '*';
            recordL.tLen = 0; 

        }
        if (countHitsR > 0)
        {
            // Compute mapqs
            unsigned count1R = countEqualHits(store, bestIdR);
            unsigned count2R = countEqualHits(store, secBestIdR);
            if (countHitsR == 1)
                mapqR = 255 - bestScoreR; //-10*log(1-pow(10, -bestScoreR/10.0)); 
            else if (count1R != 1)
                mapqR = 0;
            else
                mapqR = secBestScoreR - bestScoreR - 10*std::log10(count2R);

            // Set record entries with matepair information
            recordR.flag = 0;
            recordR.flag |= 0x008;   // Mate is unmapped (or not paired with this one)
            recordR.flag |= 0x0080;  // Is right read 
            recordL.rNextId = '*';
            recordL.pNext = '*';
            recordR.tLen = 0;  

        }
    }*/
}


template <typename TAlignedRead, typename TMInfo, typename TFragStore, typename TOptions>
inline int 
_compareAlignedReadAndMateInfo2(TAlignedRead const &a, TMInfo const &b, TFragStore const &fragStore, TOptions &options)
{
    if (a.contigId < b.contigId) return -1;
    if (a.contigId > b.contigId) return 1;

    typename TFragStore::TContigPos posA = _min(a.beginPos, a.endPos);
    if (posA < b.beginPos - options.intervalOffset)     // TODO more?
    {
        //std::cout << " posA " << posA << " b.beginPos: " << b.beginPos << std::endl;
        return -1;
    }
    if (posA > b.beginPos + options.intervalOffset) return 1;
    
    bool reversedA = (a.beginPos > a.endPos);
    if (!reversedA && b.reversed) 
    {
        //std::cout << "Reversed..." << std::endl;
        return -1;
    }
    if (reversedA && !b.reversed) return 1;

    typedef typename TFragStore::TMatePairStore     TMatePairStore;
    typedef typename Value<TMatePairStore>::Type    TMatePair;
    typename TMatePair::TId matePairIdA = TMatePair::INVALID_ID;

    if (a.readId < length(fragStore.readStore))
        matePairIdA = fragStore.readStore[a.readId].matePairId;
        
    if (matePairIdA < b.matePairId) 
    {
        //std::cout << "MAtePairId..." << std::endl;
        return -1;
    }
    if (matePairIdA > b.matePairId) return 1;
    return 0;
}

// Use pairMatchIds to store alignId of other match mate, which can still be used to access this entry in alignedReadStore (if order not touched since building!) 
template<typename TFragmentStore, typename TMatchMateInfos, typename TOptions>
inline void 
_generatePseudoPairMatchIds (
    TFragmentStore &store,
    TMatchMateInfos & matchMateInfos,
    TOptions &options)
{
    typedef typename TFragmentStore::TAlignedReadStore              TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                 TAlignedRead;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TIter;    
    typedef typename Iterator<TMatchMateInfos, Standard>::Type      TMIter;    
            
    TIter it = begin(store.alignedReadStore, Standard());
    TIter itEnd = end(store.alignedReadStore, Standard());
    TMIter mit = begin(matchMateInfos, Standard());
    TMIter mitEnd = end(matchMateInfos, Standard());

    if (it == itEnd || mit == mitEnd) return;

    // sort the aligned read store by: begin position, contig name
    //std::sort(it,  itEnd,  AlignedMateLess_<TFragmentStore>(fragStore)); // Should be sorted by position for one read anyway
    std::sort(mit, mitEnd, MatchMateInfoLess_());

    while (true)
    {
        // skip already aligned reads
        while (it->pairMatchId != TAlignedRead::INVALID_ID)
            if (++it == itEnd) return;

        int cmp = _compareAlignedReadAndMateInfo2(*it, *mit, store, options);
        if (cmp == 0)   // both are equal -> link them
        {
            (*it).pairMatchId = (*mit).pairMatchId;
            store.alignedReadStore[(*it).pairMatchId].pairMatchId = (*it).id; // Only possible if order of alignedReadStore not touched, alignId == position
        }
        if (cmp >= 0)   // MateInfo is less or equal
        {
            if (++mit == mitEnd) return;
        }
        if (cmp <= 0)   // AlignedRead is less or equal
        {
            if (++it == itEnd) return;
        }
    }
} 


template<typename TOptions, typename TModel>
inline bool
postProcessMain(TOptions &options, TModel const &)
{
    typedef Align<Dna5String,ArrayGaps> TAlign; 
    typedef FragmentStore<MyFragmentStoreConfig> TFragmentStore;

    typedef typename TFragmentStore::TReadNameStore            TReadNameStore;
    typedef NameStoreCache<TReadNameStore, CharString>          TReadNameStoreCache;

    typedef typename TFragmentStore::TContigNameStore                              TContigNameStore;
    typedef NameStoreCache<TContigNameStore, CharString>                            TContigNameStoreCache;
    typedef typename TFragmentStore::TContigStore                                  TContigStore;
    typedef typename Value<TContigStore>::Type                                      TContig;
    typedef typename TContig::TContigSeq                                            TContigSeq;
    typedef typename TFragmentStore::TReadSeq                                      TReadSeq;
    typedef typename TFragmentStore::TAlignedReadStore                             TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                                 TAlignedRead;
 
    typedef BamIOContext<TContigNameStore, TContigNameStoreCache>                   TBamIOContext;
    typedef typename TAlignedRead::TGapAnchors                                      TReadGapAnchors;
    typedef Gaps<TReadSeq, AnchorGaps<TReadGapAnchors> >                            TReadGaps;
    typedef typename TContig::TGapAnchors                                           TContigGapAnchors;
    typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> >                        TContigGaps;
    typedef typename TFragmentStore::TContigPos                                    TContigPos;
    
    typedef typename TAlignedRead::TId     TId;
    
    typedef StringSet<TContigGapAnchors> TSetContigGapAnchors;  // TODO unneccessary ?
    typedef StringSet<TContigGaps> TSetContigGaps;

    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type          TMatePairStoreElement;
    typedef MatchMateInfo_<TContigPos, TId>                                         TMatchMateInfo;
    typedef String<TMatchMateInfo>                                                  TMatchMateInfos;

    // Initialize aligment scores
    typedef double  TValue;
    typedef Score<TValue, BsTagList<BsCaseCT, TModel, Left> >           TBsScoreCTLeft;
    typedef Score<TValue, BsTagList<BsCaseCT, TModel, Right> >          TBsScoreCTRight;
    typedef Score<TValue, BsTagList<BsCaseGA, TModel, Left> >           TBsScoreGALeft;
    typedef Score<TValue, BsTagList<BsCaseGA, TModel, Right> >          TBsScoreGARight;

    BsSubstitutionMatrix<TValue, BsCaseCT, BsSimple> bsSubstitutionMatrixCT(options.globalMethRate, options.bsConversionRate);
    BsSubstitutionMatrix<TValue, BsCaseGA, BsSimple> bsSubstitutionMatrixGA(options.globalMethRate, options.bsConversionRate);
    TValue const * seqErrorFreqs; 
    TValue const * insErrorFreqs;
    TValue const * delErrorFreqs;
    if (options.nonSimpleSubstErrors) seqErrorFreqs = SeqErrorFreqs<TValue, BsNonSimple>::getData();
    else seqErrorFreqs = SeqErrorFreqs<TValue, BsSimple>::getData();
    if (options.nonSimpleInsErrors) insErrorFreqs = InsErrorFreqs<TValue, BsNonSimple>::getData();
    else insErrorFreqs = InsErrorFreqs<TValue, BsSimple>::getData();
    if (options.nonSimpleDelErrors)
    {
        delErrorFreqs = DelErrorFreqs<TValue, BsNonSimple>::getData();
        options.scalingFactorDelErrors = 3.5;       // Compute automatically...
    }
    else delErrorFreqs = DelErrorFreqs<TValue, BsSimple>::getData();

    TBsScoreCTLeft  scoringSchemeCTLeft(options, bsSubstitutionMatrixCT, seqErrorFreqs, insErrorFreqs, delErrorFreqs);
    TBsScoreCTRight scoringSchemeCTRight(options, bsSubstitutionMatrixCT, seqErrorFreqs, insErrorFreqs, delErrorFreqs);
    TBsScoreGALeft  scoringSchemeGALeft(options, bsSubstitutionMatrixGA, seqErrorFreqs, insErrorFreqs, delErrorFreqs);
    TBsScoreGARight scoringSchemeGARight(options, bsSubstitutionMatrixGA, seqErrorFreqs, insErrorFreqs, delErrorFreqs);
    std::cout << "Computation of alignment scores finished." << std::endl;
        
    // Load original reads into fragmentStore
    TFragmentStore store;
    if (options.readFileName2 == "")
        loadReadsCroppedId(store, options.readFileName);
    else
        loadReadsCroppedId(store, options.readFileName, options.readFileName2);

    TReadNameStoreCache readNameCache(store.readNameStore);
    refresh(store.readNameStoreCache);

    // Load reference into fragmentStore
    loadContigs(store, options.refFileName);
    refresh(store.contigNameStoreCache);

    // Parse SAM file ....
    std::fstream inStream(toCString(options.samFileName), std::ios::binary | std::ios::in);
    if (!inStream.good())
    {
        std::cerr << "ERROR: Could not open " << toCString(options.samFileName) << " for reading.\n";
        return 1;
    }
    RecordReader<std::fstream, seqan::SinglePass<> > reader(inStream);
    std::fstream outStream(toCString(options.outputFileName), std::ios::binary | std::ios::out);
    if (!outStream.good())
    {
        std::cerr << "ERROR: Could not open " << toCString(options.outputFileName) << " for writing.\n";
        return 1;
    }

    BamHeader header;
    BamAlignmentRecord record;
    TBamIOContext bamIOContext(store.contigNameStore, store.contigNameStoreCache);
    // Read header
    if (readRecord(header, bamIOContext, reader, Sam()) != 0)  // ATTENTIONE: adds contigName to contigNameStore if not existing already...! 
    {
        std::cerr << "ERROR: Could not read SAM header record!\n";
        return 1;
    }  
    // Write out header again. Maybe add information about BS mapping ?
    if (write2(outStream, header, bamIOContext, Sam()) != 0)
    {
        std::cerr << "ERROR: Could not write header to SAM file "  << "\n";
        return 1;
    }

    // Parse SAM file, verify and write to output
    // data structure to temporarily store information about match mates
    TMatchMateInfos matchMateInfos;

    CharString currReadName;
    TId readId;
    //TSetContigGapAnchors setContigGapAnchors;   // TODO only store anchorGaps and assign contig seq later for output
    TSetContigGaps setContigGaps;

    int helper = 0;
    while (!atEnd(reader))
    {
        ++helper;
        if (readRecord(record, bamIOContext, reader, Sam()) != 0)  // // ATTENTIONE: adds contigName to contigNameStore if not existing already...!
        {
            std::cerr << "ERROR: Could not read sam record!\n";
            return 1;
        }       
        if (record.qName != currReadName)   // First entry for read: Verify previous read and read in new read
        {
            //std::cout << "new read name: " << length(store.alignedReadStore) << std::endl;
            if (!empty(store.alignedReadStore))   // After all entries for one read are loaded
            {
                //std::cout << "Length alignedReadStore: " << length(store.alignedReadStore) << std::endl;
                if (!empty(store.matePairStore))
                {
                    // set the match mate IDs using the information stored in matchMateInfos
                    _generatePseudoPairMatchIds(store, matchMateInfos, options);
                    
                    VerifiedMates verifiedMates = verifyMates(store, options);  // For mates: Find best alignments and compute mapq, verify
                    //std::cout << "alignedReadIdL: " << verifiedMates.alignedReadIdL << "  alignedReadIdr: " << verifiedMates.alignedReadIdR << "  length(alignedReadStore) " <<  length(store.alignedReadStore) << std::endl;
                    //std::cout << "Contig gaps: "  << std::endl;
                    //std::cout << setContigGaps[verifiedMates.alignedReadIdL] << std::endl;
                    unsigned max4Errors = floor(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedMates.alignedReadIdL].readId])/100.0);
                    if (verifiedMates.mapqL >= options.minMapq && store.alignQualityStore[verifiedMates.alignedReadIdL].errors <= max4Errors && store.alignQualityStore[verifiedMates.alignedReadIdL].score <= options.maxScore)
                        writeBsAlignment(outStream, bamIOContext, store, setContigGaps[verifiedMates.alignedReadIdL], verifiedMates.alignedReadIdL, verifiedMates.mapqL, verifiedMates.recordL, options);
                    max4Errors = floor(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedMates.alignedReadIdR].readId])/100.0);
                    if (verifiedMates.mapqR >= options.minMapq && store.alignQualityStore[verifiedMates.alignedReadIdR].errors <= max4Errors && store.alignQualityStore[verifiedMates.alignedReadIdR].score <= options.maxScore)
                        writeBsAlignment(outStream, bamIOContext, store, setContigGaps[verifiedMates.alignedReadIdR], verifiedMates.alignedReadIdR, verifiedMates.mapqR, verifiedMates.recordR, options);
               }
                else
                {
                    VerifiedRead verifiedRead = verifyRead(store, options);   // Find best alignment and compute mapq, verify
                    unsigned max4Errors = floor(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedRead.alignedReadId].readId])/100.0);
                    if (verifiedRead.mapq >= options.minMapq && store.alignQualityStore[verifiedRead.alignedReadId].errors <= max4Errors && store.alignQualityStore[verifiedRead.alignedReadId].score <= options.maxScore)
                    {
                        if (store.readNameStore[verifiedRead.alignedReadId] == "read.000000860")
                        {
                            std::cout << "read.000000860" << std::endl;
                            std::cout << "max4Error: " << options.max4Error << std::endl;
                            std::cout << "max4Errors: " << max4Errors << std::endl;
                            std::cout << "errors: " << store.alignQualityStore[verifiedRead.alignedReadId].errors << std::endl;
                        }

                        writeBsAlignment(outStream, bamIOContext, store, setContigGaps[verifiedRead.alignedReadId], verifiedRead.alignedReadId, verifiedRead.mapq, verifiedRead.record, options);
                    }
               }
            }
            ////////////////////////////////////////////////////////////////////////////////////
            // Read entry for next read
            clear(store.alignedReadStore);
            clear(store.alignQualityStore);
            clear(matchMateInfos);
            clear(setContigGaps);

            currReadName = record.qName;
        }
       
        if ((record.flag & 0x4) == 1) continue;     // Read is unmapped
        // Get readId (could be curr. read or mate) -> Only Id, seq. we will get the original from readSeqStore 
        // If read name not found, skip entry
        if (!getIdByName(store.readNameStore, record.qName, readId, readNameCache)) continue;

        if ((record.flag & 0x1) == 1)    // If paired: Get readId for current mate
        {
            // if the read is in the store and paired
            // check the mate pair store if it is the same mate of the pair
            // assuming that only one flag 0x040 or 0x0080 is 1
            int inPair = 1 - ((record.flag & 0x40) >> 6);	// bit 7 is set => inPair = 0
                                                             // else inPair = 1 (even if bits 6 and 7 are not set)
            TId matePairId = store.readStore[readId].matePairId;
            if (matePairId != TMatePairStoreElement::INVALID_ID)
            {
                readId = store.matePairStore[matePairId].readId[inPair];
                if (readId == TMatePairStoreElement::INVALID_ID) continue; 
            }
        }
        // Get positions
        unsigned len;
        _getLengthInRef(record.cigar, len);         // 

        TContigPos beginPos = record.beginPos;
        TContigPos endPos = beginPos + len;
        if ( ((record.flag & 0x10) >> 4) == 1) // Reverse 
        {
            TContigPos temp = beginPos;
            beginPos = endPos;
            endPos = temp;
        }
        // Create alignedReadStore entry
        TId id = length(store.alignedReadStore);
        TReadGapAnchors readGapAnchors;
        TAlignedRead alignedRead = TAlignedRead(id, readId, record.rID, beginPos, endPos, readGapAnchors);
        appendValue(store.alignedReadStore, alignedRead);

        TReadSeq readSeq =  store.readSeqStore[readId];
        TContigSeq contigInf;
        TReadGaps readGaps;   
        // Set readGaps source and get contig infix
        TContigPos beginInf;
        TContigPos endInf;
        if (beginPos < endPos) 
        {
            setSource(readGaps, readSeq);
            if (beginPos < options.intervalOffset) beginInf = 0;    
            else beginInf = beginPos-options.intervalOffset;
            if (endPos > (long)length(store.contigStore[record.rID].seq) - options.intervalOffset) endInf = length(store.contigStore[record.rID].seq);
            else endInf = endPos+options.intervalOffset;
        }
        else 
        {
            TReadSeq readSeq = store.readSeqStore[readId];
            reverseComplement(readSeq);
            //ModifiedString<ModifiedString<TReadSeq, ModView<FunctorComplement<Dna5> > >, ModReverse> readSeqRC(readSeq);
            for (int i = length(store.readSeqStore[readId]) - 1; i >= 0; --i)   // assign rev. compl qualities 
            {
                assignQualityValue(readSeq[i], getQualityValue(store.readSeqStore[readId][length(store.readSeqStore[readId]) - 1 - i]));
            }
            assignSource(readGaps, readSeq);  
            if (endPos < options.intervalOffset) beginInf = 0;
            else beginInf = endPos-options.intervalOffset;
            if (beginPos > (long)length(store.contigStore[record.rID].seq) - options.intervalOffset) endInf = length(store.contigStore[record.rID].seq);
            else endInf = beginPos+options.intervalOffset;
        }
        // Create contigGaps
        TContigGaps contigGaps; 
        contigInf = infix(store.contigStore[record.rID].seq, beginInf, endInf);  
        assignSource(contigGaps, contigInf);        // TODO : copy necessary? 
        // Realign and set alignedReadStore entries
        reAlign4(readGaps, contigGaps, store, id, scoringSchemeCTLeft, scoringSchemeCTRight, scoringSchemeGALeft, scoringSchemeGARight, options);
        //std::cout << "align55657:" << std::endl;
        //std::cout << "Contig: " <<  contigGaps.data_cutBegin << "  " << contigGaps.data_cutEnd << "  " << contigGaps.data_viewCutBegin << "  " << contigGaps.data_viewCutEnd << std::endl;
        //std::cout << contigGaps << std::endl;
        //std::cout << readGaps << std::endl;
        //appendValue(setContigGapAnchors, _dataAnchors(contigGaps));
        appendValue(setContigGaps, contigGaps, Generous());

        // Store information about mate
        if (getMateNo(store, readId) == 0 && (record.flag & 0x8) == 0)  // store info only if read is the first mate and if the second mate is mapped
        {
            typename TMatePairStoreElement::TId matePairId = store.readStore[readId].matePairId;
            
            TMatchMateInfo matchMateInfo = {readId, record.rID, id, matePairId, (record.flag & 0x20) != 0, record.pNext};
            appendValue(matchMateInfos, matchMateInfo);
            back(store.alignedReadStore).pairMatchId = id;  // pairMatchId == alignedRead id of first mate // abuse 
        }
    }
    if(!empty(store.alignedReadStore))
    {
        // Deal with last read ...
        if (!empty(store.matePairStore))
        { 
            // set the match mate IDs using the information stored in matchMateInfos
            _generatePseudoPairMatchIds(store, matchMateInfos, options);
            VerifiedMates verifiedMates = verifyMates(store, options);  // For mates: Find best alignments and compute mapq, verify
            unsigned max4Errors = floor(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedMates.alignedReadIdL].readId])/100.0);
            if (verifiedMates.mapqL >= options.minMapq && store.alignQualityStore[verifiedMates.alignedReadIdL].errors <= max4Errors && store.alignQualityStore[verifiedMates.alignedReadIdL].score <= options.maxScore)
                writeBsAlignment(outStream, bamIOContext, store, setContigGaps[verifiedMates.alignedReadIdL], verifiedMates.alignedReadIdL, verifiedMates.mapqL, verifiedMates.recordL, options);
            max4Errors = floor(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedMates.alignedReadIdR].readId])/100.0);
            if (verifiedMates.mapqR >= options.minMapq && store.alignQualityStore[verifiedMates.alignedReadIdR].errors <= max4Errors && store.alignQualityStore[verifiedMates.alignedReadIdR].score <= options.maxScore)
                writeBsAlignment(outStream, bamIOContext, store, setContigGaps[verifiedMates.alignedReadIdR], verifiedMates.alignedReadIdR, verifiedMates.mapqR, verifiedMates.recordR, options);

        }
        else
        {
            VerifiedRead verifiedRead = verifyRead(store, options);   // Find best alignment and compute mapq, verify
            unsigned max4Errors = floor(options.max4Error*length(store.readSeqStore[store.alignedReadStore[verifiedRead.alignedReadId].readId])/100.0);
            if (verifiedRead.mapq >= options.minMapq && store.alignQualityStore[verifiedRead.alignedReadId].errors <= max4Errors && store.alignQualityStore[verifiedRead.alignedReadId].score <= options.maxScore)
                writeBsAlignment(outStream, bamIOContext, store, setContigGaps[verifiedRead.alignedReadId], verifiedRead.alignedReadId, verifiedRead.mapq, verifiedRead.record, options);
        }
    }
    return 0;
}

#endif

