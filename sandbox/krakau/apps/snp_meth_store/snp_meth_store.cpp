/*==========================================================================
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your options) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

//#define SEQAN_ENABLE_PARALLELISM
  
#define CALL_PROFILE
//#define SNPSTORE_DEBUG

#include <seqan/platform.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <seqan/parallel.h>

#ifdef PLATFORM_WINDOWS
#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
#define SEQAN_DEFAULT_TMPDIR "./"
#endif

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <seqan/arg_parse.h>

#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp>

#include "bs_alphabets.h"
#include "util.h"
#include "bs_realignment.h"
#include "snp_meth_store.h"
#include "bs_one_calling.h"

using namespace std;
using namespace seqan;


// load entire genome into memory
template <typename TGenomeSet, typename TGenomeNames>
bool loadGenomes(TGenomeSet &genomes, StringSet<CharString> &fileNameList, ::std::map<CharString,unsigned> &gIdStringToIdNumMap, TGenomeNames & genomeNames)
{
    unsigned gSeqNo = 0;
    unsigned filecount = 0;
    CharString temp;
    clear(genomeNames);
    while(filecount < length(fileNameList))
    {
        clear(temp);
        MultiFasta multiFasta;
        if (!open(multiFasta.concat, toCString(fileNameList[filecount]), OPEN_RDONLY)) return false;
        split(multiFasta, Fasta());
        
        unsigned seqCount = length(multiFasta);
        if(length(genomes) < gSeqNo+seqCount) 
            resize(genomes,gSeqNo+seqCount);
        for(unsigned i = 0; i < seqCount; ++i)
        {
            assignSeq(genomes[gSeqNo+i], multiFasta[i], Fasta());       // read Genome sequence
            assignSeqId(temp, multiFasta[i], Fasta());
            for (unsigned pos = 0; pos < length(temp); ++pos)
            {
                if(temp[pos]=='\t' || temp[pos]=='\b' || temp[pos]==' ')
                {
                    resize(temp,pos);
                    break;
                }
            }
            gIdStringToIdNumMap.insert(::std::make_pair(temp,gSeqNo+i)); // keeps the whole fasta ID including white spaces
            appendValue(genomeNames,temp);
        }
        gSeqNo += seqCount;
        ++filecount;
    }
    resize(genomes,gSeqNo);
    return (gSeqNo > 0);
}





// transform global cooridnates to coordinates relative to chromosomal segment
template<typename TFragmentStore, typename TContigPos, typename TOptions>
void 
transformCoordinates(TFragmentStore &fragmentStore, TContigPos startCoord, TOptions&)
{
    typedef typename TFragmentStore::TAlignedReadStore          TMatches;
    //typedef typename Value<TMatches>::Type                      TMatch;
    typedef typename Iterator<TMatches,Standard>::Type          TMatchIt;
    
    TMatchIt mIt        = begin(fragmentStore.alignedReadStore,Standard());
    TMatchIt mItEnd     = end(fragmentStore.alignedReadStore,Standard());
    
    while(mIt != mItEnd)
    {
        (*mIt).endPos -= startCoord;
        (*mIt).beginPos -= startCoord;
        ++mIt;
    }
    
}

// copy matches relevant for next window
template<typename TFragmentStore, typename TSetContigAnchorGaps, typename TSize, typename TContigPos, typename TOptions>
void 
copyNextWindowMatchesAndReads(TFragmentStore &fragmentStore,
                              TSetContigAnchorGaps &setContigAnchorGaps,
                              typename TFragmentStore::TReadSeqStore &tmpReads,
                              typename TFragmentStore::TReadStore &tmpRs,
                              typename TFragmentStore::TAlignedReadStore &tmpMatches,
                              typename TFragmentStore::TAlignQualityStore &tmpQualities,
                              TSetContigAnchorGaps &tmpSetContigAnchorGaps,
                              TSize ,
                              TContigPos currentWindowEnd,
                              TOptions &options)
{

    typedef typename TFragmentStore::TAlignedReadStore          TMatches;
    typedef typename Value<TMatches>::Type                      TMatch;
    typedef typename Iterator<TMatches,Standard>::Type          TMatchIt;
    typedef typename Iterator<TSetContigAnchorGaps,Standard>::Type          TSetContigGapsIter;
    typedef typename Id<TFragmentStore>::Type                   TId;
    //typedef typename Value<TReadClips>::Type                    TPair;

    SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore),length(fragmentStore.alignQualityStore));

    ::std::sort(begin(fragmentStore.alignedReadStore, Standard()), end(fragmentStore.alignedReadStore, Standard()), LessGPos<TMatch>());    

    if(options._debugLevel > 1 )::std::cout << "Copying matches overlapping more than one window ... \n";
    
    TMatchIt mIt        = end(fragmentStore.alignedReadStore,Standard());
    TMatchIt mItBegin   = begin(fragmentStore.alignedReadStore,Standard());
    --mIt;

    TSetContigGapsIter itG = end(setContigAnchorGaps);
    //TSetContigGapsIter itGBegin = begin(setContigAnchorGaps);
    --itG;

    // We will use minCoord/maxCoord to store the temporarily minimal and maximal coordinates in the window.
    int minCoord = maxValue<int>();
    int maxCoord = minValue<int>();
    //CharString str = "discBef";
    //_dumpMatches(fragmentStore, str);
    
    //std::cout << " max hit length = " << options.maxHitLength << std::endl;
    // look for matches that are inside our window of interest, copy corresponding matches,reads,qualities
    while(mIt >= mItBegin && _min((*mIt).beginPos,(*mIt).endPos) + (TContigPos)options.maxHitLength + (TContigPos)options.windowBuff >= currentWindowEnd )
    {
        if( _max((*mIt).beginPos,(*mIt).endPos) + (TContigPos)options.windowBuff > currentWindowEnd )
        {
            TId id = length(tmpMatches);
            appendValue(tmpMatches,*mIt);
            tmpMatches[id].id = id;
            tmpMatches[id].readId = id;
            appendValue(tmpReads,fragmentStore.readSeqStore[(*mIt).readId]);
            appendValue(tmpRs,fragmentStore.readStore[(*mIt).readId]);
            appendValue(tmpQualities,fragmentStore.alignQualityStore[(*mIt).readId]);
            appendValue(tmpSetContigAnchorGaps, *itG);
            maxCoord = std::max(maxCoord, (int)std::max(mIt->beginPos, mIt->endPos));
            minCoord = std::min(minCoord, (int)std::min(mIt->beginPos, mIt->endPos));
        }
        --mIt;
        --itG;
    }

    // Write minimal and maximal coordinate from reads in this window to options.minCoord/options.maxCoord.
    if (minCoord != maxValue<int>())
        options.minCoord = minCoord;
    if (maxCoord != minValue<int>())
        options.maxCoord = maxCoord;
   
    if(options._debugLevel > 1)
        std::cout << length(tmpMatches)<<" matches left over from previous window.\n"; 
}


// little helper
template<typename TMatch>
char
orientation(TMatch & match)
{
    if(match.endPos > match.beginPos) 
        return 'F';
    else
        return 'R';
}

// discard reads stacking up, give preference to high quality reads
template<typename TFragmentStore, typename TSize, typename TOptions>
void
applyPileupCorrection(TFragmentStore    &fragmentStore, 
                      TSize                         arrayBeginPos, 
                      TSize                         arrayEndPos, 
                      TOptions                      &options)
{
    typedef StringSet<String<Dna5Q>, Owner<ConcatDirect<> > >    TFalseType;
    typedef typename TFragmentStore::TAlignedReadStore  TMatches;
    //typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore TMatchQualities;
    //typedef typename Value<TMatchQualities>::Type       TMatchQuality;
    typedef typename TFragmentStore::TReadSeqStore      TReads;
    //typedef typename Value<TReads>::Type                TRead;
    //typedef typename TFragmentStore::TContigStore       TGenomeSet;
    //typedef typename Value<TGenomeSet>::Type            TGenome;
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename Iterator<TMatches,Standard>::Type  TMatchIterator;
    
    if(IsSameType<TReads, TFalseType >::VALUE)
        std::cout << "Hier stimmt was nciht. strinsetspec concat direct\n";
    
    if(options._debugLevel > 0) std::cout << arrayEndPos - arrayBeginPos  << " matches subject to pile up correction." << std::endl;

    if(options.orientationAware) 
        ::std::sort(iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard()),
                    iter(fragmentStore.alignedReadStore, arrayEndPos, Standard()), 
                    LessGStackOaMQ<TMatches,TMatchQualities>(fragmentStore.alignQualityStore)); 
    else
        ::std::sort(iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard()),
                    iter(fragmentStore.alignedReadStore, arrayEndPos, Standard()), 
                    LessGStackMQ<TMatches,TMatchQualities>(fragmentStore.alignQualityStore));
    
    
    TMatchIterator matchIt          = iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
    TMatchIterator matchRangeEnd    = iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());
    TMatchIterator matchItKeep      = matchIt;
    
    while(matchIt != matchRangeEnd)
    {
        TContigPos currentBegin = _min((*matchIt).beginPos,(*matchIt).endPos);
        TContigPos currentEnd = _max((*matchIt).beginPos,(*matchIt).endPos);
        unsigned currentSeqno = (*matchIt).contigId;
        char currentOrientation = orientation(*matchIt);
        unsigned currPile = 0;
        while(matchIt != matchRangeEnd 
              && (*matchIt).contigId == currentSeqno 
              && _min((*matchIt).beginPos,(*matchIt).endPos) == currentBegin 
              && _max((*matchIt).beginPos,(*matchIt).endPos) == currentEnd 
              && (!options.orientationAware || orientation(*matchIt) == currentOrientation) 
              && currPile < options.maxPile)
        {
            *matchItKeep = *matchIt;
            ++currPile;
            ++matchIt;
            ++matchItKeep;
        }
        //if(matchRangeEnd > matchItEnd) ::std::cerr <<"neeeeeeee\n";
        while(matchIt != matchRangeEnd 
              && (*matchIt).contigId == currentSeqno 
              && _min((*matchIt).beginPos,(*matchIt).endPos) == currentBegin 
              && _max((*matchIt).beginPos,(*matchIt).endPos) == currentEnd 
              && (!options.orientationAware || orientation(*matchIt) == currentOrientation))
            ++matchIt;
        
    }
    
    if(options._debugLevel > 0) std::cout << matchIt - matchItKeep << " matches discarded." << std::endl;
    resize(fragmentStore.alignedReadStore,matchItKeep - begin(fragmentStore.alignedReadStore,Standard()));
    
    //  ::std::cout << "numMatches = " << length(fragmentStore.alignedReadStore) << ::std::endl;
    SEQAN_ASSERT_LEQ(length(fragmentStore.alignedReadStore), length(fragmentStore.alignQualityStore));
    SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore), length(fragmentStore.alignQualityStore));
    
    //  str="pileAft";
    //  _dumpMatches(fragmentStore,str);
}

// average quality of read is kept as extra info for each match
// used for prioritization in pile up correction
template<typename TFragmentStore, typename TSize, typename TOptions>
void
addReadQualityToMatches(TFragmentStore  &fragmentStore, 
                        TSize                           arrayBeginPos,
                        TSize                           arrayEndPos,
                        TOptions &)
{
    typedef typename TFragmentStore::TAlignedReadStore      TMatches;
    typedef typename Value<TMatches>::Type                  TMatch;
    typedef typename TFragmentStore::TReadSeqStore          TReads;
    typedef typename Value<TReads>::Type                    TRead;
    typedef typename Iterator<TMatches, Standard>::Type     TIterator;
    
    TIterator it = iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
    TIterator itEnd = iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());
    
    int avgRQ;
    for (; it != itEnd; ++it) 
    {
        TRead const &read = fragmentStore.readSeqStore[(*it).readId];
        avgRQ = 0;
        for(unsigned i = 0; i < length(read); ++i)
            avgRQ += (int) getQualityValue(read[i]);
        // watch out, this is new: use mapping quality if given
        //if((fragmentStore.alignQualityStore[(*it).id]).score == 0 || (char)(avgRQ/length(read))<(fragmentStore.alignQualityStore[(*it).id]).score)
        if((char)(avgRQ/length(read))<(fragmentStore.alignQualityStore[(*it).id]).score)

            (fragmentStore.alignQualityStore[(*it).id]).score = (char)(avgRQ/length(read));
    } 
}

// checks whether an alignment has indels
template<typename TValue, typename TSpec>
bool alignmentHasIndel(Align<TValue,TSpec> &align)
{
    typedef Align<TValue,TSpec>  TAlign;
    typedef typename Row<TAlign>::Type      TAlignRow;
    typedef typename Iterator<TAlignRow>::Type  TAlignIterator;
    
    bool hasIndel = false;
    for(unsigned i = 0; !hasIndel && i < length(rows(align)); ++i)
    {
        TAlignIterator it = begin(row(align,i));
        TAlignIterator itEnd = end(row(align,i));
        while (it != itEnd)
        {
            if(isGap(it))
            {
                hasIndel = true;
                break;
            }
            ++it;
        }
    }
    return hasIndel;
}

// perform read clipping
template<typename TFragmentStore, typename TReadClips, typename TSize, typename TOptions>
void
clipReads(TFragmentStore    &fragmentStore,
          TReadClips    &readClips,
          TSize     arrayBeginPos,
          TSize     arrayEndPos,
          TOptions  &options)
{
    typedef typename TFragmentStore::TAlignedReadStore      TMatches;
    typedef typename Value<TMatches>::Type              TMatch;
    typedef typename TFragmentStore::TAlignQualityStore     TAlignQualityStore;     // TMatchQualities
    typedef typename Value<TAlignQualityStore>::Type        TAlignQuality;
    typedef typename TFragmentStore::TReadSeqStore          TReads;
    typedef typename Value<TReads>::Type                TRead;
    typedef typename Iterator<TMatches, Standard>::Type     TIterator;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;             // TGenome
    typedef typename Position<TContigSeq>::Type         TContigPos;
    
    TIterator it = iter(fragmentStore.alignedReadStore, arrayBeginPos, Standard());
    TIterator itEnd = iter(fragmentStore.alignedReadStore, arrayEndPos, Standard());
    Align<TRead, ArrayGaps> align;
    resize(rows(align), 2);
#ifdef SNPSTORE_DEBUG
    bool extraV = true;
#endif
    
    Score<int> scoreType = Score<int>(0, -999, -1001, -1000);
    if(length(readClips) < (arrayEndPos-arrayBeginPos)) ::std::cout << length(readClips) << " readclips but " << (arrayEndPos-arrayBeginPos) << " many reads.\n";
    TContigSeq &currGenome = fragmentStore.contigStore[0].seq;
    int kickout = 0;
    for (; it != itEnd; ++it) 
    {
        TMatch &m = *it;
        TRead &read = fragmentStore.readSeqStore[m.readId];
        int clipLeft = readClips[m.readId].i1;
        int clipRight = readClips[m.readId].i2;
        TContigPos beginPos = (m.beginPos < m.endPos ) ? m.beginPos : m.endPos;
        TContigPos endPos = (m.beginPos > m.endPos ) ? m.beginPos : m.endPos;
        TAlignQuality &aliQ = fragmentStore.alignQualityStore[m.id];
        
#ifdef SNPSTORE_DEBUG
        TContigPos beginBefore = beginPos;
#endif
        if(m.id != m.readId) ::std::cout << "match id != readId \n";
        if(clipLeft+clipRight > (int)length(read) || clipLeft > (int)length(read) || clipRight > (int)length(read))
        {
            if(options._debugLevel > 0)::std::cout << "WARNING: clipLeft+clipRight > readLen!!\n";
#ifdef SNPSTORE_DEBUG
            ::std::cout << "readlength = "<<length(read)<< " \n";
            ::std::cout << "readId = "<<m.readId << "id=" << m.id << " \n";
            ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
            ::std::cout << "read=" << read << std::endl;
            ::std::cout << "beginPos=" << beginPos << std::endl;
            ::std::cout << "endPos=" << endPos << std::endl; 
#endif
            clipLeft = length(read);
            clipRight = 0;
        }
#ifdef SNPSTORE_DEBUG
            ::std::cout << "readlength = "<<length(read)<< " \n";
            ::std::cout << "readId = "<<m.readId << "id=" << m.id << " \n";
            ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
            ::std::cout << "read=" << read << std::endl;
            ::std::cout << "beginPos=" << beginPos << std::endl;
            ::std::cout << "endPos=" << endPos << std::endl; 
#endif
        if(clipLeft > 0 || clipRight > 0)
        {
            //  if(extraV) ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << std::endl;
            if((int)length(read)-clipLeft-clipRight < options.minClippedLength)
            {
                if(options._debugLevel > 1 )
                    ::std::cout << "Discarded: "<<read<<" at position "<< beginPos <<"\n";
                m.endPos = m.beginPos = 0; // "mask" read
                read = "";
                ++kickout;
                continue;
            }
            // adjust read sequence
            read = infix(read,clipLeft,length(read)-clipRight);
            
            // upate avg read quality
            int avgRQ = 0;
            for(unsigned i = 0; i < length(read); ++i)
                avgRQ += (int) getQualityValue(read[i]);
            aliQ.score = (char)(avgRQ/length(read));
            //      if(extraV) ::std::cout << "aliQ.score = " << (int)aliQ.score << ::std::endl;
            
            
            //do semi-global alignment
            assignSource(row(align, 0), read);
            assignSource(row(align, 1), infix(currGenome, beginPos, endPos));
            if ((*it).endPos < (*it).beginPos)
                reverseComplement(source(row(align, 1)));
            
            int score = globalAlignment(align, scoreType, AlignConfig<false,true,true,false>(), Gotoh());       
            aliQ.errors = (unsigned) round((float)-score/1000);
            
#ifdef SNPSTORE_DEBUG
            if(extraV) ::std::cout << align << std::endl;
            if(extraV) ::std::cout << "aliQ.errors = " << (int) aliQ.errors << std::endl; 
#endif
            
            // transform first and last read character to genomic positions
            if(aliQ.pairScore == 1)
            {
                unsigned viewPosReadFirst = toViewPosition(row(align, 0), 0);
                unsigned viewPosReadLast  = toViewPosition(row(align, 0), length(read) - 1);
                unsigned genomePosReadFirst = toSourcePosition(row(align, 1), viewPosReadFirst);
                unsigned genomePosReadLast  = toSourcePosition(row(align, 1), viewPosReadLast);
                //              if(isGap(row(align,1),viewPosReadFirst))
                //              {
                //                  std::cout << "bgein gap --> do nothing " << std::endl;
                //                  
                //              }
                
                if(isGap(row(align,1),viewPosReadLast))
                {
                    genomePosReadLast--;                    
                }
#ifdef SNPSTORE_DEBUG
                if(extraV)
                {   ::std::cout << "viewPosReadFirst = " << viewPosReadFirst << ::std::endl;
                    ::std::cout << "viewPosReadLast = " << viewPosReadLast << ::std::endl;
                    ::std::cout << "genomePosReadFirst = " << genomePosReadFirst << ::std::endl;
                }
#endif
                if(m.beginPos < m.endPos) // forward
                {
                    endPos = beginPos + (genomePosReadLast + 1);
                    beginPos += genomePosReadFirst;
                }
                else
                {
                    beginPos = endPos - genomePosReadLast - 1;
                    endPos = endPos - genomePosReadFirst;
                }
                
                if(!alignmentHasIndel(align)) aliQ.pairScore = 0;
            }
            else
            {
                if(m.beginPos < m.endPos) // forward
                {
                    endPos -= clipRight;
                    beginPos += clipLeft;
                }
                else
                {
                    endPos -= clipLeft;
                    beginPos += clipRight;
                }
            }
            
            // set genomic positions
            if(m.beginPos < m.endPos) // forward
            {
                m.endPos = endPos;
                m.beginPos = beginPos;
            }
            else
            {
                m.endPos = beginPos;
                m.beginPos = endPos;
            }
            if(beginPos > 300000000 || endPos > 300000000)
            {
                ::std::cout << "bgein groesser 300mio neu beginPos = "<< beginPos << " endpos=" << endPos << ::std::endl;
#ifdef SNPSTORE_DEBUG
                ::std::cout << "WARNING: clipLeft+clipRight > readLen!!??\n";
                ::std::cout << "readId = "<<m.readId << "id=" << m.id << " \n";
                ::std::cout << "clipLeft = " << clipLeft << " clipRight = "<<clipRight << "\n";
                ::std::cout << "read=" << read << std::endl;
                ::std::cout << "beginPos before=" << beginBefore << std::endl;
                ::std::cout << "beginPos=" << beginPos << std::endl;
                ::std::cout << "endPos=" << endPos << std::endl;
#endif
            }
            
            
        }
    }
    if(options._debugLevel > 0) 
        ::std::cout << kickout <<" reads too short after clipping, discarding!\n";
    
}

template<typename TTmpReads, typename TTmpMatches, typename TTmpQualities, typename TStr>
void
_dumpMatches(TTmpReads & reads, TTmpMatches & matches, TTmpQualities & qualities, TStr str)
{
    
    std::cout << "Length of matches = " << length(matches)  << "\n";
    std::cout << "Length of reads   = " << length(reads)  << "\n";
    std::cout << "Length of matchqs = " << length(qualities)  << "\n";
    
    for(unsigned i = 0 ; i < length(matches); ++i)
    {
        char ori = (matches[i].beginPos < matches[i].endPos) ? 'F' : 'R';
        std::cout << "--"<<str<<"Match number " << i << ":\n";
        std::cout << "--"<<str<<"MatchId  = " << matches[i].id << "\n";
        std::cout << "--"<<str<<"ReadId   = " << matches[i].readId << "\n";
        std::cout << "--"<<str<<"ContigId = " << matches[i].contigId << std::flush << "\n";
        std::cout << "--"<<str<<"gBegin   = " << _min(matches[i].beginPos, matches[i].endPos) << "\n";
        std::cout << "--"<<str<<"gEnd     = " << _max(matches[i].beginPos, matches[i].endPos) << "\n";
        std::cout << "--"<<str<<"orient   = " << ori << std::flush << std::endl;
        if(length(qualities) > matches[i].id)
        {
            std::cout << "--"<<str<<"EditDist = " << (int) qualities[matches[i].id].errors << "\n";
            std::cout << "--"<<str<<"AvgQ     = " << (int) qualities[matches[i].id].score << "\n";
        }
        std::cout << "--"<<str<<"Readseq  = " << reads[matches[i].readId] << std::flush << "\n";
        
    }
}

template<typename TContigIntervals, typename TFragmentStore, typename TMethOptions>
inline void
assignIntervalsToContigs(TContigIntervals &contigIntervals, TFragmentStore &fragStore, TMethOptions &methOptions)
{
    typedef typename Value<TContigIntervals>::Type      TIntervals;
    typedef typename Value<TIntervals>::Type            TInterval;

    resize(contigIntervals, length(fragStore.contigStore));
    unsigned i = 0;
    for (unsigned contigId = 0; contigId < length(fragStore.contigStore); ++contigId)
    {
        while (i < length(methOptions.intervals) && methOptions.intervals[i].contigName == fragStore.contigNameStore[contigId])
        {
            TInterval interval;
            interval.i1 = methOptions.intervals[i].startPos;
            interval.i2 = methOptions.intervals[i].endPos;
            appendValue(contigIntervals[contigId], interval);
            ++i;
        }
    }
}

template<typename TFileStream, typename TSpec, typename TContigId, typename TReaders, typename TContexts, typename TRecords, typename TContigIntervals, typename TOptions, typename TMethOptions>
inline bool
detectSNPsForContig(TFileStream &fileStream, 
                    FragmentStore<TSpec> &fragmentStore1, 
                    TContigId &currContigId, 
                    TReaders &recordReaders, 
                    TContexts &contexts, 
                    TRecords &records, 
                    TContigIntervals &contigIntervals,  
                    TOptions &options, 
                    TMethOptions &methOptions)
{
    typedef          FragmentStore<TSpec>               TFragmentStore;
    typedef typename TFragmentStore::TContigSeq         TContigSeq;             // TGenome
    typedef typename TFragmentStore::TContigPos         TContigPos;
    typedef typename TFragmentStore::TAlignedReadStore  TAlignedReadStore;      // TMatches
    typedef typename Value<TAlignedReadStore>::Type     TAlignedRead;
    typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;     // TMatchQualities
    typedef typename TFragmentStore::TReadStore         TReadStore;             // TReadSet
    typedef typename TFragmentStore::TReadSeqStore      TReadSeqStore;              // TReadSet
    typedef typename Value<TReadSeqStore>::Type         TReadSeq;
    typedef typename TFragmentStore::TContigStore       TContigStore;           // TGenomeSet
    typedef typename Value<TContigStore>::Type          TContig;
    typedef          TContigSeq                         TGenome;
    typedef          StringSet<TGenome>                 TGenomeSet;

    typedef String<typename TFragmentStore::TContigGapAnchor>               TContigAnchorGaps;
    typedef Gaps<Dna5String, AnchorGaps<TContigAnchorGaps> >                TContigGaps;
    typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> > TReadGaps;

    typedef typename TFragmentStore::TContigNameStore                       TContigNameStore;
    typedef NameStoreCache<TContigNameStore, CharString>                    TContigNameStoreCache;
    typedef BamIOContext<TContigNameStore, TContigNameStoreCache>           TBamIOContext;
    typedef String<String<typename TFragmentStore::TContigGapAnchor> >      TSetContigAnchorGaps;
    
    if (!empty(methOptions.intervals) && empty(contigIntervals[currContigId])) return 0;
    unsigned currInterval = 0;

    // parse matches batch by batch
    TContigPos currentWindowBegin = 0;
    if(options._debugLevel > 0) ::std::cout << "Scanning genome #" << currContigId << " ..." << ::std::endl;
    
    // containers for those matches that overlap two windows    
    TAlignedReadStore tmpMatches;
    TAlignQualityStore tmpQualities;
    TReadStore tmpRs;
    TReadSeqStore tmpReads;     // Something went wrong when keeping all reads of contig, so keep it like that for the moment
    TSetContigAnchorGaps tmpSetContigAnchorGaps;
    options.minCoord = MaxValue<unsigned>::VALUE;
    options.maxCoord = 0;
    
    // snp calling is done for all positions between windowBegin and windowEnd
    // bs_change
    if (!empty(methOptions.intervals))
    {
        SEQAN_ASSERT_LT((unsigned)contigIntervals[currContigId][currInterval].i2, (unsigned)length(fragmentStore1.contigStore[currContigId].seq));
        currentWindowBegin = contigIntervals[currContigId][currInterval].i1;
    }
    while(currentWindowBegin <  (TContigPos)length(fragmentStore1.contigStore[currContigId].seq))
    {
        TContigPos currentWindowEnd = currentWindowBegin + options.windowSize; 
        if(currentWindowEnd > (TContigPos)length(fragmentStore1.contigStore[currContigId].seq)) currentWindowEnd = (TContigPos)length(fragmentStore1.contigStore[currContigId].seq);
        // bs_change
        if (!empty(methOptions.intervals) && currentWindowEnd > contigIntervals[currContigId][currInterval].i2) 
        {
            currentWindowEnd = contigIntervals[currContigId][currInterval].i2;
            //std::cout << "currInterval startPos: " << contigIntervals[currContigId][currInterval].i1 << "  endPos: " << contigIntervals[currContigId][currInterval].i2 << std::endl;
        }
        //std::cout << "currentWindowBegin: " << currentWindowBegin << "  currentWindowEnd: " << currentWindowEnd << std::endl;

        if(options._debugLevel > 0)
            std::cout << "Sequence number " << currContigId << " window " << currentWindowBegin << ".." << currentWindowEnd << "\n";
       
        TFragmentStore fragmentStoreTmp;  // Use temp. fragmentStore to store only current reads and to use only infix of contig seq for realigning etc. 
        TSetContigAnchorGaps setContigAnchorGaps;
        
        // add the matches that were overlapping with this and the last window (copied in order to avoid 2 x makeGlobal)
        if(!empty(tmpMatches))
        {
            //sumreads -=  length(tmpReads);  // count these reads only once
            resize(fragmentStoreTmp.alignQualityStore,length(tmpMatches));
            resize(fragmentStoreTmp.alignedReadStore,length(tmpMatches));
            resize(fragmentStoreTmp.readSeqStore,length(tmpMatches));
            resize(fragmentStoreTmp.readStore,length(tmpMatches));
            resize(setContigAnchorGaps, length(tmpMatches));
            
            arrayMoveForward(begin(tmpQualities,Standard()), end(tmpQualities,Standard()), begin(fragmentStoreTmp.alignQualityStore,Standard()));
            arrayMoveForward(begin(tmpMatches,Standard()), end(tmpMatches,Standard()), begin(fragmentStoreTmp.alignedReadStore,Standard()));
            arrayMoveForward(begin(tmpReads,Standard()), end(tmpReads,Standard()), begin(fragmentStoreTmp.readSeqStore,Standard()));
            arrayMoveForward(begin(tmpRs,Standard()), end(tmpRs,Standard()), begin(fragmentStoreTmp.readStore,Standard()));
            arrayMoveForward(begin(tmpSetContigAnchorGaps,Standard()), end(tmpSetContigAnchorGaps,Standard()), begin(setContigAnchorGaps,Standard()));
#ifdef SNPSTORE_DEBUG
            CharString strstr = "afterCopyInFrag";
            _dumpMatches(fragmentStoreTmp,strstr);
#endif  
        }     
        TContig contig;
        contig.seq = fragmentStore1.contigStore[currContigId].seq;
        appendValue(fragmentStoreTmp.contigStore, contig);
        appendValue(fragmentStoreTmp.contigNameStore, fragmentStore1.contigNameStore[currContigId]);

        // parse matches for current window
        if(options._debugLevel > 0) std::cout << "Parsing reads up to position " << currentWindowEnd << "...\n";
        for(unsigned j = 0; j < length(options.readFNames); ++j)
        {
            unsigned sizeBefore = length(fragmentStoreTmp.alignedReadStore);
            
            int result = readMatchesFromSamBam(setContigAnchorGaps, *recordReaders[j], contexts[j], records[j], fragmentStoreTmp, fragmentStore1,
                                               currContigId, currentWindowBegin, currentWindowEnd, options);                       

            if(result == CALLSNPS_GFF_FAILED)
            {
                std::cerr << "Failed to open read file " << options.readFNames[j] << "or reads are not sorted correctly. " << ::std::endl;
                return result;
            }
            if(result > 0) return result;
            if(options._debugLevel > 0) std::cout << "parsed reads of file " << j << "\n";
            
            // store average quality of each read if mapq is 0 or if average quality is < mapq
            // addReadQualityToMatches(fragmentStoreTmp,sizeBefore,(unsigned)length(fragmentStoreTmp.alignedReadStore),options); // bs_change: we will deal with real mapqs
            // do pile up correction if lane-based. lane-specific pileup correction not really used
            if(options.maxPile != 0 && options.laneSpecificMaxPile) 
                applyPileupCorrection(fragmentStoreTmp,sizeBefore,(unsigned)length(fragmentStoreTmp.alignedReadStore),options);
        }
        SEQAN_ASSERT_EQ(length(fragmentStoreTmp.alignedReadStore), length(setContigAnchorGaps));

        if (options._debugLevel > 1)  // number of reads currently in memory
            std::cerr << lengthSum(fragmentStoreTmp.readSeqStore) << " bps of " << length(fragmentStoreTmp.readSeqStore) << " reads in memory." << ::std::endl;
        //sumreads +=  length(fragmentStoreTmp.readSeqStore);  // total count of reads
        
        // do merged pile up correction
        if(options.maxPile != 0 && !options.laneSpecificMaxPile)
            applyPileupCorrection(fragmentStoreTmp,(unsigned)0,(unsigned)length(fragmentStoreTmp.alignedReadStore),options);
        
        // these were set while parsing matches, first and last position of parsed matches
        TContigPos startCoord = _max((int)options.minCoord-options.realignAddBorder,0);// can be < currentWindowBegin
        TContigPos endCoord = _min(options.maxCoord+options.realignAddBorder,length(fragmentStore1.contigStore[currContigId].seq)); // can be > currentWindoEnd 
         
        if(!empty(fragmentStoreTmp.alignedReadStore))
        {
            //initial values of min and max coords for next round are set here
            if(currentWindowEnd != (TContigPos)length(fragmentStore1.contigStore[currContigId].seq))
            {
                clear(tmpMatches);
                clear(tmpQualities);
                clear(tmpRs);
                clear(tmpReads);
                clear(tmpSetContigAnchorGaps);
                copyNextWindowMatchesAndReads(fragmentStoreTmp, setContigAnchorGaps, tmpReads, tmpRs, tmpMatches, tmpQualities, tmpSetContigAnchorGaps, currContigId, currentWindowEnd, options);
            }   
#ifdef SNPSTORE_DEBUG
            std::cout << "Min = " << options.minCoord << " Max = " << options.maxCoord << std::endl;
            std::cout << "startCoord = " << startCoord << " endCoord = " << endCoord << std::endl;
            std::cout << "currentWindowBegin = " << currentWindowBegin << " currentWindowEnd = " << currentWindowEnd << std::endl;
#endif
            // Aligned read coordinates are relative to current chromosomal window (segment)
            transformCoordinates(fragmentStoreTmp,startCoord,options);
            // set the current contig segment as contig sequence 
            fragmentStoreTmp.contigStore[0].seq = infix(fragmentStore1.contigStore[currContigId].seq,startCoord,endCoord);
            
            if(options.realign)
            {
                doCheckRealignCall(fragmentStoreTmp, startCoord, currentWindowBegin, currentWindowEnd, setContigAnchorGaps, fileStream, methOptions, options);
            }
            else
            {
                // Convert gaps of pairwise matches to msa
                //std::cout << "convertPairWiseToGlobalAlignment..." << std::endl;
#ifdef CALL_PROFILE 
                double timeStamp = sysTime();
#endif 
                convertPairWiseToGlobalAlignment(fragmentStoreTmp, setContigAnchorGaps);
#ifdef CALL_PROFILE
                Times::instance().time_convertPWToGlobal += (sysTime() - timeStamp);
#endif
                doSnpAndMethCalling(fragmentStoreTmp,  startCoord, currentWindowBegin, currentWindowEnd, false, fileStream, methOptions, options);    //bs
            }
        }
        if (!empty(methOptions.intervals) && currentWindowEnd == contigIntervals[currContigId][currInterval].i2)   // Reached end off current interval
        {
            if (currInterval < length(contigIntervals[currContigId])-1)   // If existent, go to next interval within the current contig
            {
                currentWindowBegin = contigIntervals[currContigId][currInterval+1].i1;
                clear(tmpReads); 
                clear(tmpRs); 
                clear(tmpMatches); 
                clear(tmpQualities);
                clear(tmpSetContigAnchorGaps);
            }
            else currentWindowBegin = length(fragmentStore1.contigStore[currContigId].seq);    // No interval int this contig to analyze anymore -> jump to end of contig
            ++currInterval;
        }
        else currentWindowBegin = currentWindowEnd;
    } 
    return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec, typename TMethOptions>
int detectSNPs(SNPCallingOptions<TSpec> &options, TMethOptions &methOptions)
{ 
    typedef FragmentStore<SnpStoreSpec_>                                        TFragmentStore;
    typedef typename TFragmentStore::TContigPos                                 TContigPos;
    typedef typename TFragmentStore::TContigNameStore                           TContigNameStore;
    typedef NameStoreCache<TContigNameStore, CharString>                        TContigNameStoreCache;
    typedef BamIOContext<TContigNameStore, TContigNameStoreCache>               TBamIOContext;
    
    SEQAN_PROTIMESTART(load_time);
    
    // Load genomes in FragmentStore
    TFragmentStore      fragmentStore1; 
    StringSet<CharString> genomeFileNameList; // possible for later to load multiple files
    appendValue(genomeFileNameList,options.genomeFName);
    loadContigs(fragmentStore1, genomeFileNameList);
    refresh(fragmentStore1.contigNameStoreCache);
   
    // Assign intervals to analyze to string which can be accessed by contigId
    String<String<Interval<TContigPos> > > contigIntervals;
    if(!empty(methOptions.intervals)) assignIntervalsToContigs(contigIntervals, fragmentStore1, methOptions);
    // Prepare genotype priors
    computeGenotypePriors(methOptions, options);
    // Store fileName for temp file for each contig, needed later to merge in correct order
    String<CharString> contigTempFileNames;
    resize(contigTempFileNames, length(fragmentStore1.contigStore));
    
#ifndef SEQAN_ENABLE_PARALLELISM
    // For each contig we need own record readers (for each given sam file)
    vector< ::std::fstream* > readFileStreams;
    readFileStreams.resize(length(options.readFNames));
    String<RecordReader<std::fstream, SinglePass< > >* > recordReaders;
    String<TBamIOContext > contexts;
    String<BamAlignmentRecord> records; 
    resize(recordReaders, length(options.readFNames));
    resize(contexts, length(options.readFNames));
    resize(records, length(options.readFNames));
    for (unsigned i = 0; i < length(options.readFNames); ++i)
    {
        readFileStreams[i] = new fstream(toCString(options.readFNames[i]), ios_base::in | ios::binary);
        if(!(*(readFileStreams[i])).is_open())
        {
            ::std::cerr << "Failed to open read file " << options.readFNames[i] << ::std::endl;
            return CALLSNPS_GFF_FAILED;
        }
        recordReaders[i] = new RecordReader<std::fstream,SinglePass< > >(*readFileStreams[i]);
        contexts[i] = TBamIOContext(fragmentStore1.contigNameStore, fragmentStore1.contigNameStoreCache);
        // Read header.
        BamHeader header;
        readRecord(header, contexts[i], *recordReaders[i], Sam()); 
        clear(records[i].qName);
    }
#endif

#ifdef SEQAN_ENABLE_PARALLELISM
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) // TODO Check if guided is faster
#endif
    for (unsigned currContigId = 0; currContigId < length(fragmentStore1.contigStore); ++currContigId)
    {
#ifdef SEQAN_ENABLE_PARALLELISM
        // For each contig we need own record readers (for each given sam file)
        vector< ::std::fstream* > readFileStreams;
        readFileStreams.resize(length(options.readFNames));
        String<RecordReader<std::fstream, SinglePass< > >* > recordReaders;
        String<TBamIOContext > contexts;
        String<BamAlignmentRecord> records; 
        resize(recordReaders, length(options.readFNames));
        resize(contexts, length(options.readFNames));
        resize(records, length(options.readFNames));
        for (unsigned i = 0; i < length(options.readFNames); ++i)
        {
            readFileStreams[i] = new fstream(toCString(options.readFNames[i]), ios_base::in | ios::binary);
            if(!(*(readFileStreams[i])).is_open())
            {
                ::std::cerr << "Failed to open read file " << options.readFNames[i] << ::std::endl;
                return CALLSNPS_GFF_FAILED;
            }
            recordReaders[i] = new RecordReader<std::fstream,SinglePass< > >(*readFileStreams[i]);
            contexts[i] = TBamIOContext(fragmentStore1.contigNameStore, fragmentStore1.contigNameStoreCache);
            // Read header.
            BamHeader header;
            readRecord(header, contexts[i], *recordReaders[i], Sam()); 
            clear(records[i].qName);
        }
#endif
        // Temp. contig output
        CharString tempFileName; 
        append(tempFileName, SEQAN_TEMP_FILENAME());      // TODO save in "/dev/shm", in arbeitspeicher, schneller, nur wenn files nicht zu groß
        std::cout << "temp file Name: " << tempFileName << std::endl;
        stringstream ss;
        ss << currContigId;
        append(tempFileName, ss.str());
        contigTempFileNames[currContigId] = tempFileName;
        std::ofstream tempFileStream; 
        tempFileStream.open(toCString(tempFileName), std::ios::binary | std::ios_base::out);
        if(!tempFileStream.is_open())
        {
            std::cerr << "ERROR: Failed to open temp. output file for current contig."  << std::endl;
        }
        detectSNPsForContig(tempFileStream, fragmentStore1, currContigId, recordReaders, contexts, records, contigIntervals, options, methOptions); 
        tempFileStream.close();
    }    

    // open out file streams and store open file pointers
    ::std::ofstream snpFileStream; 
    snpFileStream.open(toCString(options.outputSNP),::std::ios_base::out);
    if(!snpFileStream.is_open())
        return CALLSNPS_OUT_FAILED;
    snpFileStream << "#" << (options.programCall).str() << "\n";
    std::cout << "outputFormat: " << options.outputFormat << std::endl;
    if(options.outputFormat < 2)
    {
        snpFileStream << "#chr\tpos\tref\tcov\tcall";
        snpFileStream << "\tScore";
        snpFileStream << "\n";
    }
    // Append content of temp files to real output in contig order
    CharString buffer;
    resize(buffer, 1000);
    for (unsigned i = 0; i < length(fragmentStore1.contigStore); ++i)
    {
        std::fstream tempFileStream; 
        tempFileStream.open(toCString(contigTempFileNames[i]), std::ios::binary | std::ios_base::in);
        if(!tempFileStream.is_open())
        {
            std::cerr << "ERROR: Failed to open temp. contig output file." << contigTempFileNames[i] << std::endl;
            return CALLSNPS_OUT_FAILED;
        }
        while (!streamEof(tempFileStream) && seqan::streamError(tempFileStream) == 0)
        {
            int num = streamReadBlock(&buffer[0], tempFileStream, length(buffer));
            streamWriteBlock(snpFileStream, &buffer[0], num);
        }
        remove(toCString(contigTempFileNames[i]));  // Delete temp. files
    }
    snpFileStream.close();

    methOptions.statsCGMethylated = methOptions.statsCGMethylated/methOptions.countCG;
    methOptions.statsCHGMethylated = methOptions.statsCHGMethylated/methOptions.countCHG;
    methOptions.statsCHHMethylated = methOptions.statsCHHMethylated/methOptions.countCHH;
    std::cout << "Average methylation rate in CG context: " << methOptions.statsCGMethylated << std::endl;
    std::cout << "Average methylation rate in CHG context: " << methOptions.statsCHGMethylated << std::endl;
    std::cout << "Average methylation rate in CHH context: " << methOptions.statsCHHMethylated << std::endl;
    return 0;
}


template <typename TSpec, typename TMethOptions>
ArgumentParser::ParseResult
parseCommandLine(SNPCallingOptions<TSpec> & options, TMethOptions &methOptions, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("snp_meth_store");
    // Set short description, version, and date.
    setShortDescription(parser, "SnpStore");
    setVersion(parser, "1.0");
    setDate(parser, "October 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIGENOME FILE\\fP\" \"\\fIMAPPED READ FILE(S)\\fP\" ");
    addDescription(parser, "SNP and Indel Calling in Mapped Read Data.");

    // We require two arguments.
    addArgument(parser, ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, 0, "fasta fa FASTA FA");
    addArgument(parser, ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "IN"));
    setValidValues(parser, 1, "sam");

    addSection(parser, "Options: ");
    addOption(parser, ArgParseOption("o", "output", "Output file for SNPs (must be set, no default construction).", ArgParseArgument::OUTPUTFILE));
    setValidValues(parser, "o", "vcf");
    addOption(parser, ArgParseOption("if", "input-format", "Set input format: 0 for GFF format and 1 for SAM format (both must be sorted according to genome positions).", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "if", "0");
    addOption(parser, ArgParseOption("of", "output-format", "Set output format: 0 to output all candidate snps amd 1 to output successful candidate snps only testtesttest.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "of", "0");
    addOption(parser, ArgParseOption("dc", "dont-clip", "Ignore clip tags in gff. Default: off."));

    addOption(parser, ArgParseOption("mu", "multi", "Keep non-unique fragmentStore.alignedReadStore. Default: off."));
    addOption(parser, ArgParseOption("hq", "hide-qualities", "Only show coverage (no qualities) in SNP output file. Default: off."));
    addOption(parser, ArgParseOption("sqo", "solexa-qual-offset", "Base qualities are encoded as Ascii value - 64 (instead of Ascii - 33)."));
    addOption(parser, ArgParseOption("id", "indel-file", "Output file for called indels in gff format. Default: off.", ArgParseArgument::OUTPUTFILE));  
    setValidValues(parser, "id", "vcf");
    addOption(parser, ArgParseOption("m", "method", "Set method used for SNP calling: 0 for threshold method and 1 for maq method.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "m", "1");
    addOption(parser, ArgParseOption("mp", "max-pile", "Maximal number of matches allowed to pile up at the same genome position.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("mmp", "merged-max-pile", "Do pile up correction on merged lanes. Default: off."));
    addOption(parser, ArgParseOption("mc", "min-coverage", "Minimal required number of reads covering a candidate position.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("fc", "force-call", "Always call base if count is >= fc, ignore other parameters. Default: off.", ArgParseArgument::INTEGER));
    setMinValue(parser, "force-call", "1");
    addOption(parser, ArgParseOption("oa", "orientation-aware", "Distinguish between forward and reverse reads. Default: off."));
    addOption(parser, ArgParseOption("cnv", "output-cnv", "Name of CNV result file.", ArgParseArgument::OUTPUTFILE));
    setValidValues(parser, "cnv", "cnv");
    hideOption(parser, "cnv");
    addOption(parser, ArgParseOption("op", "output-positions", "Name of positions output file.", ArgParseArgument::OUTPUTFILE)); 
    setValidValues(parser, "op", "txt");
    hideOption(parser, "op");
    addOption(parser, ArgParseOption("ip", "input-positions", "Name of positions input file.", ArgParseArgument::INPUTFILE)); 
    setValidValues(parser, "ip", "txt");
    hideOption(parser, "ip");
    addOption(parser, ArgParseOption("mpr", "max-polymer-run", "Discard indels in homopolymer runs longer than mpr.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("dp", "diff-pos", "Minimal number of different read positions supporting the mutation.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("eb", "exclude-border", "Exclude read positions within eb base pairs of read borders for SNV calling.", ArgParseArgument::INTEGER)); 
    addOption(parser, ArgParseOption("su", "suboptimal", "Keep suboptimal reads."));  
    addOption(parser, ArgParseOption("re", "realign", "Realign reads around indel candidates."));   
    addOption(parser, ArgParseOption("cq", "corrected-quality", "New quality calibration factor.", ArgParseArgument::DOUBLE));
    hideOption(parser, "cq");
    addOption(parser, ArgParseOption("pws", "parse-window-size", "Genomic window size for parsing reads (concerns memory consumption, choose smaller windows for higher coverage).", ArgParseArgument::INTEGER));
    setMinValue(parser, "parse-window-size", "1");
    setMaxValue(parser, "parse-window-size", "100000");
    addOption(parser, ArgParseOption("reb", "realign-border", "Realign border.", ArgParseArgument::INTEGER));
    setMinValue(parser, "realign-border", "0");
    setMaxValue(parser, "realign-border", "10");
    hideOption(parser, "reb");

    addSection(parser, "SNP calling options: ");
    addSection(parser, " Threshold method related: ");
    addOption(parser, ArgParseOption("mm", "min-mutations", "Minimal number of observed mutations for mutation to be called.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("pt", "perc-threshold", "Minimal percentage of mutational base for mutation to be called.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mq", "min-quality", "Minimal average quality of mutational base for mutation to be called.", ArgParseArgument::DOUBLE));
    addSection(parser, " Maq method related: ");
    addOption(parser, ArgParseOption("th", "theta", "Dependency coefficient.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("hr", "hetero-rate", "Heterozygote rate.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mmq", "min-map-quality", "Minimum base call (mapping) quality for a match to be considered.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("ch", "corrected-het", "Use amplification bias corrected distribution for heterozygotes. Default: off."));
    addOption(parser, ArgParseOption("maf", "mean-alleleFreq", "Mean ref allele frequency in heterozygotes.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("ac", "amp-cycles", "Number of amplification cycles.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("ae", "amp-efficiency", "Polymerase efficiency, probability of amplification.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("in", "initial-N", "Initial allele population size.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("pht", "print-hetTable", "Print het table. Default: off.")); 
    hideOption(parser, "pht");
    addOption(parser, ArgParseOption("mec", "min-explained-column", "Minimum fraction of alignment column reads explained by genotype call.", ArgParseArgument::DOUBLE));

    addSection(parser, "Indel calling options: ");
    addOption(parser, ArgParseOption("it", "indel-threshold", "Minimal number of indel-supporting reads required for indel calling.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("ipt", "indel-perc-threshold", "Minimal ratio of indel-supporting/covering reads for indel to be called.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("iqt", "indel-quality-thresh", "Minimal average quality of inserted base/deletion-neighboring bases for indel to be called.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("bsi", "both-strands-indel", "Both strands need to be observed for indel to be called. Default: off."));
    addOption(parser, ArgParseOption("iw", "indel-window", "Overlap window used for indel calling.", ArgParseArgument::INTEGER)); 
    hideOption(parser, "iw");
    addOption(parser, ArgParseOption("ebi", "exclude-border-indel", "Same as option -eb but for indel candidates.", ArgParseArgument::INTEGER)); 
    addOption(parser, ArgParseOption("cws", "cnv-window-size", "CNV window size.", ArgParseArgument::INTEGER)); 
    hideOption(parser, "cws");
    setMaxValue(parser, "cnv-window-size", "10000");

    addSection(parser, "Bs options: ");
    addOption(parser, ArgParseOption("I", "intervals", "Intervals to analyse.",  ArgParseArgument::STRING, "STR", true));
    addOption(parser, ArgParseOption("spl", "beta-sampling", "Sample beta values (instead of newton)."));
    addOption(parser, ArgParseOption("ply", "polynomial", "Use polynomial prob. functiopn (instead of naive multiplication)."));
    addOption(parser, ArgParseOption("nsp", "normal-space", "Normal space (instead of LogSapce)."));
    addOption(parser, ArgParseOption("msc", "min-score", "Minimum score to call.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mpc", "min-prob", "Minimum genotype probability to call.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("umq", "use-mapq", "Use mapqs as weights for SNP/BS calling."));

    addOption(parser, ArgParseOption("gp", "genotype-priors", "Use non-uniform genotype prior probabilities (taking transversion etc. into account)."));
    addOption(parser, ArgParseOption("nec", "ns-errors-calling", "Use non-uniform sequencing error frequencies for calling."));

    // Realignment
    addOption(parser, ArgParseOption("nse", "ns-subst-errors", "Use non-uniform substitution error frequencies for realigning."));
    addOption(parser, ArgParseOption("nsi", "ns-ins-errors", "Use non-uniform insertion error frequencies for realigning."));
    addOption(parser, ArgParseOption("nsd", "ns-del-errors", "Use non-uniform deletion error frequencies for realigning."));

    addOption(parser, ArgParseOption("dr", "del-rate", "Genomic deletion rate.", ArgParseArgument::DOUBLE)); 
    addOption(parser, ArgParseOption("der", "del-error-rate", "Deletion error rate.", ArgParseArgument::DOUBLE)); 
    addOption(parser, ArgParseOption("ier", "ins-error-rate", "Insertion error rate.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("egs", "end-gap-score", "Simple score for end gaps (must be in the range of internal gaps) to avoid introducing of gaps.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("sl", "score-limit", "Score limit to avoid to high negative scores in profile realignment.", ArgParseArgument::DOUBLE));

    addSection(parser, "Other options: ");
    addOption(parser, ArgParseOption("lf", "log-file", "Write log file to FILE.", ArgParseArgument::OUTPUTFILE));
    setValidValues(parser, "lf", "log");
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    // Options:
    getOptionValue(options.outputSNP, parser, "output"); 
    getOptionValue(options.inputFormat, parser, "input-format");
    getOptionValue(options.outputFormat, parser, "output-format");
    if (isSet(parser, "dont-clip"))
        options.dontClip = true;
    if (isSet(parser, "multi"))
        options.keepMultiReads = true;
    if (isSet(parser, "hide-qualities"))
        options.showQualityStrings = false;
    if (isSet(parser, "solexa-qual-offset"))
        options.asciiQualOffset = 64;
    getOptionValue(options.outputIndel, parser, "indel-file");
    getOptionValue(options.method, parser, "method");
    getOptionValue(options.maxPile, parser, "max-pile");
    if (isSet(parser, "merged-max-pile"))
        options.laneSpecificMaxPile = false;
    getOptionValue(options.minCoverage, parser, "min-coverage");
    getOptionValue(options.forceCallCount, parser, "force-call");
    if (isSet(parser, "orientation-aware"))
        options.orientationAware = true;
    getOptionValue(options.outputCNV, parser, "output-cnv");
    getOptionValue(options.outputPosition, parser, "output-positions");
    getOptionValue(options.inputPositionFile, parser, "input-positions");
    getOptionValue(options.maxPolymerRun, parser, "max-polymer-run");
    getOptionValue(options.minDifferentReadPos, parser, "diff-pos");
    getOptionValue(options.excludeBorderPos, parser, "exclude-border");
    if (isSet(parser, "suboptimal"))
        options.keepSuboptimalReads = true;
    if (isSet(parser, "realign"))
        options.realign = true;
    getOptionValue(options.newQualityCalibrationFactor, parser, "corrected-quality");
    getOptionValue(options.windowSize, parser, "parse-window-size");
    getOptionValue(options.realignAddBorder, parser, "realign-border");
    // SNP Calling Options:
    getOptionValue(options.minMutT, parser, "min-mutations");
    getOptionValue(options.percentageT, parser, "perc-threshold");
    getOptionValue(options.avgQualT, parser, "min-quality");
    getOptionValue(options.theta, parser, "theta");
    getOptionValue(options.hetRate, parser, "hetero-rate");
    getOptionValue(options.minMapQual, parser, "min-map-quality");
    if (isSet(parser, "corrected-het"))
        options.correctedHetTable = true;
    getOptionValue(options.meanAlleleFrequency, parser, "mean-alleleFreq");
    getOptionValue(options.amplificationCycles, parser, "amp-cycles");
    getOptionValue(options.amplificationEfficiency, parser, "amp-efficiency");
    getOptionValue(options.initialN, parser, "initial-N");
    if (isSet(parser, "print-hetTable"))
        options.printHetTable = true;
    getOptionValue(options.minExplainedColumn, parser, "min-explained-column");
    // Indel Calling Options:
    getOptionValue(options.indelCountThreshold, parser, "indel-threshold");
    getOptionValue(options.indelPercentageT, parser, "indel-perc-threshold");
    getOptionValue(options.indelQualityThreshold, parser, "indel-quality-thresh");
    if (isSet(parser, "both-strands-indel"))
        options.bothIndelStrands = true;
    getOptionValue(options.indelWindow, parser, "indel-window");
    getOptionValue(options.indelDepthMinOverlap, parser, "exclude-border-indel");
    getOptionValue(options.cnvWindowSize, parser, "cnv-window-size");
    // Bs Options:
    clear(methOptions.intervals);
    if (isSet(parser, "intervals"))
    {
        resize(methOptions.intervals, getOptionValueCount(parser, "intervals"), Exact());
        for (unsigned i = 0; i < getOptionValueCount(parser, "intervals"); ++i)
        {
            CharString startPos;
            CharString endPos;
            CharString value;
            getOptionValue(value, parser, "intervals", i);
            unsigned j;
            for (j = 0; value[j] != ':' && j < length(value); ++j)
            {
                appendValue(methOptions.intervals[i].contigName, value[j], Generous());
            }
            ++j;
            for (; value[j] != '-' && j < length(value); ++j)
            {
                appendValue(startPos, value[j], Generous());
            }
            ++j;
            methOptions.intervals[i].startPos = lexicalCast<unsigned>(startPos);
            for (; j < length(value); ++j)
            {
                appendValue(endPos, value[j], Generous());
            }
            methOptions.intervals[i].endPos = lexicalCast<unsigned>(endPos);

            std::cout << " contigName: " << methOptions.intervals[i].contigName << "   startPos: " << methOptions.intervals[i].startPos << "  endPos: " << methOptions.intervals[i].endPos;
        }
    }
    if (isSet(parser, "beta-sampling"))
        methOptions.betaSampling = true;
    else
        methOptions.betaSampling = false;

    if (isSet(parser, "polynomial"))
        methOptions.polynomialProbFunction = true;
    else
    {
        methOptions.polynomialProbFunction = false;
    }
    if (isSet(parser, "normal-space")) methOptions.logSpace = false; 
    else methOptions.logSpace = true;
    
    getOptionValue(methOptions.minScoreToCallSnp, parser, "min-score");
    getOptionValue(methOptions.minProbToCallSnp, parser, "min-prob");

    if (isSet(parser, "use-mapq")) methOptions.useMapq = true;
    else methOptions.useMapq = false;

    if (isSet(parser, "genotype-priors")) methOptions.uniformGenPriors = false;
    if (isSet(parser, "ns-errors-calling")) methOptions.uniformSeqErrorsCalling = false;
    if (isSet(parser, "ns-subst-errors")) methOptions.nonSimpleSubstErrors = true;
    if (isSet(parser, "ns-ins-errors")) methOptions.nonSimpleInsErrors = true;
    if (isSet(parser, "ns-del-errors")) methOptions.nonSimpleDelErrors = true;

    getOptionValue(methOptions.delRate, parser, "del-rate");
    getOptionValue(methOptions.delErrorRate, parser, "del-error-rate");
    getOptionValue(methOptions.insErrorRate, parser, "ins-error-rate");
    getOptionValue(methOptions.endGapScore, parser, "end-gap-score");
    getOptionValue(methOptions.scoreLimit, parser, "score-limit");

    // Other Options:
    getOptionValue(options.outputLog, parser, "log-file");
    if (isSet(parser, "verbose"))
        options._debugLevel = max(options._debugLevel, 1);
    if (isSet(parser, "very-verbose"))
        options._debugLevel = max(options._debugLevel, 2);
    
    getArgumentValue(options.genomeFName, parser, 0);
    unsigned countFiles = getArgumentValueCount(parser, 1);
    resize(options.readFNames, countFiles);
    for (unsigned i = 0; i < countFiles; ++i)
        getArgumentValue(options.readFNames[i], parser, 1, i);

    if(options.windowSize > 1000000) 
        options.windowSize = 10000;
    
    
    if(options.runID == "")
    {
        ::std::string tempStr = toCString(options.readFNames[0]);
        size_t lastPos = tempStr.find_last_of('/') + 1;
        if (lastPos == tempStr.npos) lastPos = tempStr.find_last_of('\\') + 1;
        if (lastPos == tempStr.npos) lastPos = 0;
        options.runID = tempStr.substr(lastPos);
    }

    return ArgumentParser::PARSE_OK;
}



int main(int argc, const char *argv[]) 
{
    ArgumentParser parser;
    SNPCallingOptions<>     options;
    MethCallingOptions      methOptions;
    ArgumentParser::ParseResult res = parseCommandLine(options, methOptions, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;
    
    for(int arg = 0; arg < argc; ++arg) {
        options.programCall << argv[arg] << " ";
    }

    //////////////////////////////////////////////////////////////////////////////
    // check for variants
#ifdef CALL_PROFILE 
    double timeStamp = sysTime();
#endif  
    int result = detectSNPs(options, methOptions);
    if (result > 0)
    {
        cerr << "ERROR: Something went wrong. Try 'snpStore --help' for more information." << endl << endl;
        return 0;
    }
#ifdef CALL_PROFILE
    Times::instance().time_all = sysTime() - timeStamp;
    std::cout << "  Time needed for all: " << Times::instance().time_all/60.0 << "min" << std::endl;
    std::cout << "  Time needed for doBsCalling: " << Times::instance().time_doBsCalling/60.0 << "min" << std::endl;
    std::cout << "  Time needed for optimization: " << Times::instance().time_optimization/60.0 << "min" << std::endl;
    std::cout << "  Time needed for IO: " << Times::instance().time_IO/60.0 << "min" << std::endl;
    std::cout << "  Time needed for convertPairwiseToGlobal: " << Times::instance().time_convertPWToGlobal/60.0 << "min" << std::endl;
#endif

    std::cout << "CounteBViolated: " << methOptions.counteBViolated << std::endl; 
    std::cout << "CountPlanB (Newton, error bound violated or f'' > 0: sampling applied): " << methOptions.countPlanB << std::endl; 
    std::cout << "CountNoPlanB: " << methOptions.countNoPlanB << std::endl;
    std::cout << "CountCovTooLow: " << methOptions.countCovTooLow << std::endl;
    std::cout << "CountScoreTooLow: " << methOptions.countScoreTooLow << std::endl;


    return result;
}

