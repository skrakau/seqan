#ifndef SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_2_H
#define SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_2_H

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>


using namespace seqan;



template<typename TMethLevels, typename TContigNames, typename TFileName>
inline bool
loadMethLevels(TMethLevels &methLevelsTop, TMethLevels &methLevelsBottom, TContigNames &contigNames, TFileName &fileName)
{
    clear(contigNames);
    clear(methLevelsTop);
    clear(methLevelsBottom);
    SequenceStream streamRead(toCString(fileName)); 

    CharString id;
    CharString seq;
    bool prevTop = false;
    if (!isGood(streamRead))
    {
        std::cerr << "ERROR: Could not open methylation level FASTA file.\n";
        return 1;
    }
    while (!atEnd(streamRead))
    {
        if (readRecord(id, seq, streamRead) != 0)                                
        {
            std::cerr << "ERROR: Could not read from methylation level FASTA file!\n";
            return 1;
        }
        if (suffix(id, length(id)-3) == "TOP")
        {
            SEQAN_ASSERT(!prevTop);
            appendValue(contigNames, prefix(id, length(id)-4));
            appendValue(methLevelsTop, seq);
            prevTop = true;
        }
        else if (suffix(id, length(id)-3) == "BOT")
        {
            SEQAN_ASSERT(prevTop);
            SEQAN_ASSERT_EQ(back(contigNames), prefix(id, length(id)-4)); 
            appendValue(methLevelsBottom, seq);
            prevTop = false;
        }
        else
        {
            std::cerr << "ERROR: Something wrong with methylation level FASTA file!\n";
            return 1;
        }
    }
    return 0;
}

inline double
getMethLevel(char &c)
{
    if (c < '>')  // '>' cannot be used as value
        return ((double)(c - 33)* 0.0125);
    else
        return ((double)(c - 34)* 0.0125);
}

    
template<typename TOptions, typename TSnpFormatTag>
bool benchmark(TOptions &options, TSnpFormatTag const &, MethMason2 const &)
{
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;

    std::fstream simSnpFile(toCString(options.simSnpsFile), std::ios::binary | std::ios::in);
    TRecordReader simSnpReader(simSnpFile);

    std::fstream callFile(toCString(options.calledSnpsFile), std::ios::binary | std::ios::in);
    TRecordReader callReader(callFile);

    String<CharString> methLevelsTop;
    String<CharString> methLevelsBottom;
    String<CharString> simMethContigNames; 

    // read fastas
    loadMethLevels(methLevelsTop, methLevelsBottom, simMethContigNames, options.simMethsFile);
    if (empty(simMethContigNames) || empty(methLevelsTop) || empty(methLevelsBottom)) 
    {
        std::cerr << "ERROR: File containing simulated methylation levels is empty pr something went wrong.\n";
        return 1;
    }
    typedef Iterator<String<CharString> >::Type TMethContigIter;
    TMethContigIter itTop(begin(methLevelsTop));
    TMethContigIter itTopEnd(end(methLevelsTop));
    TMethContigIter itBottom(begin(methLevelsBottom));
    TMethContigIter itContigName(begin(simMethContigNames));
    bool simMethsAtEnd = false;

    // First skip lines with metainformation
    CharString helper;
    CharString simSnpContig;
    clear(simSnpContig);
    CharString callContig;
    clear(callContig);
    CharString simMethContig;
    clear(simMethContig);

    clear(helper);
    readNChars(helper, simSnpReader, 1);
    while(!atEnd(simSnpReader) && helper[0] == '#')
    {
        skipLine(simSnpReader);
        clear(helper);
        readUntilWhitespace(helper, simSnpReader);
    }
    simSnpContig = helper;
    clear(helper);
    readNChars(helper, callReader, 1);
    while(!atEnd(callReader) && helper[0] == '#')
    {
       skipLine(callReader);
       clear(helper);
       readUntilWhitespace(helper, callReader);
    }
    callContig = helper;
    clear(helper);
    simMethContig = simMethContigNames[0];
   
    unsigned simSnpPos = 0;
    unsigned callPos = 0;

    unsigned c_SimSnps = 0;     // TODO: struct stats
    unsigned c_CalledSnps = 0;
    unsigned c_NotListed = 0;    // Not found snps
    unsigned c_NotCalled = 0;    // 
    unsigned c_WrongCalled = 0;
    unsigned c_RightCalled = 0;
    unsigned c_FalsePositive = 0;    // False found snps
    unsigned c_FalsePositiveAtMeth = 0;    // False called snp at meth position

    String<unsigned> methDists;
    String<unsigned> methDistsAtSnps;
    resize(methDists, 20, 0);
    resize(methDistsAtSnps, 20, 0);

    unsigned maxPos = 0;

    clear(helper);
    skipWhitespaces(simSnpReader);
    readUntilWhitespace(helper, simSnpReader);
    simSnpPos = lexicalCast<size_t>(helper);

    clear(helper);
    skipWhitespaces(callReader);
    readUntilWhitespace(helper, callReader);
    callPos = lexicalCast<size_t>(helper);
    clear(helper);

    if (atEnd(simSnpReader)) simSnpPos = maxPos;
    if (atEnd(callReader)) callPos = maxPos;

    CharString currContig = callContig;
    if (atEnd(simSnpReader)) std::cout << " simSnpReader at end!" << std::endl;
    if (atEnd(callReader)) std::cout << " callReader at end!" << std::endl;
    if (atEnd(simSnpReader) || atEnd(callReader) )
    {
        std::cerr << "ERROR: One of the input files is empty. " << std::endl;
        return 1;
    }
    while(!simMethsAtEnd)
    {
        /*if (simSnpPos > 1000000 || (callPos >= 508900 && callPos < 508980))
        {
            std::cout << "ERROR: Should be at end, simSnpPos " << simSnpPos << " callPos: " << callPos <<"  maxPos: " << maxPos << std::endl;
            if (atEnd(simSnpReader)) std::cout << " simSnpReader at end!" << std::endl;
            if (atEnd(callReader)) std::cout << " callReader at end!" << std::endl;
            if (atEnd(simMethReader)) std::cout << " simMethReader at end!" << std::endl;
        }*/

        // Get highest position (documented)
        if (!atEnd(simSnpReader) && simSnpPos >= maxPos) maxPos = simSnpPos +1;
        if (!atEnd(callReader) && callPos >= maxPos) maxPos = callPos+1;

        // If reader reached end of file, increase corresponding 'position' so that its not taking into account anymore
        if (atEnd(simSnpReader)) simSnpPos = maxPos;
        if (atEnd(callReader)) callPos = maxPos;

        if (simSnpContig == callContig && simSnpContig == simMethContig)
        {
            currContig = callContig;
        }
        else
        {
            if (simSnpContig != currContig) simSnpPos = maxPos+1;   //??
            if (callContig != currContig) callPos = maxPos+1;
        }

        // TODO check if methLevels simulated but not called
        if (atEnd(simSnpReader)) std::cout << " simSnpReader at end!" << std::endl;
        if (atEnd(callReader)) std::cout << " callReader at end!" << std::endl;
        //std::cout << " simSnpPos " << simSnpPos << " callPos: " << callPos << "  maxPos: " << maxPos << std::endl;

        if (simSnpPos < callPos)     // Sim snp not called, iterate to next sim pos
        {   
            Dna5 refAlleleSimSnp;
            Dna5String genotypeSimSnp;
            if (readSimSnpRecord2(refAlleleSimSnp, genotypeSimSnp, simSnpReader, TSnpFormatTag()) == 1)
            {
                std::cout << "ERROR: b) Something went wrong with sim snp record at pos: " << simSnpPos << " callPos: " << callPos << std::endl;
                return 1;
            }
            if (refAlleleSimSnp !=  genotypeSimSnp[0] || refAlleleSimSnp !=  genotypeSimSnp[1])     // no snp simulated, something else
            {
                ++c_NotListed;
                ++c_SimSnps;
            }
            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader, TSnpFormatTag());
        }
        else if (simSnpPos > callPos &&
                callPos < maxPos )    // compare simulated meth level and called meth level, no snp simulated at this position
        {
            double methLevelTopSimMeth = -1;
            double methLevelBottomSimMeth = -1;
            methLevelTopSimMeth = getMethLevel((*itTop)[callPos]);
            methLevelBottomSimMeth = getMethLevel((*itBottom)[callPos]); 

            //std::cout << "callPos: " << callPos << "  methLevelTopSimMeth: " <<  methLevelTopSimMeth << "  methLevelBottomSimMeth: " << methLevelBottomSimMeth << std::endl;
            // callReader
            Dna5 refAlleleCall;
            Dna5String genotypeCall;
            double methLevelTopCall;
            double methLevelBottomCall;
            readCallRecord2(refAlleleCall, genotypeCall, methLevelTopCall, methLevelBottomCall, callReader);           
           
            if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
            {
                ++c_FalsePositive;
                ++c_FalsePositiveAtMeth;
                ++c_CalledSnps;
                //std::cout << "False positive at meth position, callPos: " << callPos << " ref: " <<  refAlleleCall << " called genotype: " << genotypeCall <<  std::endl; 
            }
            else 
            {
                if (methLevelTopCall < -0.5 && methLevelTopSimMeth > 0.0)  std::cerr << "Simulated top meth level, but not called: " << callPos << std::endl;
                if (methLevelBottomCall < -0.5 && methLevelBottomSimMeth > 0.0)  std::cerr << "Simulated bottom meth level but not called: " << callPos << std::endl;
                if (methLevelTopCall >= 0) incrementMethDist(methDists, methLevelTopSimMeth, methLevelTopCall);
                if (methLevelBottomCall >= 0) incrementMethDist(methDists, methLevelBottomSimMeth, methLevelBottomCall);
            }

            if (callPos + 1 >= length(*itTop))  // end of contig simMethLevels, iter to next contig meths
            {
                ++itTop;
                ++itBottom;
                ++itContigName;
                if (itTop == itTopEnd) simMethsAtEnd = true;
                else simMethContig = (*itContigName);
            }
            readCallRecord1(callContig, callPos, callReader);
        }
        else if (simSnpPos == callPos &&
                 simSnpPos < maxPos )    // compare simulated snp, simulated meth level and called snp and meth level
        {
            Dna5 refAlleleSimSnp;
            Dna5String genotypeSimSnp;
            if (readSimSnpRecord2(refAlleleSimSnp, genotypeSimSnp, simSnpReader, TSnpFormatTag()) == 1)
            {
                std::cout << "ERROR: b) Something went wrong with sim snp record at pos: " << simSnpPos << " callPos: " << callPos << std::endl;
                return 1;
            }
            double methLevelTopSimMeth = -1;
            double methLevelBottomSimMeth = -1;
            methLevelTopSimMeth = getMethLevel((*itTop)[callPos]);
            methLevelBottomSimMeth = getMethLevel((*itBottom)[callPos]); 

            // callReader
            Dna5 refAlleleCall;
            Dna5String genotypeCall;
            double methLevelTopCall;
            double methLevelBottomCall;
            readCallRecord2(refAlleleCall, genotypeCall, methLevelTopCall, methLevelBottomCall, callReader);      

            if (refAlleleSimSnp !=  genotypeSimSnp[0] || refAlleleSimSnp !=  genotypeSimSnp[1])     // no snp simulated, something else
            {
                ++c_SimSnps;
                
                CharString buffer;
                resize(buffer, 1024);
                sprintf(&buffer[0], "At position: %u", callPos);
                SEQAN_ASSERT_EQ_MSG(refAlleleCall, refAlleleSimSnp, toCString(buffer));

                if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
                {
                    ++c_CalledSnps;
                    if (genotypeCall[0] != genotypeSimSnp[0] || genotypeCall[1] != genotypeSimSnp[1])
                    {
                        ++c_WrongCalled;
                        //std::cout << "Wrong call at snp pos " << " callPos: " << callPos << " Sim: " << genotypeSimSnp[0]  << genotypeSimSnp[1] << " Call: " <<  genotypeCall[0] << genotypeCall[1] << std::endl; 
                    }
                    else
                    {
                        ++c_RightCalled;
                       if (methLevelTopCall < -0.5 && methLevelTopSimMeth > 0.0)  std::cerr << "Simulated top meth level, but not called: " << callPos << std::endl;
                       if (methLevelBottomCall < -0.5 && methLevelBottomSimMeth > 0.0)  std::cerr << "Simulated bottom meth level but not called: " << callPos << std::endl;
                       if (methLevelTopCall >= 0) incrementMethDist(methDistsAtSnps, methLevelTopSimMeth, methLevelTopCall);
                       if (methLevelBottomCall >= 0) incrementMethDist(methDistsAtSnps, methLevelBottomSimMeth, methLevelBottomCall);
                    }
                }
                else ++c_NotCalled;
            }
            else
            {
                if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
                {
                    ++c_FalsePositive;
                    ++c_FalsePositiveAtMeth;
                    ++c_CalledSnps;
                    //std::cout << "False positive at meth position, callPos: " << callPos << " ref: " <<  refAlleleCall << " called genotype: " << genotypeCall <<  std::endl; 
                }
                else 
                {
                    if (methLevelTopCall < -0.5 && methLevelTopSimMeth > 0.0)  std::cerr << "Simulated top meth level, but not called: " << callPos << std::endl;
                    if (methLevelBottomCall < -0.5 && methLevelBottomSimMeth > 0.0)  std::cerr << "Simulated bottom meth level but not called: " << callPos << std::endl;
                    if (methLevelTopCall >= 0) incrementMethDist(methDists, methLevelTopSimMeth, methLevelTopCall);
                    if (methLevelBottomCall >= 0) incrementMethDist(methDists, methLevelBottomSimMeth, methLevelBottomCall);
                }
           }

            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader, TSnpFormatTag());
            readCallRecord1(callContig, callPos, callReader);
            if (callPos + 1 >= length(*itTop))  // end of contig simMethLevels, iter to next contig meths
            {
                ++itTop;
                ++itBottom;
                ++itContigName;
                if (itTop == itTopEnd) simMethsAtEnd = true;
                else simMethContig = (*itContigName);
            }
            continue;
            
        }
        else if (simSnpPos >= maxPos && simSnpPos == callPos)
        {
            if (atEnd(simSnpReader) && atEnd(callReader)) simMethsAtEnd = true;

            if (callPos + 1 >= length(*itTop))  // end of contig simMethLevels, iter to next contig meths
            {
                ++itTop;
                ++itBottom;
                ++itContigName;
                if (itTop == itTopEnd) simMethsAtEnd = true;
                else simMethContig = (*itContigName);
            }
            continue;
        }

    }


    std::ofstream out(toCString(options.outputFile), std::ios::binary | std::ios::out);
    if (!out.good())
        std::cerr << " ERROR: Could not open output file!\n";

    out << "Simulated snps: \t" << c_SimSnps << '\n';
    out<< "Called snps: \t" << c_CalledSnps << '\n';
    out<< '\n';
    out<< "Not listed: \t" << c_NotListed << '\n';
    out<< "Not called: \t" << c_NotCalled << '\n';
    out<< "False positive: \t" << c_FalsePositive << '\n';
    out<< "False positive at meth positions: \t" << c_FalsePositiveAtMeth << '\n';
    out<< "Wrong called genotype: \t" << c_WrongCalled << '\n'; 
    out<< "Right called genotype: \t" << c_RightCalled << '\n';
    
    // Sensitivity
    out<< '\n';
    //out<< "Sensitivity: " << (double)c_RightCalled/(double)(c_RightCalled + c_FalsePositive) << '\n';
    out<< "Sensitivity: " << (double)c_RightCalled/(double)(c_RightCalled + (c_NotCalled + c_NotListed)) << '\n';   // Can be caused by to low coverage


    // Specificity
    // No. of true negatives: TODO wrong called not taken into account?
    unsigned trueNegatives = maxPos - c_FalsePositive;  
    out<< "Specificity: "  << (double)trueNegatives/(double)(trueNegatives + c_FalsePositive) <<  '\n';
    out<< "Precision: "  << (double)c_RightCalled /(double)(c_RightCalled + c_FalsePositive + c_WrongCalled) <<  '\n';



    // Methylation calls 
    out << '\n';
    out << "Methylation levels at Snp free positions:" << '\n';
    out << "Distance:\t" << "Count:\n";
    out << "<0.000001:\t" << methDists[0] << '\n';
    out << "<0.000005:\t" << methDists[1] << '\n';
    out << "<0.00001:\t"  << methDists[2] << '\n';
    out << "<0.00005:\t"  << methDists[3] << '\n';
    out << "<0.0001:\t"   << methDists[4] << '\n';
    out << "<0.0005:\t"   << methDists[5] << '\n';
    out << "<0.001:\t"    << methDists[6] << '\n';
    out << "<0.005:\t"    << methDists[7] << '\n';
    out << "<0.01:\t"     << methDists[8] << '\n';
    out << "<0.05:\t"     << methDists[9] << '\n';

    out << "<0.1:\t" << methDists[10] << '\n';
    out << "<0.2:\t" << methDists[11] << '\n';
    out << "<0.3:\t" << methDists[12] << '\n';
    out << "<0.4:\t" << methDists[13] << '\n';
    out << "<0.5:\t" << methDists[14] << '\n';
    out << "<0.6:\t" << methDists[15] << '\n';
    out << "<0.7:\t" << methDists[16] << '\n';
    out << "<0.8:\t" << methDists[17] << '\n';
    out << "<0.9:\t" << methDists[18] << '\n';
    out << "<1.0:\t" << methDists[19] << '\n';

    out << "Methylation levels at positions where Snp occured:" << '\n';
    out << "Distance:\t" << "Count:\n";

    out << "<0.000001:\t" << methDistsAtSnps[0] << '\n';
    out << "<0.000005:\t" << methDistsAtSnps[1] << '\n';
    out << "<0.00001:\t"  << methDistsAtSnps[2] << '\n';
    out << "<0.00005:\t"  << methDistsAtSnps[3] << '\n';
    out << "<0.0001:\t"   << methDistsAtSnps[4] << '\n';
    out << "<0.0005:\t"   << methDistsAtSnps[5] << '\n';
    out << "<0.001:\t"    << methDistsAtSnps[6] << '\n';
    out << "<0.005:\t"    << methDistsAtSnps[7] << '\n';
    out << "<0.01:\t"     << methDistsAtSnps[8] << '\n';
    out << "<0.05:\t"     << methDistsAtSnps[9] << '\n';

    out << "<0.1:\t" << methDistsAtSnps[10] << '\n';
    out << "<0.2:\t" << methDistsAtSnps[11] << '\n';
    out << "<0.3:\t" << methDistsAtSnps[12] << '\n';
    out << "<0.4:\t" << methDistsAtSnps[13] << '\n';
    out << "<0.5:\t" << methDistsAtSnps[14] << '\n';
    out << "<0.6:\t" << methDistsAtSnps[15] << '\n';
    out << "<0.7:\t" << methDistsAtSnps[16] << '\n';
    out << "<0.8:\t" << methDistsAtSnps[17] << '\n';
    out << "<0.9:\t" << methDistsAtSnps[18] << '\n';
    out << "<1.0:\t" << methDistsAtSnps[19] << '\n';

    return 0;
}








    

 
#endif
