#ifndef SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_PARSE_METH_LEVEL_H
#define SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_PARSE_METH_LEVEL_H

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

inline void
incrementMethLevels(String<unsigned> &methLevels, double &level)
{
    if      (level < 0.05) ++methLevels[0];
    else if (level < 0.1) ++methLevels[1];
    else if (level < 0.15) ++methLevels[2];
    else if (level < 0.2) ++methLevels[3];
    else if (level < 0.25) ++methLevels[4];
    else if (level < 0.3) ++methLevels[5];
    else if (level < 0.35) ++methLevels[6];
    else if (level < 0.4) ++methLevels[7];
    else if (level < 0.45) ++methLevels[8];
    else if (level < 0.5) ++methLevels[9];
    else if (level < 0.55) ++methLevels[10];
    else if (level < 0.6) ++methLevels[11];
    else if (level < 0.65) ++methLevels[12];
    else if (level < 0.7) ++methLevels[13];
    else if (level < 0.75) ++methLevels[14];
    else if (level < 0.8) ++methLevels[15];
    else if (level < 0.85) ++methLevels[16];
    else if (level < 0.9) ++methLevels[17];
    else if (level < 0.95) ++methLevels[18];
    else ++ methLevels[19];
}

template<typename TPos, typename TReader>
inline void
readBedRecord1(CharString &contigName, TPos &pos, TReader &reader)
{
    clear(contigName);
    readUntilWhitespace(contigName, reader);  // contigName
    skipWhitespaces(reader);

    CharString posStr;
    readUntilWhitespace(posStr, reader);     // pos
    pos = lexicalCast<TPos>(posStr);
    //--pos;      // TODO check 0, 1 based for all
}


template<typename TReader>
inline bool
readBedRecord2(double &methLevelTop, double &methLevelBottom, unsigned &cov, TReader &reader)
{
    methLevelTop = -1.0;
    methLevelBottom = -1.0;
    CharString helper;
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // End pos
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);        // Meth. level
    CharString methLevel = helper;

    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);     // Coverage
    cov = lexicalCast<unsigned>(helper);
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);     // Orientation
    if (helper == "+")
        methLevelTop = lexicalCast<double>(methLevel)/100.0;
    else if (helper == "-")
        methLevelBottom = lexicalCast<double>(methLevel)/100.0;
    else
        std::cerr << " ERROR: Something wrong with BED file (strand expected, but " << helper << " given!\n";

    skipLine(reader);

    return 0;
}




template<typename TOptions>
bool parse_methLevel_My(TOptions &options)
{
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;

    // bisSnp
    CharString fileNameIn = options.calledSnpsFile;
    std::fstream fileIn(toCString(fileNameIn), std::ios::binary | std::ios::in);
    TRecordReader reader(fileIn);

    // output
    CharString fileNamePlot = options.calledSnpsFile;
    append(fileNamePlot, ".levels");
    std::ofstream outPlot(toCString(fileNamePlot), std::ios::binary | std::ios::out);
    if (!outPlot.good())
        std::cerr << " ERROR: Could not open output plot file!\n";

    // output stats
    CharString fileNameStats = options.outputFile;
    std::fstream outFile(toCString(fileNameStats), std::ios::binary | std::ios::out | std::ios::app);
    if (!outFile.good())
        std::cerr << " ERROR: Could not open output file!\n";

    String<CharString> methLevelsTop;
    String<CharString> methLevelsBottom;
    String<CharString> simMethContigNames; 

    // read fastas
    loadMethLevels(methLevelsTop, methLevelsBottom, simMethContigNames, options.simMethsFile);
    if (empty(simMethContigNames) || empty(methLevelsTop) || empty(methLevelsBottom)) 
    {
        std::cerr << "ERROR: File containing simulated methylation levels is empty or something went wrong.\n";
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
    CharString simMethContig;
    clear(simMethContig);
    CharString contig;
    clear(contig);

    clear(helper);
    simMethContig = simMethContigNames[0];
    readNChars(helper, reader, 1);
    while(!atEnd(reader) && helper[0] == '#')
    { 
        skipLine(reader);
        clear(helper);
        readUntilWhitespace(helper, reader);
    }
    contig = helper;

    unsigned pos = 0;
    String<unsigned> methDists;
    resize(methDists, 20, 0);
    String<unsigned> methLevelsSim;
    resize(methLevelsSim, 20, 0);

    //unsigned maxPos = 0;

    clear(helper);
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);
    pos = lexicalCast<size_t>(helper);
    clear(helper);
    CharString currContig = contig;

    unsigned countNotCalled = 0;
    unsigned countMethCalls = 0;
    if (atEnd(reader)) std::cout << " reader at end!" << std::endl;
    while(!simMethsAtEnd && !atEnd(reader))
    {
        if (pos < length(*itTop) && contig == currContig) 
        {
            double methLevelTopSimMeth = -1;
            double methLevelBottomSimMeth = -1;
            methLevelTopSimMeth = getMethLevel((*itTop)[pos]);
            methLevelBottomSimMeth = getMethLevel((*itBottom)[pos]); 
            Dna5 refAllele;
            Dna5String genotypeCall;
            double methLevelTopCall;
            double methLevelBottomCall;
            unsigned cov;
            readCallRecord2(refAllele, genotypeCall, methLevelTopCall, methLevelBottomCall, cov, reader); 

            if (cov >= options.minCovMeth)
            {
                if (methLevelTopCall < -0.5 && methLevelTopSimMeth > 0.0)
                {
                    ++countNotCalled;    
                    //std::cerr << "Simulated top meth level, but not called: " << pos << std::endl;
                    continue;
                }
                if (methLevelBottomCall < -0.5 && methLevelBottomSimMeth > 0.0) 
                {
                    ++countNotCalled;  
                    //std::cerr << "Simulated bottom meth level but not called: " << pos << std::endl;
                    continue;
                }
                
                ++countMethCalls; 

                if (methLevelTopCall >= 0) incrementMethDist(methDists, methLevelTopSimMeth, methLevelTopCall);
                if (methLevelBottomCall >= 0) incrementMethDist(methDists, methLevelBottomSimMeth, methLevelBottomCall);

                // Output for dotplot (only if meth level called
                if (methLevelTopCall > -0.1) outPlot << methLevelTopSimMeth << '\t' << methLevelTopCall << '\n';
                if (methLevelBottomCall > -0.1) outPlot << methLevelBottomSimMeth << '\t' << methLevelBottomCall << '\n';

                if (refAllele == 'C') incrementMethLevels(methLevelsSim, methLevelTopSimMeth);
                if (refAllele == 'G') incrementMethLevels(methLevelsSim, methLevelBottomSimMeth);
            }
            readCallRecord1(contig, pos, reader);
        }
        else if (contig != currContig)
        {
            ++itTop;
            ++itBottom;
            ++itContigName;
            if (itTop == itTopEnd) simMethsAtEnd = true;
            else simMethContig = (*itContigName);
        }
        else
            std::cerr << " ERROR: This should not happen!\n";
    }
   
    // Methylation calls 
    std::cout << '\n';
    std::cout << "Methylation levels at Snp free positions:" << '\n';
    std::cout << "Distance:\t" << "Count:\n";
    std::cout << "<0.000001:\t" << methDists[0] << '\n';
    std::cout << "<0.000005:\t" << methDists[1] << '\n';
    std::cout << "<0.00001:\t"  << methDists[2] << '\n';
    std::cout << "<0.00005:\t"  << methDists[3] << '\n';
    std::cout << "<0.0001:\t"   << methDists[4] << '\n';
    std::cout << "<0.0005:\t"   << methDists[5] << '\n';
    std::cout << "<0.001:\t"    << methDists[6] << '\n';
    std::cout << "<0.005:\t"    << methDists[7] << '\n';
    std::cout << "<0.01:\t"     << methDists[8] << '\n';
    std::cout << "<0.05:\t"     << methDists[9] << '\n';

    std::cout << "<0.1:\t" << methDists[10] << '\n';
    std::cout << "<0.2:\t" << methDists[11] << '\n';
    std::cout << "<0.3:\t" << methDists[12] << '\n';
    std::cout << "<0.4:\t" << methDists[13] << '\n';
    std::cout << "<0.5:\t" << methDists[14] << '\n';
    std::cout << "<0.6:\t" << methDists[15] << '\n';
    std::cout << "<0.7:\t" << methDists[16] << '\n';
    std::cout << "<0.8:\t" << methDists[17] << '\n';
    std::cout << "<0.9:\t" << methDists[18] << '\n';
    std::cout << "<1.0:\t" << methDists[19] << '\n';

    // Simulated methylation levels
    std::cout << "Methylation levels simulated:" << '\n';
    std::cout << "Level:\t" << "Count:\n";
    std::cout << "<0.05:\t" << methLevelsSim[0] << '\n';
    std::cout << "<0.1:\t"  << methLevelsSim[1] << '\n';
    std::cout << "<0.15:\t" << methLevelsSim[2] << '\n';
    std::cout << "<0.2:\t"  << methLevelsSim[3] << '\n';
    std::cout << "<0.25:\t" << methLevelsSim[4] << '\n';
    std::cout << "<0.3:\t"  << methLevelsSim[5] << '\n';
    std::cout << "<0.35:\t" << methLevelsSim[6] << '\n';
    std::cout << "<0.4:\t"  << methLevelsSim[7] << '\n';
    std::cout << "<0.45:\t" << methLevelsSim[8] << '\n';
    std::cout << "<0.5:\t"  << methLevelsSim[9] << '\n';
    std::cout << "<0.55:\t" << methLevelsSim[10] << '\n';
    std::cout << "<0.6:\t"  << methLevelsSim[11] << '\n';
    std::cout << "<0.65:\t" << methLevelsSim[12] << '\n';
    std::cout << "<0.7:\t"  << methLevelsSim[13] << '\n';
    std::cout << "<0.75:\t" << methLevelsSim[14] << '\n';
    std::cout << "<0.8:\t"  << methLevelsSim[15] << '\n';
    std::cout << "<0.85:\t" << methLevelsSim[16] << '\n';
    std::cout << "<0.9:\t"  << methLevelsSim[17] << '\n';
    std::cout << "<0.95:\t" << methLevelsSim[18] << '\n';
    std::cout << "<1.00:\t" << methLevelsSim[19] << '\n';

    outFile << "snp_meth_store: No. of methylation calls: " << countMethCalls << '\n';
    outFile << "snp_meth_store: No. of positions where meth simulated and listed but no meth called: " << countNotCalled << '\n';

    return 0;
}

template<typename TOptions>
bool parse_methLevel_BisSNP(TOptions &options)
{
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;

    // bisSnp
    CharString fileNameIn = options.bedFile;
    std::fstream fileIn(toCString(fileNameIn), std::ios::binary | std::ios::in);
    TRecordReader reader(fileIn);

    // output
    CharString fileNamePlot = options.bedFile;
    append(fileNamePlot, ".levels");
    std::ofstream outPlot(toCString(fileNamePlot), std::ios::binary | std::ios::out);
    if (!outPlot.good())
        std::cerr << " ERROR: Could not open output plot file!\n";

    // output stats
    CharString fileNameStats = options.outputFile;
    std::fstream outFile(toCString(fileNameStats), std::ios::binary | std::ios::out | std::ios::app);
    if (!outFile.good())
        std::cerr << " ERROR: Could not open output file!\n";

    String<CharString> methLevelsTop;
    String<CharString> methLevelsBottom;
    String<CharString> simMethContigNames; 

    // read fastas
    loadMethLevels(methLevelsTop, methLevelsBottom, simMethContigNames, options.simMethsFile);
    if (empty(simMethContigNames) || empty(methLevelsTop) || empty(methLevelsBottom)) 
    {
        std::cerr << "ERROR: File containing simulated methylation levels is empty or something went wrong.\n";
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
    CharString simMethContig;
    clear(simMethContig);
    CharString contig;
    clear(contig);

    clear(helper);
    simMethContig = simMethContigNames[0];
    readUntilWhitespace(helper, reader); 
    while(!atEnd(reader) && (helper[0] == '#' || helper == "track"))
    { 
        skipLine(reader);
        clear(helper);
        readUntilWhitespace(helper, reader);
    }
    contig = helper;

    unsigned pos = 0;
    String<unsigned> methDists;
    resize(methDists, 20, 0);
    String<unsigned> methLevelsSim;
    resize(methLevelsSim, 20, 0);

    //unsigned maxPos = 0;

    clear(helper);
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);
    pos = lexicalCast<size_t>(helper);
    clear(helper);
    CharString currContig = contig;

    unsigned countNotCalled = 0;
    unsigned countMethCalls = 0;
    unsigned countListed = 0;
    unsigned countCovered = 0;
    unsigned countCalled = 0;

    if (atEnd(reader)) std::cout << " bed reader at end!" << std::endl;
    while(!simMethsAtEnd && !atEnd(reader))
    {
        if (pos < length(*itTop) && contig == currContig) 
        {
            //std::cout << "contigName: " << contig << " pos: " << pos << std::endl;
            ++countMethCalls; 
            double methLevelTopSimMeth = -1;
            double methLevelBottomSimMeth = -1;
            methLevelTopSimMeth = getMethLevel((*itTop)[pos]);
            methLevelBottomSimMeth = getMethLevel((*itBottom)[pos]); 
            Dna5 refAllele;
            Dna5String genotypeCall;
            double methLevelTopCall;
            double methLevelBottomCall;
            unsigned cov;
            readBedRecord2(methLevelTopCall, methLevelBottomCall, cov, reader);
            readBedRecord1(contig, pos, reader);
            ++countListed;
            if (cov >= options.minCovMeth/2)
            {
                //std::cout << "Call: " << methLevelTopCall << '\t' << methLevelBottomCall << '\n';
                //std::cout << " Sim: " << methLevelTopSimMeth << '\t' << methLevelBottomSimMeth << "\n\n";
                ++countCovered;
                if (methLevelTopCall < -0.5 && methLevelTopSimMeth > 0.0)
                {
                    ++countNotCalled;    
                    //std::cerr << "Simulated top meth level, but not called: " << pos << std::endl;
                    continue;
                }
                if (methLevelBottomCall < -0.5 && methLevelBottomSimMeth > 0.0) 
                {
                    ++countNotCalled;  
                    //std::cerr << "Simulated bottom meth level but not called: " << pos << std::endl;
                    continue;
                }
                ++countCalled;
                if (methLevelTopCall >= 0) incrementMethDist(methDists, methLevelTopSimMeth, methLevelTopCall);
                if (methLevelBottomCall >= 0) incrementMethDist(methDists, methLevelBottomSimMeth, methLevelBottomCall);
                
                // Output for dotplot (only if meth level called
                if (methLevelTopCall > -0.1) outPlot << methLevelTopSimMeth << '\t' << methLevelTopCall << '\n';
                if (methLevelBottomCall > -0.1) outPlot << methLevelBottomSimMeth << '\t' << methLevelBottomCall <<'\n';

                incrementMethLevels(methLevelsSim, methLevelTopSimMeth);
                incrementMethLevels(methLevelsSim, methLevelBottomSimMeth);
            }
        }
        else if (contig != currContig)
        {
            ++itTop;
            ++itBottom;
            ++itContigName;
            if (itTop == itTopEnd) simMethsAtEnd = true;
            else simMethContig = (*itContigName);
        }
        else
            std::cerr << " ERROR: This should not happen!\n";
    }
   
    // Methylation calls 
    std::cout << '\n';
    std::cout << "Methylation levels at Snp free positions:" << '\n';
    std::cout << "Distance:\t" << "Count:\n";
    std::cout << "<0.000001:\t" << methDists[0] << '\n';
    std::cout << "<0.000005:\t" << methDists[1] << '\n';
    std::cout << "<0.00001:\t"  << methDists[2] << '\n';
    std::cout << "<0.00005:\t"  << methDists[3] << '\n';
    std::cout << "<0.0001:\t"   << methDists[4] << '\n';
    std::cout << "<0.0005:\t"   << methDists[5] << '\n';
    std::cout << "<0.001:\t"    << methDists[6] << '\n';
    std::cout << "<0.005:\t"    << methDists[7] << '\n';
    std::cout << "<0.01:\t"     << methDists[8] << '\n';
    std::cout << "<0.05:\t"     << methDists[9] << '\n';

    std::cout << "<0.1:\t" << methDists[10] << '\n';
    std::cout << "<0.2:\t" << methDists[11] << '\n';
    std::cout << "<0.3:\t" << methDists[12] << '\n';
    std::cout << "<0.4:\t" << methDists[13] << '\n';
    std::cout << "<0.5:\t" << methDists[14] << '\n';
    std::cout << "<0.6:\t" << methDists[15] << '\n';
    std::cout << "<0.7:\t" << methDists[16] << '\n';
    std::cout << "<0.8:\t" << methDists[17] << '\n';
    std::cout << "<0.9:\t" << methDists[18] << '\n';
    std::cout << "<1.0:\t" << methDists[19] << '\n';

    // Simulated methylation levels
    std::cout << "Methylation levels simulated:" << '\n';
    std::cout << "Level:\t" << "Count:\n";
    std::cout << "<0.05:\t" << methLevelsSim[0] << '\n';
    std::cout << "<0.1:\t"  << methLevelsSim[1] << '\n';
    std::cout << "<0.15:\t" << methLevelsSim[2] << '\n';
    std::cout << "<0.2:\t"  << methLevelsSim[3] << '\n';
    std::cout << "<0.25:\t" << methLevelsSim[4] << '\n';
    std::cout << "<0.3:\t"  << methLevelsSim[5] << '\n';
    std::cout << "<0.35:\t" << methLevelsSim[6] << '\n';
    std::cout << "<0.4:\t"  << methLevelsSim[7] << '\n';
    std::cout << "<0.45:\t" << methLevelsSim[8] << '\n';
    std::cout << "<0.5:\t"  << methLevelsSim[9] << '\n';
    std::cout << "<0.55:\t" << methLevelsSim[10] << '\n';
    std::cout << "<0.6:\t"  << methLevelsSim[11] << '\n';
    std::cout << "<0.65:\t" << methLevelsSim[12] << '\n';
    std::cout << "<0.7:\t"  << methLevelsSim[13] << '\n';
    std::cout << "<0.75:\t" << methLevelsSim[14] << '\n';
    std::cout << "<0.8:\t"  << methLevelsSim[15] << '\n';
    std::cout << "<0.85:\t" << methLevelsSim[16] << '\n';
    std::cout << "<0.9:\t"  << methLevelsSim[17] << '\n';
    std::cout << "<0.95:\t" << methLevelsSim[18] << '\n';
    std::cout << "<1.00:\t" << methLevelsSim[19] << '\n';

    outFile << "Bis-SNP: No. of methylation calls: " << countMethCalls << '\n';
    outFile << "Bis-SNP: Not called: " << countNotCalled<< '\n';
    outFile << "Bis-SNP: listed: " << countListed<< '\n';
    outFile << "Bis-SNP: covered: " << countCovered<< '\n';
    outFile << "Bis-SNP: called: " << countCalled<< '\n';



    return 0;
}



#endif
