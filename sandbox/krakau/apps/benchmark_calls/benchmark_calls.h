#ifndef SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_H
#define SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_H

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>


using namespace seqan;

struct SnpMason2_;
typedef Tag<SnpMason2_> const SnpMason2;

struct SnpCustom_;
typedef Tag<SnpCustom_> const SnpCustom;

struct MethCustom_;
typedef Tag<MethCustom_> const MethCustom;

struct MethMason2_;
typedef Tag<MethMason2_> const MethMason2;


inline void
incrementMethDist(String<unsigned> &dists, double &methLevel1, double&methLevel2)
{
    double dist = std::abs(methLevel1 - methLevel2);
    if      (dist < 0.000001) ++dists[0];
    else if (dist < 0.000005) ++dists[1];
    else if (dist < 0.00001) ++dists[2];
    else if (dist < 0.00005) ++dists[3];
    else if (dist < 0.0001) ++dists[4];
    else if (dist < 0.0005) ++dists[5];
    else if (dist < 0.001) ++dists[6];
    else if (dist < 0.005) ++dists[7];
    else if (dist < 0.01) ++dists[8];
    else if (dist < 0.05) ++dists[9];
    else if (dist < 0.1) ++dists[10];
    else if (dist < 0.2) ++dists[11];
    else if (dist < 0.3) ++dists[12];
    else if (dist < 0.4) ++dists[13];
    else if (dist < 0.5) ++dists[14];
    else if (dist < 0.6) ++dists[15];
    else if (dist < 0.7) ++dists[16];
    else if (dist < 0.8) ++dists[17];
    else if (dist < 0.9) ++dists[18];
    else ++ dists[19];
}

template<typename TPos, typename TReader>
inline bool 
readCallRecord1(CharString &contigName, TPos &pos, TReader &reader)
{
    clear(contigName);
    readUntilWhitespace(contigName, reader);    // contigName
    skipWhitespaces(reader);

    CharString posStr;
    readUntilWhitespace(posStr, reader);       // pos
    pos = lexicalCast<TPos>(posStr);
    return 0;
}

template <typename TReader>
inline bool
readCallRecord2(Dna5 &refAllele, Dna5String &genotype, double &methLevelTop, double &methLevelBottom, TReader &reader)
{
    CharString helper;
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);        // ref
    refAllele = helper[0];
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qA
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qC
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qG
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qT
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qA
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qC
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qG
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qT
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // cov
    skipWhitespaces(reader);
    clear(genotype);
    readUntilWhitespace(genotype, reader);        // genotype called
    if (genotype == ".")
    {
        genotype[0] = refAllele;
        appendValue(genotype, refAllele);
    }
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);  // meth level called

    if (helper[0] != '.')
    {
        if (genotype[0] == 'C' && genotype[1] == 'G')
        {
            unsigned i = 0;
            while (i < length(helper) && helper[i] != ':')
                ++i;
            CharString pre = prefix(helper, i);
            methLevelTop = lexicalCast<double>(pre);
            CharString suff = suffix(helper, i+1);
            methLevelBottom = lexicalCast<double>(suff); 
        }
        else if (genotype[0] == 'C' || genotype[1] == 'C')
        {
            methLevelTop = lexicalCast<double>(helper);
            methLevelBottom = -1.0;
        }
        else if (genotype[0] == 'G' || genotype[1] == 'G')
        {
            methLevelTop = -1.0;
            methLevelBottom = lexicalCast<double>(helper);        
        }
    }
    else
        methLevelTop = -1.0, methLevelBottom = -1.0;

    skipLine(reader);
    return 0;
}

template <typename TReader>
inline bool
readCallRecord2(Dna5 &refAllele, Dna5String &genotype, double &methLevelTop, double &methLevelBottom,  unsigned &cov, TReader &reader)
{
    CharString helper;
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);        // ref
    refAllele = helper[0];
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qA
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qC
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qG
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qT
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qA
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qC
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qG
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // qT
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);        // cov
    cov = lexicalCast<unsigned>(helper);
    skipWhitespaces(reader);
    clear(genotype);
    readUntilWhitespace(genotype, reader);        // genotype called
    if (genotype == ".")
    {
        genotype[0] = refAllele;
        appendValue(genotype, refAllele);
    }
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);  // meth level called

    if (helper[0] != '.')
    {
        if (genotype[0] == 'C' && genotype[1] == 'G')
        {
            unsigned i = 0;
            while (i < length(helper) && helper[i] != ':')
                ++i;
            CharString pre = prefix(helper, i);
            methLevelTop = lexicalCast<double>(pre);
            CharString suff = suffix(helper, i+1);
            methLevelBottom = lexicalCast<double>(suff); 
        }
        else if (genotype[0] == 'C' || genotype[1] == 'C')
        {
            methLevelTop = lexicalCast<double>(helper);
            methLevelBottom = -1.0;
        }
        else if (genotype[0] == 'G' || genotype[1] == 'G')
        {
            methLevelTop = -1.0;
            methLevelBottom = lexicalCast<double>(helper);        
        }
    }
    else
        methLevelTop = -1.0, methLevelBottom = -1.0;

    skipLine(reader);
    return 0;
}


template<typename TPos, typename TReader>
inline void
readSimMethRecord1(CharString &contigName, TPos &pos, TReader &reader)
{
    clear(contigName);
    readUntilWhitespace(contigName, reader);  // contigName
    skipWhitespaces(reader);

    CharString posStr;
    readUntilWhitespace(posStr, reader);     // pos
    pos = lexicalCast<TPos>(posStr);
}

template<typename TReader>
inline void
readSimMethRecord2(double &methLevelTop, double &methLevelBottom, TReader &reader)
{
    skipWhitespaces(reader);
    CharString strand;
    readUntilWhitespace(strand, reader);     // +, -, +:-
    skipWhitespaces(reader);                
    CharString helper; 
    readUntilWhitespace(helper, reader);     // methylation level called
    if (strand == "+")
    {
        methLevelTop = lexicalCast<double>(helper);
        methLevelBottom = -1.0;
    }
    else if (strand == "-")
    {
        methLevelTop = -1.0;
        methLevelBottom = lexicalCast<double>(helper);
    }
    else if (strand == "+:-")
    {
        unsigned i = 0;
        while (i < length(helper) && helper[i] != ':')
            ++i;
        CharString pre = prefix(helper, i);
        methLevelTop = lexicalCast<double>(pre);
        CharString suff = suffix(helper, i+1);
        methLevelBottom = lexicalCast<double>(suff);
    }

    skipLine(reader);                        
}


template<typename TPos, typename TReader, typename TSnpFormatTag>
inline void
readSimSnpRecord1(CharString &contigName, TPos &pos, TReader &reader, TSnpFormatTag const &)
{
    clear(contigName);
    readUntilWhitespace(contigName, reader);  // contigName
    skipWhitespaces(reader);

    CharString posStr;
    readUntilWhitespace(posStr, reader);     // pos
    pos = lexicalCast<TPos>(posStr);
}

template<typename TReader>
inline bool
readSimSnpRecord2(Dna5 &refAllele, Dna5String &genotype, TReader &reader, SnpCustom const &)
{
    skipWhitespaces(reader);
    skipChar(reader, '.');
    skipWhitespaces(reader);
    CharString helper;
    readUntilWhitespace(helper, reader);     // ref
    if (length(helper) < 1)
    {
        return 1;
    }
    refAllele = helper[0];
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);     // genotype called
    resize(genotype, 2);
    if (length(helper) == 1)
    {
        genotype[0] = helper[0];
        skipWhitespaces(reader);
        skipChar(reader, '.');
        skipWhitespaces(reader);
        skipUntilWhitespace(reader);
        skipWhitespaces(reader);
        clear(helper);
        readUntilWhitespace(helper, reader);
        if (suffix(helper, 3) == "50")          // Create hetero genotype 
        {
            if (genotype[0] < refAllele)
                genotype[1] = refAllele;
            else
            {
                genotype[1] = genotype[0];
                genotype[0] = refAllele;
            }
        }
        else if (suffix(helper, 3) == "100")    // Create homo genotype
            genotype[1] = genotype[0];                   
    }
    else if (length(helper) == 3)
    {
        genotype[0] = helper[0];
        genotype[1] = helper[2];
        if (genotype[0] > genotype[1])
        {
            genotype[0] = helper[2];
            genotype[1] = helper[0];
        }
    }
    skipLine(reader);

    return 0;
}

template<typename TReader>
inline bool
readSimSnpRecord2(Dna5 &refAllele, Dna5String &genotype, TReader &reader, SnpMason2 const &)
{
    CharString helper;
    resize(genotype, 2);
    skipWhitespaces(reader);
    //skipChar(reader, '.');
    readUntilWhitespace(helper, reader); 
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);     // REF
    if (length(helper) > 1)
    {
        refAllele = 'N';
        genotype[0] = 'N';
        genotype[1] = 'N';
        skipLine(reader);
        return 0;
    }
    refAllele = helper[0];
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);     // genotype called ALT
    if ((length(helper) > 1 && helper[1] != ',') || length(helper) > 3)
    {
        refAllele = 'N';
        genotype[0] = 'N';
        genotype[1] = 'N';
        skipLine(reader);
        return 0;
    }
    if (length(helper) == 1)
    {
        genotype[0] = helper[0];
        skipWhitespaces(reader);
        skipChar(reader, '.');          // QUAL
        skipWhitespaces(reader);
        skipChar(reader, '.');    // FILTER
        skipWhitespaces(reader);
        skipChar(reader, '.');    // INFO
        skipWhitespaces(reader);
        skipChar(reader, '.');   // FORMAT
        skipWhitespaces(reader);
        clear(helper);
        readUntilWhitespace(helper, reader);    // simulated haplotype
        if (helper[0] == '0' || helper[2] == '0')         
        {
            if (genotype[0] < refAllele)
            {
                genotype[1] = refAllele;
            }
            else
            {
                genotype[1] = genotype[0];
                genotype[0] = refAllele;
            }
        }
        else 
        {
            genotype[1] = genotype[0];
        }                   
    }
    else if (length(helper) == 3)
    {
        genotype[0] = helper[0];
        genotype[1] = helper[2];
    }
    skipLine(reader);

    return 0;
}
    
template<typename TOptions, typename TSnpFormatTag>
bool benchmark(TOptions &options, TSnpFormatTag const &, MethCustom const &)
{
    std::cout << " Custom snp and meth format." << std::endl;
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;

    std::fstream simSnpFile(toCString(options.simSnpsFile), std::ios::binary | std::ios::in);
    TRecordReader simSnpReader(simSnpFile);

    std::fstream callFile(toCString(options.calledSnpsFile), std::ios::binary | std::ios::in);
    TRecordReader callReader(callFile);

    std::fstream simMethFile(toCString(options.simMethsFile), std::ios::binary | std::ios::in);
    TRecordReader simMethReader(simMethFile);

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

    readNChars(helper, simMethReader, 1);
    while(!atEnd(simMethReader) && helper[0] == '#')
    { 
        skipLine(simMethReader);
        clear(helper);
        readUntilWhitespace(helper, simMethReader);
    }
    simMethContig = helper;

   
    unsigned simSnpPos = 0;
    unsigned callPos = 0;
    unsigned simMethPos = 0;

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

    skipWhitespaces(simMethReader);
    readUntilWhitespace(helper, simMethReader);
    simMethPos = lexicalCast<size_t>(helper);

    if (atEnd(simSnpReader)) simSnpPos = maxPos;
    if (atEnd(simMethReader)) simMethPos = maxPos;
    if (atEnd(callReader)) callPos = maxPos;

    CharString currContig = callContig;
    if (atEnd(simSnpReader)) std::cout << " simSnpReader at end!" << std::endl;
    if (atEnd(callReader)) std::cout << " callReader at end!" << std::endl;
    if (atEnd(simMethReader)) std::cout << " simMethReader at end!" << std::endl;

    while(!atEnd(simSnpReader) || !atEnd(callReader) || !atEnd(simMethReader))
    {
        //if (callPos > 200) return 1;
        /*if (simSnpPos > 1000000 || (callPos >= 508900 && callPos < 508980))
        {
            std::cout << "ERROR: Should be at end, simSnpPos " << simSnpPos << " callPos: " << callPos << " simMethPos: " << simMethPos << "  maxPos: " << maxPos << std::endl;
            if (atEnd(simSnpReader)) std::cout << " simSnpReader at end!" << std::endl;
            if (atEnd(callReader)) std::cout << " callReader at end!" << std::endl;
            if (atEnd(simMethReader)) std::cout << " simMethReader at end!" << std::endl;
        }*/

        // Get highest position (documented)
        if (!atEnd(simSnpReader) && simSnpPos > maxPos) maxPos = simSnpPos +1;
        if (!atEnd(simMethReader) && simMethPos > maxPos) maxPos = simMethPos+1;
        if (!atEnd(callReader) && callPos > maxPos) maxPos = callPos+1;

        // If reader reached end of file, increase corresponding 'position' so that its not taking into account anymore
        if (atEnd(simSnpReader)) simSnpPos = maxPos;
        if (atEnd(simMethReader)) simMethPos = maxPos;
        if (atEnd(callReader)) callPos = maxPos;

        if (simSnpContig == callContig && simSnpContig == simMethContig)
        {
            currContig = callContig;
        }
        else
        {
            if (simSnpContig != currContig) simSnpPos = maxPos+1;   //??
            if (simMethContig != currContig) simMethPos = maxPos+1;
            if (callContig != currContig) callPos = maxPos+1;
        }

        //std::cout << " simSnpPos " << simSnpPos << " callPos: " << callPos << " simMethPos: " << simMethPos << "  maxPos: " << maxPos << std::endl;

        if (simSnpPos < callPos)     // Sim snp not called, iterate to next sim pos
        {
            ++c_NotListed;
            if (c_NotListed < 100)
            {
                //std::cout << "Not listed at simSnpPos: " << simSnpPos << std::endl;
            }
            ++c_SimSnps;

            skipLine(simSnpReader);                        // skip
            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader, TSnpFormatTag());
            continue;
        }
        if (simMethPos < callPos)     // Sim meth state not called, iterate to next sim meth pos
        {
            // ++?
            skipLine(simMethReader);                        // skip
            readSimMethRecord1(simMethContig, simMethPos, simMethReader);
            continue;
        }
        else if (simSnpPos > callPos && simMethPos > callPos)       // Something called (prob. a snp) but not simulated
        {
            Dna5 refAlleleCall;
            Dna5String genotypeCall;
            double methLevelTopCall;
            double methLevelBottomCall;

            if (readCallRecord2(refAlleleCall, genotypeCall, methLevelTopCall, methLevelBottomCall, callReader) == 1)  
            {
                std::cout << "ERROR: Something went wrong with call snp record at pos: " << " callPos: " << callPos  << std::endl;
                return 1;
            }

            
            if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
            {
                ++c_FalsePositive;
                ++c_CalledSnps;
                std::cout << "False positive at non-meth position, callPos: " << callPos << std::endl; 
            }

            if (readCallRecord1(callContig, callPos, callReader) == 1)
            {
                std::cout << "ERROR: Something went wrong with call snp record at pos: " << " callPos: " << callPos  << std::endl;
                return 1;
            }

            continue;
        }

        if (simSnpPos == callPos && callPos < simMethPos)     // Compare simulated snp and called snp only, no sim meth level given at that pos (no C or G in simulated genotype)
        {
            //std::cout << " Compare simulated snp and called snp only, no sim meth level given at that pos " << std::endl;
            // simSnpReader
            Dna5 refAlleleSimSnp;
            Dna5String genotypeSimSnp;
            if (readSimSnpRecord2(refAlleleSimSnp, genotypeSimSnp, simSnpReader, TSnpFormatTag()) == 1)
            {
                std::cout << "ERROR: Something went wrong with sim snp record at pos: " << simSnpPos << " callPos: " << callPos << " simMethPos: " << simMethPos << std::endl;
                return 1;
            }

            ++c_SimSnps;

            // callReader
            Dna5 refAlleleCall;
            Dna5String genotypeCall;
            double methLevelTopCall;
            double methLevelBottomCall;
            readCallRecord2(refAlleleCall, genotypeCall, methLevelTopCall, methLevelBottomCall, callReader);           

            SEQAN_ASSERT_EQ(refAlleleCall, refAlleleSimSnp);

            if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
            {
                ++c_CalledSnps;
                if (genotypeCall[0] != genotypeSimSnp[0] || genotypeCall[1] != genotypeSimSnp[1]) ++c_WrongCalled;
                else ++c_RightCalled;
            }
            else ++c_NotCalled;

            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader, TSnpFormatTag());
            readCallRecord1(callContig, callPos, callReader);
            continue;
        }
        else if (simSnpPos > callPos && simMethPos == callPos )    // compare simulated meth level and called meth level, no snp simulated at this position
        {
            //std::cout << " compare simulated meth level and called meth level, no snp simulated at this position" << std::endl;
            // simMethReader
            double methLevelTopSimMeth;
            double methLevelBottomSimMeth;
            readSimMethRecord2(methLevelTopSimMeth, methLevelBottomSimMeth, simMethReader);

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
            }
            else
            {
                if (methLevelTopSimMeth < -0.5) SEQAN_ASSERT_LEQ(methLevelTopCall, -0.5);
                if (methLevelBottomSimMeth < -0.5){
                    if (methLevelBottomCall > -0.5)
                    {
                        std::cout << "Pos: " << simMethPos << "  methLevelBottomSimMeth: " << methLevelBottomSimMeth << "  methLevelBottomCall: " << methLevelBottomCall << std::endl;
                    }
                    SEQAN_ASSERT_LEQ(methLevelBottomCall,-0.5);
                }
                if (methLevelTopCall < -0.5) SEQAN_ASSERT_LEQ(methLevelTopSimMeth, -0.5);
                if (methLevelBottomCall < -0.5) SEQAN_ASSERT_LEQ(methLevelBottomSimMeth, -0.5);
                if (methLevelTopSimMeth >= 0)
                {
                    double dist = std::abs(methLevelTopSimMeth - methLevelTopCall);

                    if      (dist < 0.000001) ++methDists[0];
                    else if (dist < 0.000005) ++methDists[1];
                    else if (dist < 0.00001) ++methDists[2];
                    else if (dist < 0.00005) ++methDists[3];
                    else if (dist < 0.0001) ++methDists[4];
                    else if (dist < 0.0005) ++methDists[5];
                    else if (dist < 0.001) ++methDists[6];
                    else if (dist < 0.005) ++methDists[7];
                    else if (dist < 0.01) ++methDists[8];
                    else if (dist < 0.05) ++methDists[9];
                    else if (dist < 0.1) ++methDists[10];
                    else if (dist < 0.2) ++methDists[11];
                    else if (dist < 0.3) ++methDists[12];
                    else if (dist < 0.4) ++methDists[13];
                    else if (dist < 0.5) ++methDists[14];
                    else if (dist < 0.6) ++methDists[15];
                    else if (dist < 0.7) ++methDists[16];
                    else if (dist < 0.8) ++methDists[17];
                    else if (dist < 0.9) ++methDists[18];
                    else
                    {
                        ++ methDists[9];
                        std::cout << "Dist: " << dist << "  methLevelTopSimMeth: " << methLevelTopSimMeth << " methLevelTopCall: " << methLevelTopCall << " pos: " <<  callPos <<  std::endl; 
                    }
                }
                if (methLevelBottomSimMeth >= 0)
                {
                    double dist = std::abs(methLevelBottomSimMeth - methLevelBottomCall);
                    if      (dist < 0.000001) ++methDists[0];
                    else if (dist < 0.000005) ++methDists[1];
                    else if (dist < 0.00001) ++methDists[2];
                    else if (dist < 0.00005) ++methDists[3];
                    else if (dist < 0.0001) ++methDists[4];
                    else if (dist < 0.0005) ++methDists[5];
                    else if (dist < 0.001) ++methDists[6];
                    else if (dist < 0.005) ++methDists[7];
                    else if (dist < 0.01) ++methDists[8];
                    else if (dist < 0.05) ++methDists[9];
                    else if (dist < 0.1) ++methDists[10];
                    else if (dist < 0.2) ++methDists[11];
                    else if (dist < 0.3) ++methDists[12];
                    else if (dist < 0.4) ++methDists[13];
                    else if (dist < 0.5) ++methDists[14];
                    else if (dist < 0.6) ++methDists[15];
                    else if (dist < 0.7) ++methDists[16];
                    else if (dist < 0.8) ++methDists[17];
                    else if (dist < 0.9) ++methDists[18];
                    else
                    {
                        ++ methDists[19];
                        std::cout << "Dist: " << dist << "  methLevelBottomSimMeth: " << methLevelBottomSimMeth << " methLevelBottomCall: " << methLevelBottomCall << " pos: " <<  callPos << std::endl; 
                    }
                }
            }

            readSimMethRecord1(simMethContig, simMethPos, simMethReader);
            //if (callPos >= 508953) std::cout << " 1 callPos: " << callPos << std::endl;
            readCallRecord1(callContig, callPos, callReader);

            //if (callPos >= 508953) std::cout << " 2 callPos: " << callPos << std::endl;
            continue;
        }
        else if (simSnpPos == callPos && simMethPos == callPos )    // compare simulated snp, simulated meth level and called snp and meth level
        {
            //std::cout << " test: compare simulated snp, simulated meth level and called snp and meth level " << std::endl; 
            // simSnpReader
            Dna5 refAlleleSimSnp;
            Dna5String genotypeSimSnp;
            if (readSimSnpRecord2(refAlleleSimSnp, genotypeSimSnp, simSnpReader, TSnpFormatTag()) == 1)
            {
                std::cout << "ERROR: Something went wrong with sim snp record at pos: " << simSnpPos << " callPos: " << callPos << " simMethPos: " << simMethPos << std::endl;
                return 1;
            }
            ++c_SimSnps;

            // simMethReader
            double methLevelTopSimMeth;
            double methLevelBottomSimMeth;
            readSimMethRecord2(methLevelTopSimMeth, methLevelBottomSimMeth, simMethReader);

            // callReader
            Dna5 refAlleleCall;
            Dna5String genotypeCall;
            double methLevelTopCall;
            double methLevelBottomCall;
            readCallRecord2(refAlleleCall, genotypeCall, methLevelTopCall, methLevelBottomCall, callReader);           

            CharString buffer;
            resize(buffer, 1024);
            sprintf(&buffer[0], "At position: %u", callPos);
            SEQAN_ASSERT_EQ_MSG(refAlleleCall, refAlleleSimSnp, toCString(buffer));

            if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
            {
                ++c_CalledSnps;
                if (genotypeCall[0] != genotypeSimSnp[0] || genotypeCall[1] != genotypeSimSnp[1]) ++c_WrongCalled;
                else
                {
                    ++c_RightCalled;
                    if (methLevelTopSimMeth < -0.5) SEQAN_ASSERT_LT_MSG(methLevelTopCall, -0.5, toCString(buffer));
                    if (methLevelBottomSimMeth < -0.5) SEQAN_ASSERT_LT_MSG(methLevelBottomCall, -0.5, toCString(buffer));
                    if (methLevelTopCall < -0.5) SEQAN_ASSERT_LT_MSG(methLevelTopSimMeth, -0.5, toCString(buffer));
                    if (methLevelBottomCall < -0.5) SEQAN_ASSERT_LT_MSG(methLevelBottomSimMeth, -0.5, toCString(buffer));
                    if (methLevelTopSimMeth >= 0)
                    {
                        double dist = abs(methLevelTopSimMeth - methLevelTopCall);
                        if      (dist < 0.000001) ++methDistsAtSnps[0];
                        else if (dist < 0.000005) ++methDistsAtSnps[1];
                        else if (dist < 0.00001) ++methDistsAtSnps[2];
                        else if (dist < 0.00005) ++methDistsAtSnps[3];
                        else if (dist < 0.0001) ++methDistsAtSnps[4];
                        else if (dist < 0.0005) ++methDistsAtSnps[5];
                        else if (dist < 0.001) ++methDistsAtSnps[6];
                        else if (dist < 0.005) ++methDistsAtSnps[7];
                        else if (dist < 0.01) ++methDistsAtSnps[8];
                        else if (dist < 0.05) ++methDistsAtSnps[9];

                        else if (dist < 0.1) ++methDistsAtSnps[10];
                        else if (dist < 0.2) ++methDistsAtSnps[11];
                        else if (dist < 0.3) ++methDistsAtSnps[12];
                        else if (dist < 0.4) ++methDistsAtSnps[13];
                        else if (dist < 0.5) ++methDistsAtSnps[14];
                        else if (dist < 0.6) ++methDistsAtSnps[15];
                        else if (dist < 0.7) ++methDistsAtSnps[16];
                        else if (dist < 0.8) ++methDistsAtSnps[17];
                        else if (dist < 0.9) ++methDistsAtSnps[18];
                        else ++ methDistsAtSnps[19];
                    }
                    if (methLevelBottomSimMeth >= 0)
                    {
                        double dist = abs(methLevelBottomSimMeth - methLevelBottomCall);
                        if (dist < 0.000001) ++methDistsAtSnps[0];
                        else if (dist < 0.000005) ++methDistsAtSnps[1];
                        else if (dist < 0.00001) ++methDistsAtSnps[2];
                        else if (dist < 0.00005) ++methDistsAtSnps[3];
                        else if (dist < 0.0001) ++methDistsAtSnps[4];
                        else if (dist < 0.0005) ++methDistsAtSnps[5];
                        else if (dist < 0.001) ++methDistsAtSnps[6];
                        else if (dist < 0.005) ++methDistsAtSnps[7];
                        else if (dist < 0.01) ++methDistsAtSnps[8];
                        else if (dist < 0.05) ++methDistsAtSnps[9];

                        else if (dist < 0.1) ++methDistsAtSnps[10];
                        else if (dist < 0.2) ++methDistsAtSnps[11];
                        else if (dist < 0.3) ++methDistsAtSnps[12];
                        else if (dist < 0.4) ++methDistsAtSnps[13];
                        else if (dist < 0.5) ++methDistsAtSnps[14];
                        else if (dist < 0.6) ++methDistsAtSnps[15];
                        else if (dist < 0.7) ++methDistsAtSnps[16];
                        else if (dist < 0.8) ++methDistsAtSnps[17];
                        else if (dist < 0.9) ++methDistsAtSnps[18];
                        else ++ methDistsAtSnps[19];
                   }
                }
            }
            else
            {
                ++c_NotCalled;
                std::cout << "Not called at pos: "  << " callPos: " << callPos <<  std::endl;
            }

            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader, TSnpFormatTag());
            readSimMethRecord1(simMethContig, simMethPos, simMethReader);
            readCallRecord1(callContig, callPos, callReader);
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
