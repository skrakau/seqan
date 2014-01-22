#ifndef SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_VS_BISSNP_H
#define SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_VS_BISSNP_H

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>


using namespace seqan;

template<typename TPos, typename TReader>
inline void
readBisSNPRecord1(CharString &contigName, TPos &pos, TReader &reader)
{
    clear(contigName);
    readUntilWhitespace(contigName, reader);  // contigName
    skipWhitespaces(reader);

    CharString posStr;
    readUntilWhitespace(posStr, reader);     // pos
    pos = lexicalCast<TPos>(posStr);
    --pos;      // TODO check 0, 1 based for all
}

template<typename TReader>
inline bool
readBisSNPRecord2(Dna5 &refAllele, Dna5String &genotype, TReader &reader)
{
    resize(genotype, 2);
    CharString helper;
    skipWhitespaces(reader);
    readUntilWhitespace(helper, reader);        // ID
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);        // REF
    refAllele = helper[0];
    skipWhitespaces(reader);
    clear(helper);
    readUntilWhitespace(helper, reader);     // genotype called ALT
    if (helper == '.' || length(helper) > 3)
    {
        std::cerr << "ERROR: Look at vcf entries of Bis-'SNP! No snp here, but something else!" << std::endl;
    }
    if (length(helper) == 1)
    {
        genotype[0] = helper[0];
        skipWhitespaces(reader);
        readUntilWhitespace(helper, reader);     // QUAL
        skipWhitespaces(reader);
        readUntilWhitespace(helper, reader);     // FILTER
        skipWhitespaces(reader);
        readUntilWhitespace(helper, reader);     // INFO
        skipWhitespaces(reader);
        readUntilWhitespace(helper, reader);    // FORMAT
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
        if ( helper[0] <  helper[2])
        {
            genotype[0] = helper[0];
            genotype[1] = helper[2];
        }
        else
        {
            genotype[0] = helper[2];
            genotype[1] = helper[0];
        }
    }
    skipLine(reader);

    return 0;
}


template<typename TOptions>
bool benchmark_bisSNP(TOptions &options)
{
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;

    std::fstream simSnpFile(toCString(options.simSnpsFile), std::ios::binary | std::ios::in);
    TRecordReader simSnpReader(simSnpFile);

    std::fstream callFile(toCString(options.calledSnpsFile), std::ios::binary | std::ios::in);
    TRecordReader callReader(callFile);

    // bisSnp
    std::fstream bisSNPFile(toCString(options.bisSnpFile), std::ios::binary | std::ios::in);
    TRecordReader bisSNPReader(bisSNPFile);

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
    CharString simSnpContig;
    clear(simSnpContig);
    CharString callContig;
    clear(callContig);
    CharString simMethContig;
    clear(simMethContig);
    CharString bisSNPContig;
    clear(bisSNPContig);

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
    readNChars(helper, bisSNPReader, 1);
    while(!atEnd(bisSNPReader) && helper[0] == '#')
    { 
        skipLine(bisSNPReader);
        clear(helper);
        readUntilWhitespace(helper, bisSNPReader);
    }
    bisSNPContig = helper;

   
    unsigned simSnpPos = 0;
    unsigned callPos = 0;
    unsigned bisSNPPos = 0;

    unsigned c_SimSnps = 0;     // TODO: struct stats
    
    unsigned c_CalledSnps = 0;
    unsigned c_NotListed = 0;    // Not found snps
    unsigned c_NotCalled = 0;    // 
    unsigned c_WrongCalled = 0;
    unsigned c_RightCalled = 0;
    unsigned c_FalsePositive = 0;    // False found snps
    unsigned c_FalsePositiveAtMeth = 0;    // False called snp at meth position

    unsigned c_bisSNPNotListed = 0;
    unsigned c_bisSNPNotCalled = 0;
    unsigned c_bisSNPWrongCalled = 0;
    unsigned c_bisSNPRightCalled = 0;
    unsigned c_bisSNPFalsePositive = 0;
    unsigned c_bisSNPFalsePositiveAtMeth = 0;

    unsigned c_BothNotListed = 0;
    unsigned c_BothNotCalled = 0;
    unsigned c_BothWrongCalled = 0;
    unsigned c_BothRightCalled = 0;
    unsigned c_BothFalsePositive = 0;
    unsigned c_BothFalsePositiveAtMeth = 0;

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
    skipWhitespaces(bisSNPReader);
    readUntilWhitespace(helper, bisSNPReader);
    bisSNPPos = lexicalCast<size_t>(helper);
    clear(helper);

    if (atEnd(simSnpReader)) simSnpPos = maxPos;
    if (atEnd(callReader)) callPos = maxPos;
    if (atEnd(bisSNPReader)) bisSNPPos = maxPos;


    CharString currContig = callContig;
    if (atEnd(simSnpReader)) std::cout << " simSnpReader at end!" << std::endl;
    if (atEnd(callReader)) std::cout << " callReader at end!" << std::endl;
    if (atEnd(bisSNPReader)) std::cout << " bisSNPReader at end!" << std::endl;
    
    if (atEnd(simSnpReader) || atEnd(callReader) || atEnd(bisSNPReader))
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
        if (!atEnd(bisSNPReader) && bisSNPPos >= maxPos) maxPos = bisSNPPos+1;


        // If reader reached end of file, increase corresponding 'position' so that its not taking into account anymore
        if (atEnd(simSnpReader)) simSnpPos = maxPos;
        if (atEnd(callReader)) callPos = maxPos;
        if (atEnd(bisSNPReader)) bisSNPPos = maxPos;

        if (simSnpContig == callContig && simSnpContig == simMethContig && bisSNPContig == simMethContig)
        {
            currContig = callContig;
        }
        else
        {
            if (simSnpContig != currContig) simSnpPos = maxPos+1;   //??
            if (callContig != currContig) callPos = maxPos+1;
            if (bisSNPContig != currContig) bisSNPPos = maxPos+1;
        }
        if (atEnd(simSnpReader)) std::cout << " simSnpReader at end!" << std::endl;
        if (atEnd(callReader)) std::cout << " callReader at end!" << std::endl;
        if (atEnd(bisSNPReader)) std::cout << " bisSNPReader at end!" << std::endl;

        if (simSnpPos < callPos && simSnpPos < bisSNPPos)     // Sim snp not called, iterate to next sim pos
        {
            Dna5 refAlleleSimSnp;
            Dna5String genotypeSimSnp;
            if (readSimSnpRecord2(refAlleleSimSnp, genotypeSimSnp, simSnpReader, SnpMason2()) == 1)
            {
                std::cout << "ERROR: b) Something went wrong with sim snp record at pos: " << simSnpPos << " callPos: " << callPos << std::endl;
                return 1;
            }
            if (refAlleleSimSnp !=  genotypeSimSnp[0] || refAlleleSimSnp !=  genotypeSimSnp[1])     // no snp simulated, something else
            {
                ++c_BothNotListed;
                ++c_NotListed;
                ++c_bisSNPNotListed;
                ++c_SimSnps;
            }
            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader, SnpMason2());
            continue;
        }
        if (bisSNPPos < callPos && bisSNPPos < simSnpPos)     // Sim snp not called, iterate to next sim pos
        {
            Dna5 refAlleleBisSNP;
            Dna5String genotypeBisSNP;
            readBisSNPRecord2(refAlleleBisSNP, genotypeBisSNP, bisSNPReader);           
           
            if (genotypeBisSNP[0] != refAlleleBisSNP || genotypeBisSNP[1] != refAlleleBisSNP)
            {
                ++c_bisSNPFalsePositive;
                if (refAlleleBisSNP == 'C' || refAlleleBisSNP == 'G') ++c_bisSNPFalsePositiveAtMeth;
            }


            if (bisSNPPos + 1 >= length(*itTop))  // end of contig simMethLevels, iter to next contig meths
            {
                ++itTop;
                ++itBottom;
                ++itContigName;
                if (itTop == itTopEnd) simMethsAtEnd = true;
                else simMethContig = (*itContigName);
            }
            readBisSNPRecord1(bisSNPContig, bisSNPPos, bisSNPReader);
            continue;
        }
        if (callPos < simSnpPos && callPos < bisSNPPos &&
                callPos < maxPos )    // compare simulated meth level and called meth level, no snp simulated at this position
        {
            double methLevelTopSimMeth = -1;
            double methLevelBottomSimMeth = -1;
            methLevelTopSimMeth = getMethLevel((*itTop)[callPos]);
            methLevelBottomSimMeth = getMethLevel((*itBottom)[callPos]); 
            Dna5 refAlleleCall;
            Dna5String genotypeCall;
            double methLevelTopCall;
            double methLevelBottomCall;
            readCallRecord2(refAlleleCall, genotypeCall, methLevelTopCall, methLevelBottomCall, callReader);           
           
            if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
            {
                ++c_FalsePositive;
                if (refAlleleCall == 'C' || refAlleleCall == 'G') ++c_FalsePositiveAtMeth;
                ++c_CalledSnps;
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
            continue;
        }
        if (callPos < simSnpPos && callPos == bisSNPPos && // called and bisSNP called, no snp simulated at this position
                callPos < maxPos )   
        {
            double methLevelTopSimMeth = -1;
            double methLevelBottomSimMeth = -1;
            methLevelTopSimMeth = getMethLevel((*itTop)[callPos]);
            methLevelBottomSimMeth = getMethLevel((*itBottom)[callPos]); 
            Dna5 refAlleleCall;
            Dna5 refAlleleBisSNP;
            Dna5String genotypeCall;
            Dna5String genotypeBisSNP;
            double methLevelTopCall;
            double methLevelBottomCall;
            readCallRecord2(refAlleleCall, genotypeCall, methLevelTopCall, methLevelBottomCall, callReader);           
            readBisSNPRecord2(refAlleleBisSNP, genotypeBisSNP, bisSNPReader); 
            if (refAlleleCall != refAlleleBisSNP) 
            {
                std::cerr << "Reference allele at call different than at bisSNP at pos: " << callPos << std::endl;
                return 1;
            }
            bool both = false;
            if (genotypeBisSNP[0] != refAlleleBisSNP || genotypeBisSNP[1] != refAlleleBisSNP)
            {
                ++c_bisSNPFalsePositive;
                if (refAlleleBisSNP == 'C' || refAlleleBisSNP == 'G') ++c_bisSNPFalsePositiveAtMeth;
                both = true;
            }
            if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
            {
                ++c_FalsePositive;
                if (refAlleleCall == 'C' || refAlleleCall == 'G') ++c_FalsePositiveAtMeth;
                if (both) ++c_BothFalsePositive;
                if (both && (refAlleleCall == 'C' || refAlleleCall == 'G')) ++c_BothFalsePositiveAtMeth;
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
            readBisSNPRecord1(bisSNPContig, bisSNPPos, bisSNPReader);
            continue;
        }
        else if (simSnpPos == callPos && callPos < bisSNPPos &&
                 simSnpPos < maxPos )    // compare simulated snp, simulated meth level and called snp and meth level
        {
            ++c_bisSNPNotListed;
            //std::cout << " test: compare simulated snp, simulated meth level and called snp and meth level " << std::endl; 
            // simSnpReader
            Dna5 refAlleleSimSnp;
            Dna5String genotypeSimSnp;
            if (readSimSnpRecord2(refAlleleSimSnp, genotypeSimSnp, simSnpReader, SnpMason2()) == 1)
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
                    if (refAlleleCall == 'C' || refAlleleCall == 'G') ++c_FalsePositiveAtMeth;
                    ++c_CalledSnps;
                }
                else 
                {
                    if (methLevelTopCall < -0.5 && methLevelTopSimMeth > 0.0)  std::cerr << "Simulated top meth level, but not called: " << callPos << std::endl;
                    if (methLevelBottomCall < -0.5 && methLevelBottomSimMeth > 0.0)  std::cerr << "Simulated bottom meth level but not called: " << callPos << std::endl;
                    if (methLevelTopCall >= 0) incrementMethDist(methDists, methLevelTopSimMeth, methLevelTopCall);
                    if (methLevelBottomCall >= 0) incrementMethDist(methDists, methLevelBottomSimMeth, methLevelBottomCall);
                }
            }

            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader, SnpMason2());
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
        else if (simSnpPos == bisSNPPos && simSnpPos < callPos &&
                 simSnpPos < maxPos )    
        {
            ++c_NotListed;
            //std::cout << " test: compare simulated snp, simulated meth level and called snp and meth level " << std::endl; 
            // simSnpReader
            Dna5 refAlleleSimSnp;
            Dna5String genotypeSimSnp;
            if (readSimSnpRecord2(refAlleleSimSnp, genotypeSimSnp, simSnpReader, SnpMason2()) == 1)
            {
                std::cout << "ERROR: b) Something went wrong with sim snp record at pos: " << simSnpPos << " callPos: " << callPos << std::endl;
                return 1;
            }

            Dna5 refAlleleBisSNP;
            Dna5String genotypeBisSNP;
            readBisSNPRecord2(refAlleleBisSNP, genotypeBisSNP, bisSNPReader);           

            if (refAlleleSimSnp !=  genotypeSimSnp[0] || refAlleleSimSnp !=  genotypeSimSnp[1])     // no snp simulated, something else
            {
                ++c_SimSnps;
                CharString buffer;
                resize(buffer, 1024);
                sprintf(&buffer[0], "At position: %u", bisSNPPos);
                SEQAN_ASSERT_EQ_MSG(refAlleleBisSNP, refAlleleSimSnp, toCString(buffer));

                if (genotypeBisSNP[0] != refAlleleBisSNP || genotypeBisSNP[1] != refAlleleBisSNP)   // if no snp, field would be '.' ?
                {
                    if (genotypeBisSNP[0] != genotypeSimSnp[0] || genotypeBisSNP[1] != genotypeSimSnp[1])
                    {
                        ++c_bisSNPWrongCalled;
                    }
                    else
                    {
                        ++c_bisSNPRightCalled;
                    }
                }
                else ++c_bisSNPNotCalled;
            }
            else
            {
                if (genotypeBisSNP[0] != refAlleleBisSNP || genotypeBisSNP[1] != refAlleleBisSNP)
                {
                    ++c_bisSNPFalsePositive;
                    if (refAlleleBisSNP == 'C' || refAlleleBisSNP == 'G') ++c_bisSNPFalsePositiveAtMeth;
                }
            }

            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader,  SnpMason2());
            readBisSNPRecord1(bisSNPContig, bisSNPPos, bisSNPReader);
            if (bisSNPPos + 1 >= length(*itTop))  // end of contig simMethLevels, iter to next contig meths
            {
                ++itTop;
                ++itBottom;
                ++itContigName;
                if (itTop == itTopEnd) simMethsAtEnd = true;
                else simMethContig = (*itContigName);
            }
            continue;
        } 
        else if (simSnpPos == callPos && callPos == bisSNPPos &&
                 simSnpPos < maxPos )    // compare simulated snp, simulated meth level and called snp and meth level
        {
            //std::cout << " test: compare simulated snp, simulated meth level and called snp and meth level " << std::endl; 
            // simSnpReader
            Dna5 refAlleleSimSnp;
            Dna5String genotypeSimSnp;
            if (readSimSnpRecord2(refAlleleSimSnp, genotypeSimSnp, simSnpReader,  SnpMason2()) == 1)
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
            Dna5 refAlleleBisSNP;
            Dna5String genotypeCall;
            Dna5String genotypeBisSNP;

            double methLevelTopCall;
            double methLevelBottomCall;
            readCallRecord2(refAlleleCall, genotypeCall, methLevelTopCall, methLevelBottomCall, callReader);           
            readBisSNPRecord2(refAlleleBisSNP, genotypeBisSNP, bisSNPReader); 

            if (refAlleleSimSnp !=  genotypeSimSnp[0] || refAlleleSimSnp !=  genotypeSimSnp[1])     // no snp simulated, something else
            {
                ++c_SimSnps;
                CharString buffer;
                resize(buffer, 1024);
                sprintf(&buffer[0], "At position: %u", callPos);
                SEQAN_ASSERT_EQ_MSG(refAlleleCall, refAlleleSimSnp, toCString(buffer));
                // Check my call
                bool bothRight = false;
                bool bothWrong = false;
                if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
                {
                    ++c_CalledSnps;
                    if (genotypeCall[0] != genotypeSimSnp[0] || genotypeCall[1] != genotypeSimSnp[1])
                    {
                        ++c_WrongCalled;
                        bothWrong = true;
                    }
                    else
                    {
                        ++c_RightCalled;
                        bothRight = true;
                       if (methLevelTopCall < -0.5 && methLevelTopSimMeth > 0.0)  std::cerr << "Simulated top meth level, but not called: " << callPos << std::endl;
                       if (methLevelBottomCall < -0.5 && methLevelBottomSimMeth > 0.0)  std::cerr << "Simulated bottom meth level but not called: " << callPos << std::endl;
                       if (methLevelTopCall >= 0) incrementMethDist(methDistsAtSnps, methLevelTopSimMeth, methLevelTopCall);
                       if (methLevelBottomCall >= 0) incrementMethDist(methDistsAtSnps, methLevelBottomSimMeth, methLevelBottomCall);
                    }
                }
                else ++c_NotCalled;
                // Check bisSNP call
                if (genotypeBisSNP[0] != refAlleleBisSNP || genotypeBisSNP[1] != refAlleleBisSNP)
                {
                    if (genotypeBisSNP[0] != genotypeSimSnp[0] || genotypeBisSNP[1] != genotypeSimSnp[1])
                    {
                        ++c_bisSNPWrongCalled;
                        if (bothWrong) ++c_BothWrongCalled;
                    }
                    else
                    {
                        ++c_bisSNPRightCalled;
                        if (bothRight) ++c_BothRightCalled;
                    }
                }
                else ++c_bisSNPNotCalled;
            }
            else
            {
                if (refAlleleCall != refAlleleBisSNP) 
                {
                    std::cerr << "Reference allele at call different than at bisSNP at pos: " << callPos << std::endl;
                    return 1;
                }
                bool both = false;
                if (genotypeBisSNP[0] != refAlleleBisSNP || genotypeBisSNP[1] != refAlleleBisSNP)
                {
                    ++c_bisSNPFalsePositive;
                    if (refAlleleBisSNP == 'C' || refAlleleBisSNP == 'G') ++c_bisSNPFalsePositiveAtMeth;
                    both = true;
                }
                if (genotypeCall[0] != refAlleleCall || genotypeCall[1] != refAlleleCall)
                {
                    ++c_FalsePositive;
                    if (refAlleleCall == 'C' || refAlleleCall == 'G') ++c_FalsePositiveAtMeth;
                    if (both) ++c_BothFalsePositive;
                    if (both && (refAlleleCall == 'C' || refAlleleCall == 'G')) ++c_BothFalsePositiveAtMeth;
                }
                else 
                {
                    if (methLevelTopCall < -0.5 && methLevelTopSimMeth > 0.0)  std::cerr << "Simulated top meth level, but not called: " << callPos << std::endl;
                    if (methLevelBottomCall < -0.5 && methLevelBottomSimMeth > 0.0)  std::cerr << "Simulated bottom meth level but not called: " << callPos << std::endl;
                    if (methLevelTopCall >= 0) incrementMethDist(methDists, methLevelTopSimMeth, methLevelTopCall);
                    if (methLevelBottomCall >= 0) incrementMethDist(methDists, methLevelBottomSimMeth, methLevelBottomCall);
                }
            }

            readSimSnpRecord1(simSnpContig, simSnpPos, simSnpReader,  SnpMason2());
            readCallRecord1(callContig, callPos, callReader);
            readBisSNPRecord1(bisSNPContig, bisSNPPos, bisSNPReader);
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

    out << "Bis-SNP:" << '\n';
    unsigned bisSNP_trueNegatives = maxPos - c_bisSNPFalsePositive;  
    out<< "Specificity: "  << (double)bisSNP_trueNegatives/(double)(bisSNP_trueNegatives + c_bisSNPFalsePositive) <<  '\n';
    out<< "Recall: " << (double)c_bisSNPRightCalled/(double)(c_bisSNPRightCalled + c_bisSNPWrongCalled + c_bisSNPNotCalled + c_bisSNPNotListed - - options.simulatedAtRefN) << '\n';
    out<< "Precision: "  << (double)c_bisSNPRightCalled /(double)(c_bisSNPRightCalled + c_bisSNPFalsePositive + c_bisSNPWrongCalled) <<  '\n';

    out<< '\n';
     out << "Not listed: \t" << c_bisSNPNotListed << '\n';
     out << "Not called: \t" << c_bisSNPNotCalled << '\n';
     out << "False positive: \t" << c_bisSNPFalsePositive << '\n';
     out << "False positive at meth positions: \t" << c_bisSNPFalsePositiveAtMeth << '\n';
     out << "Wrong called genotype: \t" << c_bisSNPWrongCalled << '\n'; 
     out << "Right called genotype: \t" << c_bisSNPRightCalled << '\n';
     out << '\n';

    out<< '\n';
    out << "Me:" << '\n';
    unsigned trueNegatives = maxPos - c_FalsePositive; 
    out<< "Specificity: "  << (double)trueNegatives/(double)(trueNegatives + c_FalsePositive) <<  '\n';
    out<< "Recall: " << (double)c_RightCalled/(double)(c_RightCalled + c_WrongCalled + c_NotCalled + c_NotListed) << '\n';
    out<< "Precision: "  << (double)c_RightCalled /(double)(c_RightCalled + c_FalsePositive + c_WrongCalled) <<  '\n';

    out << "Simulated snps: \t" << c_SimSnps << '\n';
    out<< "Called snps: \t" << c_CalledSnps << '\n';
    out<< '\n';
    out<< "Not listed: \t" << c_NotListed << '\n';
    out<< "Not called: \t" << c_NotCalled << '\n';
    out<< "False positive: \t" << c_FalsePositive << '\n';
    out<< "False positive at meth positions: \t" << c_FalsePositiveAtMeth << '\n';
    out<< "Wrong called genotype: \t" << c_WrongCalled << '\n'; 
    out<< "Right called genotype: \t" << c_RightCalled << '\n';

    out<< '\n';
    out<< '\n';

     out << '\n';
     out << "Both:" << '\n';
     out << "Not listed: \t" << c_BothNotListed << '\n';
     out << "Not called: \t" << c_BothNotCalled << '\n';
     out << "False positive: \t" << c_BothFalsePositive << '\n';
     out << "Wrong called genotype: \t" << c_BothWrongCalled << '\n'; 
     out << "Right called genotype: \t" << c_BothRightCalled << '\n';

    
    std::cout << " maxPos: " << maxPos << std::endl;

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
