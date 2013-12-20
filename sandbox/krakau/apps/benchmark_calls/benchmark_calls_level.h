#ifndef SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_LEVEL_H
#define SANDBOX_KRAKAU_APPS_BENCHMARK_CALLS_BENCHMARK_CALLS_LEVEL_H

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>


using namespace seqan;


// called snps vs. simulated snps
template<typename TOptions>
bool mainBenchmark_level(TOptions &options)
{
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
   
    CharString simSnpPos;
    CharString callPos;
    CharString simMethPos;

    char simMethO;
    CharString simMethLevel;
    double callMethLevel;
    CharString callCoverage;
    
    Dna5 simRef;
    Dna5 callRef;

    Dna5 simA1;
    Dna5 simA2;
    Dna5 callA1;
    Dna5 callA2;

    unsigned c_SimSnps = 0;
    unsigned c_CalledSnps = 0;
    unsigned c_NotListed = 0;    // Not found snps
    unsigned c_NotCalled = 0;    // 
    unsigned c_WrongCalled = 0;
    unsigned c_RightCalled = 0;
    unsigned c_FalsePositive = 0;    // False found snps
    unsigned c_FalsePositiveAtMeth = 0;    // False called snp at meth position
    unsigned c_Fuzzy = 0;       // at the moment: just positions without snp
    unsigned c_FPFuzzy = 0;



    unsigned c_rightMeths = 0;
    unsigned c_falseMeths = 0;
    unsigned c_weirdSnpMeths = 0;
    
    unsigned maxPos = 0;
    double methThres = 0.2;

    clear(simSnpPos);
    skipWhitespaces(simSnpReader);
    readUntilWhitespace(simSnpPos, simSnpReader);

    clear(callPos);
    skipWhitespaces(callReader);
    readUntilWhitespace(callPos, callReader);

    clear(simMethPos);
    skipWhitespaces(simMethReader);
    readUntilWhitespace(simMethPos, simMethReader);

    CharString currContig = simSnpContig;
    // get rid of record reader, doesnt make sense, only for sequence files etc...
    while(!atEnd(simSnpReader) || !atEnd(callReader) || !atEnd(simMethReader))
    {
        // Get highest position (documented)
        if (!atEnd(simSnpReader) && lexicalCast<size_t>(simSnpPos) > maxPos) maxPos = lexicalCast<size_t>(simSnpPos);
        if (!atEnd(simMethReader) && lexicalCast<size_t>(simMethPos) > maxPos) maxPos = lexicalCast<size_t>(simMethPos);
        if (!atEnd(callReader) && lexicalCast<size_t>(callPos) > maxPos) maxPos = lexicalCast<size_t>(callPos);

        // If reader reached end of file, increase corresponding 'position' so that its not taking into account anymore
        CharString buffer;
        resize(buffer, 1024);
        sprintf(&buffer[0], "%u", maxPos+1);
        if (atEnd(simSnpReader)) simSnpPos = buffer;
        if (atEnd(simMethReader)) simMethPos = buffer;
        if (atEnd(callReader)) callPos = buffer;

        
        if (simSnpContig == callContig && simSnpContig == simMethContig)
        {
            currContig = callContig;
        }
        else
        {
            if (simSnpContig != currContig) simSnpPos = buffer;
            if (simMethContig != currContig) simMethPos = buffer;
            if (callContig != currContig) callPos = buffer;
        }
       

        if (lexicalCast<size_t>(simSnpPos) < lexicalCast<size_t>(callPos))     // Sim snp not called, iterate to next sim pos
        {
            ++c_NotListed;
            if (c_NotListed < 100)
            {
                //std::cout << "Not listed at simSnpPos: " << simSnpPos << std::endl;
            }
            ++c_SimSnps;

            skipLine(simSnpReader);                        // skip
            clear(simSnpContig);
            readUntilWhitespace(simSnpContig, simSnpReader);  // chrom
            skipWhitespaces(simSnpReader);
            clear(simSnpPos);
            readUntilWhitespace(simSnpPos, simSnpReader);     // pos
            continue;
        }
        if (lexicalCast<size_t>(simMethPos) < lexicalCast<size_t>(callPos))     // Sim meth state not called, iterate to next sim meth pos
        {
            // ++?

            skipLine(simMethReader);                        // skip
            clear(simMethContig);
            readUntilWhitespace(simMethContig, simMethReader);  // chrom
            skipWhitespaces(simMethReader);
            clear(simMethPos);
            readUntilWhitespace(simMethPos, simMethReader);     // pos
            continue;
        }
        else if (lexicalCast<size_t>(simSnpPos) > lexicalCast<size_t>(callPos) && lexicalCast<size_t>(simMethPos) > lexicalCast<size_t>(callPos))       // Something called (prob. a snp) but not simulated
        {
            //std::cout << " Something called (prob. a snp) but not simulated " << std::endl; 
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // ref
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qA
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qC
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qG
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qT
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qA
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qC
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qG
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qT
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // cov
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // genotype called
            if (helper != '.')
            {
                callA1 = helper[0];
                callA2 = helper[1];
                ++c_FalsePositive;
                ++c_CalledSnps;
                // ++ Meth?
                std::cout << "False positive at non-meth position, callPos: " << callPos << std::endl; 
            }
            // TODO: maybe snp and then methlevel ? check?

            skipLine(callReader);                           // skip
            clear(callContig);
            readUntilWhitespace(callContig, callReader);    // chrom
            skipWhitespaces(callReader);
            clear(callPos);
            readUntilWhitespace(callPos, callReader);       // pos
            continue;
        }

        if (lexicalCast<size_t>(simSnpPos) == lexicalCast<size_t>(callPos) && lexicalCast<size_t>(callPos) < lexicalCast<size_t>(simMethPos))     // Compare simulated snp and called snp only, no sim meth level given at that pos
        {
            //std::cout << " Compare simulated snp and called snp only, no sim meth level given at that pos " << std::endl;
            // simSnpReader
            skipWhitespaces(simSnpReader);
            skipChar(simSnpReader, '.');
            skipWhitespaces(simSnpReader);
            clear(helper);
            readUntilWhitespace(helper, simSnpReader);     // ref
            simRef = helper[0];
            skipWhitespaces(simSnpReader);
            clear(helper);
            readUntilWhitespace(helper, simSnpReader);     // genotype called
            if (length(helper) == 1)
            {
                simA1 = helper[0];
                // Check if hetero or homo snp
                skipWhitespaces(simSnpReader);
                skipChar(simSnpReader, '.');
                skipWhitespaces(simSnpReader);
                skipUntilWhitespace(simSnpReader);
                skipWhitespaces(simSnpReader);
                clear(helper);
                readUntilWhitespace(helper, simSnpReader);
                if (prefix(helper, 3) == "AF=")
                {
                    if (suffix(helper, 3) == "50")          // Create hetero genotype 
                    {
                        if (simA1 < simRef)
                            simA2 = simRef;
                        else
                        {
                            simA2 = simA1;
                            simA1 = simRef;
                        }
                    }
                    else if (suffix(helper, 3) == "100")    // Create homo genotype
                        simA2 = simA1;                   
                    else
                    {
                        std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simSnpPos << std::endl;
                        return 1;
                    }
                }
                else
                {
                    std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simSnpPos << std::endl;
                    return 1;
                }
            }
            else if (length(helper) == 3)
            {
                simA1 = helper[0];
                simA2 = helper[2];
                if (simA1 > simA2)
                {
                    simA1 = helper[2];
                    simA2 = helper[0];
                }
            }
            else if (length(helper) != 0)
            {
                std::cerr << "ERROR: Something wrong in simulated SNPs file with ALT entry at position: " << simSnpPos << std::endl;
                return 1;
            }

            ++c_SimSnps;
            skipLine(simSnpReader);                        // skip
            clear(simSnpContig);
            readUntilWhitespace(simSnpContig, simSnpReader);  // chrom
            skipWhitespaces(simSnpReader);
            clear(simSnpPos);
            readUntilWhitespace(simSnpPos, simSnpReader);     // pos

            // callReader
            skipWhitespaces(callReader);
            clear(helper);
             readUntilWhitespace(helper, callReader);        // ref
            Dna5 ref = helper[0];
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qA
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qC
            unsigned countC_F = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qG
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qT
            unsigned countT_F = length(helper)-2;
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qA  (T on reverse)
            unsigned countT_R = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qC
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qG  (C on reverse)
            unsigned countC_R = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qT
            skipWhitespaces(callReader);
            clear(callCoverage);
            readUntilWhitespace(callCoverage, callReader);        // cov
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // genotype called 
            bool snpCalled = false;
            if (helper != '.')
            {
                callA1 = helper[0];
                callA2 = helper[1];
                ++c_CalledSnps;
                snpCalled = true;
                // compare called genotypes
                if (simA1 != callA1 || simA2 != callA2)
                {
                    ++c_WrongCalled;
                    if (options.verbosity >= 4)
                    {
                        std::cout << "Wrong called at position: " << callPos << std::endl;
                        std::cout << "simA1: " << simA1 << " simA2: " << simA2 << " callA1: " << callA1 << " callA2: " << callA2 << std::endl;
                    }
                }
                else
                {
                    ++c_RightCalled;
                   // std::cout << "Right called at callPos: " << callPos << std::endl;
                }
            }
            else
            {
                ++c_NotCalled;
                if (options.verbosity >= 4)
                {
                    std::cout << "Not called at position: " << callPos << " coverage: " << callCoverage << std::endl;
                }
            }
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);  // meth level called

            if (helper[0] != '.')
            {
                if (!snpCalled)
                {
                    callA1 = ref;
                    callA2 = ref;
                }
                // if genotype CG -> methLevel1:methLevel2
                if (callA1 == 'C' && callA2 == 'G')
                {
                    unsigned i = 0;
                    while (i < length(helper) && helper[i] != ':')
                        ++i;
                    std::cout << std::endl;

                    CharString pre = prefix(helper, i);
                    callMethLevel = lexicalCast<double>(pre); 
                    if (countC_F + countT_F > 0)
                    {
                        double guess = (double)countC_F / (double)(countC_F + countT_F);
                        // Check methlevel roughly and count if too different               // We don't know whats right and wrong, just to check if weird things happen!!!
                        if (abs(guess - callMethLevel) > methThres)
                        {
                            ++c_falseMeths;
                            std::cout << "     Snp: Weird meth called at pos: " << callPos << std::endl;
                        }
                       if (abs(guess - callMethLevel) < methThres)
                            ++c_rightMeths;
                    }
                    CharString suff = suffix(helper, i+1);
                    callMethLevel = lexicalCast<double>(suff); 
                    if (countC_R + countT_R > 0)
                    {
                        double guess = (double)countC_R / (double)(countC_R + countT_R);
                        // Check methlevel roughly and count if too different
                        if (abs(guess - callMethLevel) > methThres)
                        {
                            ++c_falseMeths;
                            std::cout << "     Snp: Weird meth called at pos: " << callPos << std::endl;
                        }
                       if (abs(guess - callMethLevel) < methThres)
                            ++c_rightMeths;
                    }
                }
                else
                {   
                    if (callA1 == 'C' || callA2 == 'C')
                    {
                        callMethLevel = lexicalCast<double>(helper); 
                        if (countC_F + countT_F > 0)
                        {
                            double guess = (double)countC_F / (double)(countC_F + countT_F);
                            // Check methlevel roughly and count if too different               // We don't know whats right and wrong, just to check if weird things happen!!!
                            if (abs(guess - callMethLevel) > methThres)
                             {
                                ++c_falseMeths;
                                std::cout << "     Snp: Weird meth called at pos: " << callPos << std::endl;
                            }
                           if (abs(guess - callMethLevel) < methThres)
                                ++c_rightMeths;
                        }
                    }
                    else if (callA1 == 'G' || callA2 == 'G')
                    {
                        callMethLevel = lexicalCast<double>(helper); 
                        if (countC_R + countT_R > 0)
                        {
                            double guess = (double)countC_R / (double)(countC_R + countT_R);
                            // Check methlevel roughly and count if too different
                            if (abs(guess - callMethLevel) > methThres)
                            {
                                ++c_falseMeths;
                                std::cout << "     Snp: Weird meth called at pos: " << callPos << std::endl;
                            }
                           if (abs(guess - callMethLevel) < methThres)
                                ++c_rightMeths;
                        }
                    }
                }
            }
            skipWhitespaces(callReader);
            readUntilOneOf(helper, callReader, ':');    // genotype prob
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);    // genotype score

            skipLine(callReader);
            clear(callContig);
            readUntilWhitespace(callContig, callReader);    // chrom
            skipWhitespaces(callReader);
            clear(callPos);
            readUntilWhitespace(callPos, callReader);       // pos

            continue;
        }
        else if (lexicalCast<size_t>(simSnpPos) > lexicalCast<size_t>(callPos) && lexicalCast<size_t>(simMethPos) == lexicalCast<size_t>(callPos) )    // compare simulated meth level and called meth level, no snp simulated at this position
        {
            //std::cout << " compare simulated meth level and called meth level, no snp simulated at this position" << std::endl;
            // simMethReader
            skipWhitespaces(simMethReader);
            clear(helper);
            readUntilWhitespace(helper, simMethReader);          // strand, '+' or '-' -> Use it to check?
            simMethO = helper[0];
            skipWhitespaces(simMethReader);                             
            readUntilWhitespace(simMethLevel, simMethReader);          // methylation level simulated

            skipLine(simMethReader);                        // skip
            clear(simMethContig);
            readUntilWhitespace(simMethContig, simMethReader);  // chrom
            skipWhitespaces(simMethReader);
            clear(simMethPos);
            readUntilWhitespace(simMethPos, simMethReader);     // pos

            // callReader
            skipWhitespaces(callReader);
            clear(helper);
             readUntilWhitespace(helper, callReader);        // ref
            Dna5 ref = helper[0];
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qA
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qC
            unsigned countC_F = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qG
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qT
            unsigned countT_F = length(helper)-2;
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qA  (T on reverse)
            unsigned countT_R = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qC
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qG  (C on reverse)
            unsigned countC_R = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qT
            skipWhitespaces(callReader);
            clear(callCoverage);
            readUntilWhitespace(callCoverage, callReader);        // cov
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // genotype called
            bool snpCalled = false;
            double guessC = (double)countC_F / (double)(countC_F + countT_F); 
            double guessG = (double)countC_R / (double)(countC_R + countT_R);
            if (ref == 'C' && guessC < 0.7 && guessC > 0.3)
                ++c_Fuzzy;
            else if (ref == 'G' && guessG < 0.7 && guessG > 0.3)
                ++c_Fuzzy;

            if (helper != '.')
            {
                callA1 = helper[0];
                callA2 = helper[1];
                ++c_FalsePositive;
                std::cout << "False positive at meth position, callPos: " << callPos << std::endl;

                if (ref == 'C' && guessC < 0.7 && guessC > 0.3)
                    ++c_FPFuzzy;
                else if (ref == 'G' && guessG < 0.7 && guessG > 0.3)
                    ++c_FPFuzzy;

                ++c_FalsePositiveAtMeth;
                ++c_CalledSnps;
                snpCalled = true;
            }

            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);  // meth level called
            if (helper[0] != '.')
            {
                if (!snpCalled)
                {
                    callA1 = ref;
                    callA2 = ref;
                }
                // if genotype CG -> methLevel1:methLevel2
                if (callA1 == 'C' && callA2 == 'G')
                {
                    unsigned i = 0;
                    while (i < length(helper) && helper[i] != ':')
                        ++i;
                    CharString pre = prefix(helper, i);
                    callMethLevel = lexicalCast<double>(pre);

                    if (countC_F + countT_F > 0)
                    {
                        double guess = (double)countC_F / (double)(countC_F + countT_F);
                        // Check methlevel roughly and count if too different               // We don't know whats right and wrong, just to check if weird things happen!!!
                        if (abs(guess - callMethLevel) > methThres)
                        {
                            ++c_falseMeths;
                            std::cout << "    No Snp simulated: Weird meth called at pos: " << callPos << std::endl;
                        }
                       if (abs(guess - callMethLevel) < methThres)
                            ++c_rightMeths;
                    }
                    CharString suff = suffix(helper, i+1);
                    callMethLevel = lexicalCast<double>(suff); 

                    if (countC_R + countT_R > 0)
                    {
                        double guess = (double)countC_R / (double)(countC_R + countT_R);
                        // Check methlevel roughly and count if too different
                        if (abs(guess - callMethLevel) > methThres)
                        {
                            ++c_falseMeths;
                            std::cout << "     No Snp simulated: Weird meth called at pos: " << callPos << std::endl;
                        }
                       if (abs(guess - callMethLevel) < methThres)
                            ++c_rightMeths;
                    }
                }
                else
                {   
                    if (callA1 == 'C' || callA2 == 'C')
                    {
                        callMethLevel = lexicalCast<double>(helper); 
                        if (countC_F + countT_F > 0)
                        {
                            double guess = (double)countC_F / (double)(countC_F + countT_F);
                            // Check methlevel roughly and count if too different               // We don't know whats right and wrong, just to check if weird things happen!!!
                            if (abs(guess - callMethLevel) > methThres)
                            {
                                ++c_falseMeths;
                                std::cout << "     No Snp simulated: Weird meth called at pos: " << callPos << std::endl;
                            }
                            if (abs(guess - callMethLevel) < methThres)
                                ++c_rightMeths;
                        }
                    }
                    else if (callA1 == 'G' || callA2 == 'G')
                    {
                        callMethLevel = lexicalCast<double>(helper); 
                        if (countC_R + countT_R > 0)
                        {
                            double guess = (double)countC_R / (double)(countC_R + countT_R);
                            // Check methlevel roughly and count if too different
                            if (abs(guess - callMethLevel) > methThres)
                             {
                                ++c_falseMeths;
                                std::cout << "     No Snp simulated: Weird meth called at pos: " << callPos << std::endl;
                            }
                           if (abs(guess - callMethLevel) < methThres)
                                ++c_rightMeths;
                        }
                    }
                }
            }
            skipWhitespaces(callReader);
            readUntilOneOf(helper, callReader, ':');    // genotype prob
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);    // genotype score

            skipLine(callReader);
            clear(callContig);
            readUntilWhitespace(callContig, callReader);    // chrom
            skipWhitespaces(callReader);
            clear(callPos);
            readUntilWhitespace(callPos, callReader);       // pos

            continue;
        }
        else if (lexicalCast<size_t>(simSnpPos) == lexicalCast<size_t>(callPos) && lexicalCast<size_t>(simMethPos) == lexicalCast<size_t>(callPos) )    // compare simulated snp, simulated meth level and called snp and meth level
        {
            //std::cout << " test: compare simulated snp, simulated meth level and called snp and meth level " << std::endl; 
            // simSnpReader
            skipWhitespaces(simSnpReader);
            skipChar(simSnpReader, '.');
            skipWhitespaces(simSnpReader);
            clear(helper);
            readUntilWhitespace(helper, simSnpReader);     // ref
            simRef = helper[0];
            skipWhitespaces(simSnpReader);
            clear(helper);
            readUntilWhitespace(helper, simSnpReader);     // genotype called
            if (length(helper) == 1)
            {
                simA1 = helper[0];
                // Check if hetero or homo snp
                skipWhitespaces(simSnpReader);
                skipChar(simSnpReader, '.');
                skipWhitespaces(simSnpReader);
                skipUntilWhitespace(simSnpReader);
                skipWhitespaces(simSnpReader);
                clear(helper);
                readUntilWhitespace(helper, simSnpReader);
                if (prefix(helper, 3) == "AF=")
                {
                    if (suffix(helper, 3) == "50")          // Create hetero genotype 
                    {
                        if (simA1 < simRef)
                            simA2 = simRef;
                        else
                        {
                            simA2 = simA1;
                            simA1 = simRef;
                        }
                    }
                    else if (suffix(helper, 3) == "100")    // Create homo genotype
                        simA2 = simA1;                   
                    else
                    {
                        std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simSnpPos << std::endl;
                        return 1;
                    }
                }
                else
                {
                    std::cerr << "ERROR: Something wrong in simulated SNPs file with INFO entry at position: " << simSnpPos << std::endl;
                    return 1;
                }
            }
            else if (length(helper) == 3)
            {
                simA1 = helper[0];
                simA2 = helper[2];
                if (simA1 > simA2)
                {
                    simA1 = helper[2];
                    simA2 = helper[0];
                }
            }
            else if (length(helper) != 0)
            {
                std::cerr << "ERROR: Something wrong in simulated SNPs file with ALT entry at position: " << simSnpPos << std::endl;
                return 1;
            }

            ++c_SimSnps;
            skipLine(simSnpReader);                        // skip
            clear(simSnpContig);
            readUntilWhitespace(simSnpContig, simSnpReader);  // chrom
            skipWhitespaces(simSnpReader);
            clear(simSnpPos);
            readUntilWhitespace(simSnpPos, simSnpReader);     // pos

            // simMethReader
            skipWhitespaces(simMethReader);
            clear(helper);
            readUntilWhitespace(helper, simMethReader);     // strand, '+' or '-' -> Use it to check?
            simMethO = helper[0];
            skipWhitespaces(simMethReader);                
            clear(helper);
            readUntilWhitespace(helper, simMethReader);     // methylation level called
            if (helper[0] != '.')
            {
                simMethLevel = helper[0];
            }

            skipLine(simMethReader);                        // skip
            clear(simMethContig);
            readUntilWhitespace(simMethContig, simMethReader);  // chrom
            skipWhitespaces(simMethReader);
            clear(simMethPos);
            readUntilWhitespace(simMethPos, simMethReader);     // pos

            // callReader
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // ref
            Dna5 ref = helper[0];
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qA
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qC
            unsigned countC_F = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qG
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qT
            unsigned countT_F = length(helper)-2;
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qA  (T on reverse)
            unsigned countT_R = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qC
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // qG  (C on reverse)
            unsigned countC_R = length(helper)-2;
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);        // qT
            skipWhitespaces(callReader);
            clear(callCoverage);
            readUntilWhitespace(callCoverage, callReader);        // cov
            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);        // genotype called 
            bool snpCalled = false;
            if (helper != '.')
            {
                callA1 = helper[0];
                callA2 = helper[1];
                ++c_CalledSnps;
                snpCalled = true;
                // compare called genotypes
                if (simA1 != callA1 || simA2 != callA2)
                {
                    ++c_WrongCalled;
                    if (options.verbosity >= 4)
                    {
                        std::cout << "Wrong called at position: " << callPos << std::endl;
                        std::cout << "simA1: " << simA1 << " simA2: " << simA2 << " callA1: " << callA1 << " callA2: " << callA2 << std::endl;
                    }
                }
                else
                {
                    ++c_RightCalled;
                   // std::cout << "Right called at callPos: " << callPos << std::endl;
                }
            }
            else
            {
                ++c_NotCalled;
                if (options.verbosity >= 4)
                {
                    std::cout << "Not called at position: " << callPos << " coverage: " << callCoverage << std::endl;
                }
            }

            skipWhitespaces(callReader);
            clear(helper);
            readUntilWhitespace(helper, callReader);    // meth level called
            if (helper[0] != '.')
            {
                if (!snpCalled)
                {
                    callA1 = ref;
                    callA2 = ref;
                }
                // if genotype CG -> methLevel1:methLevel2
                if (callA1 == 'C' && callA2 == 'G')
                {
                    unsigned i = 0;
                    while (i < length(helper) && helper[i] != ':') 
                        ++i;

                    CharString pre = prefix(helper, i);
                    callMethLevel = lexicalCast<double>(pre);

                    if (countC_F + countT_F > 0)
                    {
                        double guess = (double)countC_F / (double)(countC_F + countT_F);
                        // Check methlevel roughly and count if too different               // We don't know whats right and wrong, just to check if weird things happen!!!
                        if (abs(guess - callMethLevel) > 0.1)
                        {
                            ++c_falseMeths;
                            std::cout << "      Weird meth called at pos: " << callPos << std::endl;
                        }
                        if (abs(guess - callMethLevel) < 0.1)
                            ++c_rightMeths;
                    }
                    CharString suff = suffix(helper, i+1);
                    callMethLevel = lexicalCast<double>(suff); 

                    if (countC_R + countT_R > 0)
                    {
                        double guess = (double)countC_R / (double)(countC_R + countT_R);
                        // Check methlevel roughly and count if too different
                        if (abs(guess - callMethLevel) > methThres)
                        {
                            ++c_falseMeths;
                            std::cout << "      Weird meth called at pos: " << callPos << std::endl;
                        }
                       if (abs(guess - callMethLevel) < methThres)
                            ++c_rightMeths;
                    }
                }
                else
                {   
                    if (callA1 == 'C' || callA2 == 'C')
                    {
                        callMethLevel = lexicalCast<double>(helper); 
                        if (countC_F + countT_F > 0)
                        {
                            double guess = (double)countC_F / (double)(countC_F + countT_F);
                            // Check methlevel roughly and count if too different               // We don't know whats right and wrong, just to check if weird things happen!!!
                            if (abs(guess - callMethLevel) > methThres)
                            {
                                ++c_falseMeths;
                                std::cout << "      Weird meth called at pos: " << callPos << std::endl;
                            }
                           if (abs(guess - callMethLevel) < methThres)
                                ++c_rightMeths;
                        }
                    }
                    else if (callA1 == 'G' || callA2 == 'G')
                    {
                        callMethLevel = lexicalCast<double>(helper); 
                        if (countC_R + countT_R > 0)
                        {
                            double guess = (double)countC_R / (double)(countC_R + countT_R);
                            // Check methlevel roughly and count if too different
                            if (abs(guess - callMethLevel) > methThres)
                            {
                                ++c_falseMeths;
                                std::cout << "      Weird meth called at pos: " << callPos << std::endl;
                            }

                            if (abs(guess - callMethLevel) < methThres)
                                ++c_rightMeths;
                        }
                    }
                }
            }
            skipWhitespaces(callReader);
            readUntilOneOf(helper, callReader, ':');    // genotype prob
            skipWhitespaces(callReader);
            readUntilWhitespace(helper, callReader);    // genotype score

            skipLine(callReader);
            clear(callContig);
            readUntilWhitespace(callContig, callReader);    // chrom
            skipWhitespaces(callReader);
            clear(callPos);
            readUntilWhitespace(callPos, callReader);       // pos

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
    out<< "Positions with fuzzy C-T counts (no snp positions): \t" << c_Fuzzy << '\n';
    out<< "False positive at fuzzy sim. level: \t" << c_FPFuzzy << '\n';


    out<< "Wrong called genotype: \t" << c_WrongCalled << '\n'; 
    out<< "Right called genotype: \t" << c_RightCalled << '\n';
    


    // Sensitivity
    out<< '\n';
    out<< "Sensitivity: " << (double)c_RightCalled/(double)(c_RightCalled + c_FalsePositive) << '\n';

    // Specificity
    // No. of true negatives: TODO wrong called not taken into account?
    unsigned trueNegatives = maxPos - c_FalsePositive;  
    out<< "Specificity: "  << (double)trueNegatives/(double)(trueNegatives + c_FalsePositive) <<  '\n';



    // Methylation calls 
    out<< '\n';
    out<< "Methylation level calls: " <<  '\n';
    out<< "  Called like expected: " << c_rightMeths << '\n';
    out<< "  Weird called: "  << c_falseMeths << '\n';
    out<< "Weird snp meths called: "  << c_weirdSnpMeths << '\n';



    return 0;
}








    

 
#endif

